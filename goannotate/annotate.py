"""Per-gene annotation: clustering, consensus name/symbol selection, GO term collection."""

import logging
import math
import multiprocessing
import os
import re
import statistics
import time
from collections import Counter, defaultdict
from difflib import SequenceMatcher

from goannotate.blast import merge_hits_by_gene
from goannotate.filters import BadNameFilter
from goannotate.go import GOHierarchy, ELECTRONIC_CODES, GO_ROOT_TERMS
from goannotate.mappings import SymbolResolver, collect_go_evidence, get_go_terms_for_hits
from goannotate.models import QueryAnnotation

logger = logging.getLogger(__name__)

# Hit tuples handled here: (BlastHit, cleaned_name, symbol, symbol_origin, source_transcript)
HIT, NAME, SYMBOL, ORIGIN, SOURCE_TX = range(5)

_TOKEN = re.compile(r'\w+')


def format_symbol(symbol: str, case: str) -> str:
    """Apply the output casing convention: upper, title (Avp), lower, or source (as found)."""
    if not symbol:
        return ""
    if case == "upper":
        return symbol.upper()
    if case == "lower":
        return symbol.lower()
    if case == "title":
        return symbol[0].upper() + symbol[1:].lower()
    return symbol


############################
# CLUSTERING AND SIMILARITY
############################

def normalize_gene_name(name: str) -> str:
    """
    Normalize a gene name for comparison: lowercase, collapse whitespace,
    and drop the words the bad-names file marks as comparison noise
    (normalize_regex, e.g. "homolog"). Output names are not affected.
    """
    if not name:
        return ""
    if _worker_bad_name_filter is not None:
        return _worker_bad_name_filter.normalize_for_comparison(name)
    return ' '.join(name.lower().split())


def _best_name_index(indices: list, hits_with_names: list, bad_name_filter: BadNameFilter,
                     min_fraction: float = 0.0) -> int:
    """
    Index of the hit whose name becomes the product: the best-scoring hit
    whose name is not family-level (class C), falling back to the best hit
    when every name is family-level. With min_fraction > 0 a specific name
    is only preferred over a family-level one when its bitscore is at least
    that fraction of the best bitscore.
    """
    ranked = sorted(indices, key=lambda i: hits_with_names[i][HIT].bitscore, reverse=True)
    best = ranked[0]
    if bad_name_filter is None or not bad_name_filter.is_generic_name(hits_with_names[best][NAME]):
        return best
    floor = min_fraction * hits_with_names[best][HIT].bitscore
    for i in ranked[1:]:
        if hits_with_names[i][HIT].bitscore < floor:
            break
        if not bad_name_filter.is_generic_name(hits_with_names[i][NAME]):
            return i
    return best


####################
# WORD WEIGHTS
####################

# Words that carry no information about what a protein is (after the AHRD
# blacklist and the PANNZER2 uninformative-word list). They are ignored when
# names are compared; the shipped list is goannotate/data/stop_words.txt.
_DEFAULT_STOP_WORDS = frozenset("""
protein proteins putative predicted probable possible hypothetical uncharacterized uncharacterised
conserved homolog homologue homologous similar like family member fragment partial isoform variant
precursor domain containing type product gene cdna clone orf unnamed unknown expressed sequence complete
cds genome genomic shotgun chromosome contig novel related low quality of and the with in to a an or for
from by at on
""".split())


def default_stop_words_path():
    from pathlib import Path
    return Path(__file__).resolve().parent / "data" / "stop_words.txt"


def load_stop_words(path: str = None) -> frozenset:
    """Stop words from a file (one per line, # comments): the given one, else the shipped list."""
    if not path:
        shipped = default_stop_words_path()
        if not shipped.exists():
            return _DEFAULT_STOP_WORDS
        path = str(shipped)
    words = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                words.add(line.lower())
    return frozenset(words)


class WordWeights:
    """
    Weights of name tokens for the weighted Jaccard similarity.

    scheme "idf": every stop word weighs 0; every other token weighs
    log((N + 1) / (df + 1)), where N is the number of distinct normalized
    names in the run and df the number of names containing the token, so a
    word that appears in most names (receptor, kinase, subunit) counts for
    less than a designator that appears in few (4K36, FOXO3). A token not
    seen in the run (reference protein names used later) takes the largest
    weight.
    scheme "stoplist": stop words weigh 0, all other tokens 1.
    """

    def __init__(self, scheme: str = "idf", stop_words: frozenset = None):
        self.scheme = scheme
        self.stop_words = stop_words if stop_words is not None else _DEFAULT_STOP_WORDS
        self.weights = {}
        self.default = 1.0
        self.n_names = 0

    def fit(self, normalized_names) -> "WordWeights":
        if self.scheme != "idf":
            return self
        df = Counter()
        n = 0
        for name in set(normalized_names):
            n += 1
            for tok in set(_TOKEN.findall(name)):
                df[tok] += 1
        self.n_names = n
        if n:
            self.weights = {tok: math.log((n + 1) / (c + 1)) for tok, c in df.items()}
            self.default = math.log((n + 1) / 2)
        return self

    def weight(self, token: str) -> float:
        if token in self.stop_words:
            return 0.0
        if self.scheme == "idf":
            return self.weights.get(token, self.default)
        return 1.0

    def summary(self, top: int = 15) -> str:
        if self.scheme != "idf" or not self.weights:
            return f"scheme {self.scheme}, {len(self.stop_words)} stop words"
        light = sorted(((w, t) for t, w in self.weights.items() if t not in self.stop_words))[:top]
        return (f"scheme idf over {self.n_names:,} distinct names, {len(self.stop_words)} stop words; "
                f"lightest non-stop words: " + ", ".join(f"{t} {w:.2f}" for w, t in light))


_worker_word_weights = None


def name_similarity(norm1: str, norm2: str) -> float:
    """
    Name-only similarity between two normalized protein names:
    max(Jaccard token similarity, SequenceMatcher ratio). With word weights
    set (annotate_queries), the Jaccard similarity is weighted: stop words
    count 0 and common words count less than rare ones; when every shared
    token is a stop word the unweighted value is used.
    """
    tokens1 = set(_TOKEN.findall(norm1))
    tokens2 = set(_TOKEN.findall(norm2))
    if tokens1 and tokens2:
        if _worker_word_weights is not None:
            w = _worker_word_weights.weight
            union = sum(w(t) for t in tokens1 | tokens2)
            if union > 0:
                jaccard = sum(w(t) for t in tokens1 & tokens2) / union
            else:
                jaccard = len(tokens1 & tokens2) / len(tokens1 | tokens2)
        else:
            jaccard = len(tokens1 & tokens2) / len(tokens1 | tokens2)
    else:
        jaccard = 0.0

    seq_ratio = SequenceMatcher(None, norm1, norm2).ratio()
    return max(jaccard, seq_ratio)


def pair_similarity(name1: str, symbol1: str, name2: str, symbol2: str) -> float:
    """
    Similarity between two hits using both protein name and gene symbol.
    When both hits have symbols, symbol agreement boosts similarity and
    symbol disagreement caps it (prevents paralog merging).
    """
    norm1 = normalize_gene_name(name1)
    norm2 = normalize_gene_name(name2)
    name_sim = name_similarity(norm1, norm2)

    if symbol1 and symbol2:
        sym_sim = SequenceMatcher(None, symbol1.lower(), symbol2.lower()).ratio()
        if sym_sim == 1.0:
            return min(1.0, name_sim + 0.2)
        if sym_sim < 0.8:
            return min(name_sim, 0.3)
        return name_sim

    return name_sim


def cluster_hits_agglomerative(hits_with_names: list, similarity_threshold: float = 0.5) -> list:
    """
    Average-linkage agglomerative clustering on hits.

    hits_with_names: list of hit tuples (name at index 1, symbol at index 2)
    Returns a list of clusters, each a list of indices into hits_with_names.

    The summed similarity between every pair of clusters is kept in a
    table and updated when two clusters merge (the sums of the merged
    cluster are the sums of its two parts), so a merge costs one pass over
    the clusters instead of a pass over all pairs of hits. The average
    linkage, the merge order and the stopping rule are those of the
    original implementation.
    """
    n = len(hits_with_names)
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    total = [[0.0] * n for _ in range(n)]   # summed similarity between clusters, by representative
    for i in range(n):
        name_i, sym_i = hits_with_names[i][1], hits_with_names[i][2]
        row = total[i]
        for j in range(i + 1, n):
            v = pair_similarity(name_i, sym_i, hits_with_names[j][1], hits_with_names[j][2])
            row[j] = v
            total[j][i] = v

    members = {i: [i] for i in range(n)}
    active = list(range(n))      # cluster representatives, in the order of the original dict

    while len(active) > 1:
        best_sim = -1.0
        best_pair = None
        for ia in range(len(active)):
            a = active[ia]
            row = total[a]
            size_a = len(members[a])
            for ib in range(ia + 1, len(active)):
                b = active[ib]
                avg = row[b] / (size_a * len(members[b]))
                if avg > best_sim:
                    best_sim = avg
                    best_pair = (a, b)

        if best_sim < similarity_threshold:
            break

        a, b = best_pair
        row_a, row_b = total[a], total[b]
        for c in active:
            if c != a and c != b:
                v = row_a[c] + row_b[c]
                row_a[c] = v
                total[c][a] = v
        members[a] = members[a] + members[b]
        del members[b]
        active.remove(b)

    return [members[a] for a in active]


def collapse_gene_hits(hits: list, max_hits: int = 300) -> tuple:
    """
    One entry per UniProt subject for a gene: isoforms of a gene hit largely
    the same proteins, so the pooled hits of all isoforms repeat each subject
    once per isoform. The best-scoring entry is kept, with the tuple of all
    isoforms that hit the subject as its source; the result is ordered by
    bitscore and capped at max_hits (0 for no cap).

    Returns (collapsed hits, number of hits before collapsing).
    """
    best = {}
    sources = {}
    for t in hits:
        sid = t[HIT].subject_id
        src = t[SOURCE_TX]
        sources.setdefault(sid, set()).update(src if isinstance(src, tuple) else (src,))
        if sid not in best or t[HIT].bitscore > best[sid][HIT].bitscore:
            best[sid] = t
    out = [t[:SOURCE_TX] + (tuple(sorted(sources[sid])),) for sid, t in best.items()]
    out.sort(key=lambda t: -t[HIT].bitscore)
    if max_hits and len(out) > max_hits:
        out = out[:max_hits]
    return out, len(hits)


def _isoforms(t) -> tuple:
    src = t[SOURCE_TX]
    return src if isinstance(src, tuple) else (src,)


################################
# CLUSTER SCORING AND SELECTION
################################

def select_winning_cluster(clusters: list, hits_with_names: list) -> tuple:
    """
    Score clusters and select the winner. Always returns the best-scoring
    cluster; there is no fallback to all hits.

    Score = sum(bitscores) * sqrt(cluster_size) * (1 + 0.1 * (n_isoforms - 1))

    Returns (winning_cluster_indices, cluster_fraction).
    """
    if not clusters:
        return [], 0.0

    total_hits = len(hits_with_names)
    best_score = -1.0
    best_cluster = None

    for cluster in clusters:
        bitscores_sum = sum(hits_with_names[i][HIT].bitscore for i in cluster)
        size = len(cluster)
        n_isoforms = len(set(tx for i in cluster for tx in _isoforms(hits_with_names[i])))
        score = bitscores_sum * math.sqrt(size) * (1 + 0.1 * (n_isoforms - 1))
        if score > best_score:
            best_score = score
            best_cluster = cluster

    cluster_fraction = len(best_cluster) / total_hits if total_hits > 0 else 0.0
    return best_cluster, cluster_fraction


def symbol_group_score(bitscores: list, method: str) -> float:
    """
    Score of a symbol group: the mean bitscore of its hits ("mean"), or the
    mean times the square root of the group size ("mean_sqrt"), which lets
    agreement among many species count once symbols are in one nomenclature.
    """
    mean = sum(bitscores) / len(bitscores)
    if method == "mean_sqrt":
        return mean * math.sqrt(len(bitscores))
    return mean


def select_name_and_symbol(winning_indices: list,
                           hits_with_names: list,
                           bad_name_filter: BadNameFilter,
                           resolver: SymbolResolver = None,
                           score_method: str = "mean") -> tuple:
    """
    Select the protein name and gene symbol from the winning cluster.

    Hits are grouped by gene symbol. The highest-scoring symbol group
    (see symbol_group_score) determines both the symbol and the protein
    name (from that group's best hit). A single-hit group can only win when
    no other group has 2+ hits, and a lone symbol representing <10% of a
    cluster of 10+ hits is rejected.

    Returns (protein_name, gene_symbol, symbol_origin, winning_symbol_indices).
    """
    if not winning_indices:
        return "", "", "", []

    symbol_groups = defaultdict(list)
    for idx in winning_indices:
        symbol = hits_with_names[idx][SYMBOL]
        key = symbol.lower() if symbol else ""
        symbol_groups[key].append(idx)

    symbol_groups.pop("", [])

    # With ortholog resolution, a symbol is "in nomenclature" when it was
    # resolved to the reference species, comes from a reference-species
    # hit, or is itself a reference-species symbol in any case (zebrafish
    # il15 is IL15). Groups in nomenclature are tried first; when they give
    # no symbol (a lone resolved hit among many, rejected by the guard
    # below), all groups are tried, so that a well-supported symbol of
    # another species (slc6a4b on 45 of 104 hits) is not lost to a single
    # hit that happens to resolve.
    rounds = [symbol_groups]
    if symbol_groups and resolver is not None and resolver.has_orthologs:
        def in_nomenclature(i):
            t = hits_with_names[i]
            return (t[ORIGIN] == "ortholog"
                    or (t[ORIGIN] == "uniprot" and t[HIT].taxid == resolver.reference_taxid)
                    or resolver.is_reference_symbol(t[SYMBOL]))
        resolved_groups = {s: idxs for s, idxs in symbol_groups.items()
                           if any(in_nomenclature(i) for i in idxs)}
        if resolved_groups and len(resolved_groups) < len(symbol_groups):
            rounds = [resolved_groups, symbol_groups]

    for groups in rounds:
        if not groups:
            continue
        multi_hit = {s: idxs for s, idxs in groups.items() if len(idxs) >= 2}
        single_hit = {s: idxs for s, idxs in groups.items() if len(idxs) == 1}
        candidates = multi_hit if multi_hit else single_hit

        best_symbol_key = None
        best_score = -1.0
        for sym, idxs in candidates.items():
            score = symbol_group_score([hits_with_names[i][HIT].bitscore for i in idxs], score_method)
            if score > best_score:
                best_score = score
                best_symbol_key = sym

        winning_sym_indices = candidates[best_symbol_key]

        if not (len(winning_sym_indices) == 1
                and len(winning_indices) >= 10
                and len(winning_sym_indices) / len(winning_indices) < 0.1):
            best_idx = max(winning_sym_indices, key=lambda i: hits_with_names[i][HIT].bitscore)
            # Within one symbol group every hit names the same gene, so the
            # most specific name in the group is the product ("Olfactory
            # receptor 4K36" over "Olfactory receptor"), whatever its rank.
            name_idx = _best_name_index(winning_sym_indices, hits_with_names, bad_name_filter)
            protein_name = hits_with_names[name_idx][NAME]
            gene_symbol = hits_with_names[best_idx][SYMBOL]
            origin = hits_with_names[best_idx][ORIGIN]
            return protein_name, gene_symbol, origin, winning_sym_indices

    # No-symbol path (also reached when the single-symbol guard fires above).
    # Without a symbol tying the hits together, a specific name replaces a
    # family-level one only when it scores at least 80% of the best hit.
    best_idx = max(winning_indices, key=lambda i: hits_with_names[i][HIT].bitscore)
    name_idx = _best_name_index(winning_indices, hits_with_names, bad_name_filter, min_fraction=0.8)
    protein_name = hits_with_names[name_idx][NAME]
    gene_symbol = ""
    origin = ""

    if resolver is not None and resolver.has_ncbi:
        ncbi_symbol = resolver.ncbi_symbol(hits_with_names[best_idx][HIT].accession)
        if ncbi_symbol and not bad_name_filter.is_bad_symbol(ncbi_symbol):
            gene_symbol = ncbi_symbol
            origin = "ncbi"

    return protein_name, gene_symbol, origin, list(winning_indices)


####################
# CONFIDENCE
####################

def compute_confidence(annot: QueryAnnotation) -> tuple:
    """
    Per-gene confidence in the assigned annotation.

    Score in [0, 1] combining:
      0.45 * name concordance across all merged hits
      0.25 * min(1, winning cluster size / 10)
      0.15 * symbol present
      0.15 * symbol concordance across all merged hits
    minus 0.1 when no hit of the gene passed the hit coverage filter and
    the annotation rests on partial alignments (coverage status fallback).

    Labels: none (no hits), low (<0.4), medium (<0.7), high (>=0.7).
    Returns (label, score).
    """
    if not annot.hits:
        return "none", 0.0
    score = (0.45 * annot.name_concordance
             + 0.25 * min(1.0, annot.n_cluster_hits / 10.0)
             + 0.15 * (1.0 if annot.consensus_symbol else 0.0)
             + 0.15 * annot.symbol_concordance)
    if annot.coverage_status == "fallback":
        score = max(0.0, score - 0.1)
    if score >= 0.7:
        label = "high"
    elif score >= 0.4:
        label = "medium"
    else:
        label = "low"
    return label, round(score, 3)


###########################
# PER-GENE WORKER
###########################

_worker_go_mapping = None
_worker_go_hierarchy = None
_worker_bad_name_filter = None
_worker_resolver = None
_worker_forbidden_go = frozenset()
_worker_go_iea_min_hits = 1


def _init_worker(go_mapping, go_hierarchy, bad_name_filter, resolver, word_weights=None,
                 forbidden_go=None, go_iea_min_hits=1):
    """Initialize shared read-only data in each worker process."""
    global _worker_go_mapping, _worker_go_hierarchy, _worker_bad_name_filter, _worker_resolver
    global _worker_word_weights, _worker_forbidden_go, _worker_go_iea_min_hits
    _worker_go_mapping = go_mapping
    _worker_go_hierarchy = go_hierarchy
    _worker_bad_name_filter = bad_name_filter
    _worker_resolver = resolver
    _worker_word_weights = word_weights
    _worker_forbidden_go = frozenset(forbidden_go or ())
    _worker_go_iea_min_hits = go_iea_min_hits


def _apply_go_evidence_rules(term_codes: dict, term_support: dict, iea_min_hits: int) -> set:
    """
    Terms kept under the evidence rules: a term whose only evidence is
    electronic (IEA) must be carried by at least iea_min_hits of the source
    hits; terms with any other evidence, or with no recorded codes (older
    databases), need one hit. Root terms are never assigned.
    """
    kept = set()
    for term, codes in term_codes.items():
        if term in GO_ROOT_TERMS:
            continue
        if iea_min_hits > 1 and codes and codes <= ELECTRONIC_CODES and term_support.get(term, 0) < iea_min_hits:
            continue
        kept.add(term)
    return kept


def annotate_single_gene(gene_id: str, hits_with_names: list,
                         consensus_threshold: float, namespaces: list,
                         symbol_case: str = "upper", symbol_score: str = "mean") -> tuple:
    """
    Full annotation of a single gene (runs in a worker process).

    Returns (gene_id, result_dict).
    """
    try:
        go_mapping = _worker_go_mapping
        go_hierarchy = _worker_go_hierarchy
        bad_name_filter = _worker_bad_name_filter
        resolver = _worker_resolver

        result = {
            'hits': [],
            'filtered_hits': [],
            'consensus_name': '',
            'consensus_symbol': '',
            'go_terms': set(),
            'specific_go_terms': set(),
            'name_concordance': 0.0,
            'symbol_concordance': 0.0,
            'cluster_fraction': 0.0,
            'had_hits': False,
            'had_protein_name': False,
            'had_symbol': False,
            'symbol_source': '',
            'n_isoforms': 0,
            'n_cluster_hits': 0,
            'winning_cluster_names': [],
            'winning_cluster_symbols': [],
            'mean_cluster_bitscore': 0.0,
            'go_source': '',
            'go_evidence': {},
            'go_support': {},
            'go_taxon_removed': set(),
            'best_query_coverage': 0.0,
            'best_subject_coverage': 0.0,
        }

        if not hits_with_names:
            return (gene_id, result)

        result['had_hits'] = True
        result['n_isoforms'] = len(set(tx for t in hits_with_names for tx in _isoforms(t)))
        result['hits'] = sorted([t[HIT] for t in hits_with_names],
                                key=lambda h: h.bitscore, reverse=True)

        # Cluster and select the winning cluster
        clusters = cluster_hits_agglomerative(hits_with_names,
                                              similarity_threshold=consensus_threshold)
        winning_indices, cluster_fraction = select_winning_cluster(clusters, hits_with_names)
        result['cluster_fraction'] = cluster_fraction
        result['n_cluster_hits'] = len(winning_indices)
        result['winning_cluster_names'] = [hits_with_names[i][NAME] for i in winning_indices]
        result['winning_cluster_symbols'] = [hits_with_names[i][SYMBOL] for i in winning_indices
                                             if hits_with_names[i][SYMBOL]]
        if winning_indices:
            result['mean_cluster_bitscore'] = (
                sum(hits_with_names[i][HIT].bitscore for i in winning_indices) / len(winning_indices))

        # Select name and symbol
        protein_name, gene_symbol, symbol_origin, winning_symbol_indices = select_name_and_symbol(
            winning_indices, hits_with_names, bad_name_filter, resolver=resolver,
            score_method=symbol_score)

        result['consensus_name'] = protein_name
        result['consensus_symbol'] = format_symbol(gene_symbol, symbol_case)

        if protein_name:
            result['had_protein_name'] = True
        if gene_symbol:
            result['had_symbol'] = True
            result['symbol_source'] = symbol_origin

        # Concordance against the full merged hit pool
        if protein_name and hits_with_names:
            final_name_norm = normalize_gene_name(protein_name)
            name_match_count = sum(
                1 for t in hits_with_names
                if name_similarity(final_name_norm, normalize_gene_name(t[NAME])) >= consensus_threshold)
            result['name_concordance'] = name_match_count / len(hits_with_names)

        if gene_symbol and hits_with_names:
            final_sym_lower = gene_symbol.lower()
            symbol_match_count = sum(1 for t in hits_with_names
                                     if t[SYMBOL] and t[SYMBOL].lower() == final_sym_lower)
            result['symbol_concordance'] = symbol_match_count / len(hits_with_names)

        result['filtered_hits'] = sorted([hits_with_names[i][HIT] for i in winning_symbol_indices],
                                         key=lambda h: h.bitscore, reverse=True)
        if result['filtered_hits']:
            result['best_query_coverage'] = max(h.query_coverage for h in result['filtered_hits'])
            result['best_subject_coverage'] = max(h.subject_coverage for h in result['filtered_hits'])

        # GO terms from the winning symbol group (top 5 by bitscore), with
        # the evidence codes behind each term and the number of hits that
        # carry it
        go_source_indices = sorted(winning_symbol_indices,
                                   key=lambda i: hits_with_names[i][HIT].bitscore, reverse=True)[:5]
        go_hits = [hits_with_names[i][HIT] for i in go_source_indices]
        term_codes, term_support = collect_go_evidence(go_hits, go_mapping)
        go_source = "cluster_top5" if term_codes else ""

        # SwissProt GO supplement
        winning_name_norm = normalize_gene_name(protein_name) if protein_name else ""
        winning_sym_lower = gene_symbol.lower() if gene_symbol else ""
        go_source_set = set(go_source_indices)

        for idx in winning_indices:
            if idx in go_source_set:
                continue
            hit, cleaned_name, symbol = (hits_with_names[idx][HIT], hits_with_names[idx][NAME],
                                         hits_with_names[idx][SYMBOL])
            if not hit.is_swissprot:
                continue
            sym_match = (winning_sym_lower and symbol and symbol.lower() == winning_sym_lower)
            name_match = False
            if not sym_match and winning_name_norm and cleaned_name:
                name_norm = normalize_gene_name(cleaned_name)
                sim = pair_similarity(winning_name_norm, "", name_norm, "")
                name_match = sim >= consensus_threshold
            if sym_match or name_match:
                sp_codes, sp_support = collect_go_evidence([hit], go_mapping)
                if sp_codes:
                    if any(t not in term_codes for t in sp_codes):
                        go_source = "cluster_top5+swissprot_supplement"
                    for t, codes in sp_codes.items():
                        term_codes.setdefault(t, set()).update(codes)
                        term_support[t] = term_support.get(t, 0) + sp_support[t]

        result['go_source'] = go_source

        go_terms = _apply_go_evidence_rules(term_codes, term_support, _worker_go_iea_min_hits)

        if namespaces and go_terms:
            go_terms = go_hierarchy.filter_by_namespace(go_terms, namespaces)

        if go_terms and _worker_forbidden_go:
            removed = go_terms & _worker_forbidden_go
            if removed:
                go_terms = go_terms - removed
                result['go_taxon_removed'] = removed

        result['go_terms'] = go_terms
        if go_terms:
            result['specific_go_terms'] = go_hierarchy.filter_to_specific(go_terms)
            result['go_evidence'] = {t: ",".join(sorted(term_codes[t])) for t in go_terms}
            result['go_support'] = {t: term_support.get(t, 0) for t in go_terms}

        return (gene_id, result)

    except Exception as e:
        raise RuntimeError(f"Error annotating gene {gene_id}: {e}") from e


###########################
# PIPELINE
###########################

def _annotate_task(task):
    return annotate_single_gene(*task)


def annotate_queries(resolved_transcript_hits: dict,
                     transcript_to_gene: dict,
                     go_mapping: dict,
                     go_hierarchy: GOHierarchy,
                     bad_name_filter: BadNameFilter,
                     all_gene_ids: list,
                     consensus_threshold: float = 0.5,
                     namespaces: list = None,
                     resolver: SymbolResolver = None,
                     threads: int = None,
                     symbol_case: str = "upper",
                     coding_gene_ids: set = None,
                     symbol_score: str = "mean",
                     max_gene_hits: int = 300,
                     word_weights: WordWeights = None,
                     forbidden_go: set = None,
                     go_iea_min_hits: int = 1,
                     coverage_status: dict = None) -> dict:
    """
    Main annotation pipeline: merge isoform hits, cluster, select name/symbol,
    collect GO terms. Per-gene annotation runs in parallel via multiprocessing,
    followed by symbol propagation and confidence assignment.

    coding_gene_ids, when given, is the set of genes that can have hits at
    all (protein-coding genes); coverage statistics are reported against
    it. Without it every gene is treated as coding.
    word_weights: token weights for the name similarity (None: unweighted).
    forbidden_go: GO terms that the taxon constraints exclude for the query
    organism; they are removed from every gene.
    go_iea_min_hits: hits that must carry a term whose only evidence is
    electronic (IEA) before it is assigned (1: no requirement).
    coverage_status: transcript -> "pass" / "fallback" from the hit coverage
    filter; a gene is "pass" when any of its transcripts passed.

    Returns: dict of gene_id -> QueryAnnotation
    """
    gene_hits = merge_hits_by_gene(resolved_transcript_hits, transcript_to_gene)
    # Collapse repeated subjects across isoforms and cap the pool: clustering
    # compares all pairs of a gene's hits, so a gene with 147 isoforms and 50
    # hits each would otherwise cluster 7,350 hits.
    n_before = n_after = n_capped = 0
    largest = 0
    for gid in list(gene_hits):
        collapsed, before = collapse_gene_hits(gene_hits[gid], max_gene_hits)
        unique = len({t[HIT].subject_id for t in gene_hits[gid]})
        if max_gene_hits and unique > max_gene_hits:
            n_capped += 1
        largest = max(largest, before)
        n_before += before
        n_after += len(collapsed)
        gene_hits[gid] = collapsed
    if gene_hits:
        logger.info(f"Hits pooled per gene across isoforms: {n_before:,} -> {n_after:,} after keeping each "
                    f"UniProt subject once (largest pool {largest:,} hits); "
                    f"{n_capped:,} genes capped at {max_gene_hits} unique hits")
    if resolver is None:
        resolver = SymbolResolver()

    if word_weights is not None:
        # Comparison-only normalization is applied before the weights are fitted,
        # so the tokens match those the similarity sees.
        _init_worker(go_mapping, go_hierarchy, bad_name_filter, resolver)
        word_weights.fit(normalize_gene_name(t[NAME]) for hits in gene_hits.values() for t in hits)
        logger.info(f"Word weights for name comparison: {word_weights.summary()}")
    if forbidden_go:
        logger.info(f"GO taxon constraints exclude {len(forbidden_go):,} terms for the query organism")
    if go_iea_min_hits > 1:
        logger.info(f"GO terms with electronic evidence only need {go_iea_min_hits} supporting hits")

    n_workers = threads if threads is not None else os.cpu_count()
    n_genes_with_hits = sum(1 for gid in all_gene_ids if gid in gene_hits)
    n_workers = min(n_workers, max(1, n_genes_with_hits))

    logger.info(f"Annotating {len(all_gene_ids):,} genes "
                f"({n_genes_with_hits:,} have filtered hits) "
                f"using {n_workers} workers...")

    tasks = [(gene_id, gene_hits.get(gene_id, []), consensus_threshold, namespaces, symbol_case, symbol_score)
             for gene_id in all_gene_ids]

    start_time = time.time()

    # The main process uses the same globals (normalize_gene_name reads the
    # bad-name filter for comparison-only normalization during propagation).
    init_args = (go_mapping, go_hierarchy, bad_name_filter, resolver, word_weights, forbidden_go, go_iea_min_hits)
    _init_worker(*init_args)
    if n_workers == 1:
        results = [annotate_single_gene(*task) for task in tasks]
    else:
        chunksize = max(1, len(tasks) // (n_workers * 4))
        chunksize = min(chunksize, 50)
        results = []
        next_report = time.time() + 60
        with multiprocessing.Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=init_args
        ) as pool:
            for res in pool.imap(_annotate_task, tasks, chunksize=chunksize):
                results.append(res)
                if time.time() >= next_report:
                    logger.info(f"  {len(results):,} of {len(tasks):,} genes annotated "
                                f"({time.time() - start_time:.0f}s)")
                    next_report = time.time() + 60

    elapsed = time.time() - start_time
    genes_per_sec = len(all_gene_ids) / elapsed if elapsed > 0 else 0
    logger.info(f"  Per-gene annotation completed in {elapsed:.1f}s ({genes_per_sec:.0f} genes/sec)")

    annotations = {}
    cluster_fractions = []
    stats = {
        'total_genes': len(all_gene_ids),
        'coding_genes': len(coding_gene_ids) if coding_gene_ids is not None else len(all_gene_ids),
        'genes_with_hits': 0,
        'genes_with_protein_name': 0,
        'genes_with_go_terms': 0,
        'symbols_from_uniprot': 0,
        'symbols_from_ortholog': 0,
        'symbols_from_ncbi': 0,
        'genes_with_taxon_removed': 0,
        'terms_taxon_removed': 0,
    }

    for gene_id, result_dict in results:
        annot = QueryAnnotation(query_id=gene_id)
        annot.hits = result_dict['hits']
        annot.filtered_hits = result_dict['filtered_hits']
        annot.consensus_name = result_dict['consensus_name']
        annot.consensus_symbol = result_dict['consensus_symbol']
        annot.go_terms = result_dict['go_terms']
        annot.specific_go_terms = result_dict['specific_go_terms']
        annot.name_concordance = result_dict['name_concordance']
        annot.symbol_concordance = result_dict['symbol_concordance']
        annot.n_isoforms = result_dict['n_isoforms']
        annot.n_cluster_hits = result_dict['n_cluster_hits']
        annot.winning_cluster_names = result_dict['winning_cluster_names']
        annot.winning_cluster_symbols = result_dict['winning_cluster_symbols']
        annot.mean_cluster_bitscore = result_dict['mean_cluster_bitscore']
        annot.go_source = result_dict['go_source']
        annot.go_evidence = result_dict['go_evidence']
        annot.go_support = result_dict['go_support']
        annot.go_taxon_removed = result_dict['go_taxon_removed']
        annot.best_query_coverage = result_dict['best_query_coverage']
        annot.best_subject_coverage = result_dict['best_subject_coverage']
        if result_dict['had_hits']:
            annot.symbol_status = "vote" if result_dict['had_symbol'] else "none"
            if coverage_status:
                statuses = {coverage_status.get(tx, "") for t in gene_hits.get(gene_id, []) for tx in _isoforms(t)}
                if "pass" in statuses:
                    annot.coverage_status = "pass"
                elif "fallback" in statuses:
                    annot.coverage_status = "fallback"
        annotations[gene_id] = annot

        if result_dict['go_taxon_removed']:
            stats['genes_with_taxon_removed'] += 1
            stats['terms_taxon_removed'] += len(result_dict['go_taxon_removed'])
        if result_dict['had_hits']:
            stats['genes_with_hits'] += 1
        if result_dict['had_protein_name']:
            stats['genes_with_protein_name'] += 1
        if result_dict['had_symbol']:
            if result_dict['symbol_source'] == 'uniprot':
                stats['symbols_from_uniprot'] += 1
            elif result_dict['symbol_source'] == 'ortholog':
                stats['symbols_from_ortholog'] += 1
            elif result_dict['symbol_source'] == 'ncbi':
                stats['symbols_from_ncbi'] += 1
        if result_dict['go_terms']:
            stats['genes_with_go_terms'] += 1
        if result_dict['cluster_fraction'] > 0:
            cluster_fractions.append(result_dict['cluster_fraction'])

    # Symbol propagation: genes sharing an identical normalized protein name
    # receive the most common symbol among them when they lack one.
    name_to_genes = defaultdict(list)
    for gene_id, annot in annotations.items():
        if annot.consensus_name:
            name_to_genes[normalize_gene_name(annot.consensus_name)].append(gene_id)

    symbols_propagated = 0
    propagated_gene_ids = []
    for norm_name, gene_id_list in name_to_genes.items():
        if len(gene_id_list) < 2:
            continue
        existing_symbols = [annotations[gid].consensus_symbol for gid in gene_id_list
                            if annotations[gid].consensus_symbol]
        if not existing_symbols:
            continue
        genes_without = [gid for gid in gene_id_list if not annotations[gid].consensus_symbol]
        if not genes_without:
            continue
        best_symbol = Counter(existing_symbols).most_common(1)[0][0]
        best_lower = best_symbol.lower()
        for gid in genes_without:
            # Only when the gene's own hits carry the symbol: identical
            # product names alone ("disks large 5") are not evidence that
            # two genes are copies of the same gene.
            if not any(t[SYMBOL] and t[SYMBOL].lower() == best_lower for t in gene_hits.get(gid, [])):
                continue
            annotations[gid].consensus_symbol = best_symbol
            symbols_propagated += 1
            propagated_gene_ids.append(gid)

    stats['symbols_from_propagation'] = symbols_propagated
    if symbols_propagated > 0:
        logger.info(f"  Symbol propagation: {symbols_propagated:,} genes received symbols "
                    f"from genes with identical protein names")

    for gid in propagated_gene_ids:
        annot = annotations[gid]
        annot.symbol_status = "propagated"
        hits_pool = gene_hits.get(gid, [])
        if annot.consensus_symbol and hits_pool:
            final_sym_lower = annot.consensus_symbol.lower()
            match_count = sum(1 for t in hits_pool
                              if t[SYMBOL] and t[SYMBOL].lower() == final_sym_lower)
            annot.symbol_concordance = match_count / len(hits_pool)

    # Confidence
    confidence_counts = Counter()
    for annot in annotations.values():
        annot.confidence, annot.confidence_score = compute_confidence(annot)
        confidence_counts[annot.confidence] += 1

    # Statistics (coverage is reported against coding genes)
    coding = max(1, stats['coding_genes'])
    logger.info("Annotation statistics:")
    logger.info(f"  Total genes: {stats['total_genes']:,}; coding genes: {stats['coding_genes']:,}")
    logger.info(f"  Coding genes with filtered hits: {stats['genes_with_hits']:,} "
                f"({100*stats['genes_with_hits']/coding:.1f}%)")
    logger.info(f"  Coding genes with protein name: {stats['genes_with_protein_name']:,} "
                f"({100*stats['genes_with_protein_name']/coding:.1f}%)")
    total_symbols = (stats['symbols_from_uniprot'] + stats['symbols_from_ortholog']
                     + stats['symbols_from_ncbi'] + stats['symbols_from_propagation'])
    logger.info(f"  Coding genes with gene symbol: {total_symbols:,} ({100*total_symbols/coding:.1f}%)")
    if stats['symbols_from_ortholog'] > 0:
        logger.info(f"    - From reference-species ortholog: {stats['symbols_from_ortholog']:,}")
    if stats['symbols_from_uniprot'] > 0:
        logger.info(f"    - From UniProt GN= field: {stats['symbols_from_uniprot']:,}")
    if stats['symbols_from_ncbi'] > 0:
        logger.info(f"    - From NCBI cross-reference: {stats['symbols_from_ncbi']:,}")
    if stats['symbols_from_propagation'] > 0:
        logger.info(f"    - From name-based propagation: {stats['symbols_from_propagation']:,}")
    logger.info(f"  Coding genes with GO terms: {stats['genes_with_go_terms']:,} "
                f"({100*stats['genes_with_go_terms']/coding:.1f}%)")
    if stats['terms_taxon_removed']:
        logger.info(f"  GO taxon constraints removed {stats['terms_taxon_removed']:,} term assignments "
                    f"from {stats['genes_with_taxon_removed']:,} genes")
    if coverage_status:
        cov = Counter(a.coverage_status for a in annotations.values() if a.hits)
        logger.info(f"  Hit coverage: {cov.get('pass', 0):,} genes annotated from hits passing the coverage "
                    f"filter, {cov.get('fallback', 0):,} from unfiltered hits (no hit passed)")
    logger.info("  Confidence: " + ", ".join(
        f"{label}={confidence_counts.get(label, 0):,}" for label in ("high", "medium", "low", "none")))

    if cluster_fractions:
        logger.info(f"  Winning cluster fraction: mean={statistics.mean(cluster_fractions):.2f}, "
                    f"median={statistics.median(cluster_fractions):.2f}, "
                    f"min={min(cluster_fractions):.2f}, max={max(cluster_fractions):.2f}")

    name_concordances = [a.name_concordance for a in annotations.values() if a.consensus_name]
    symbol_concordances = [a.symbol_concordance for a in annotations.values() if a.consensus_symbol]
    if name_concordances:
        logger.info(f"  Name concordance: mean={statistics.mean(name_concordances):.2f}, "
                    f"median={statistics.median(name_concordances):.2f}, "
                    f"min={min(name_concordances):.2f}, max={max(name_concordances):.2f}")
    if symbol_concordances:
        logger.info(f"  Symbol concordance: mean={statistics.mean(symbol_concordances):.2f}, "
                    f"median={statistics.median(symbol_concordances):.2f}, "
                    f"min={min(symbol_concordances):.2f}, max={max(symbol_concordances):.2f}")

    return annotations
