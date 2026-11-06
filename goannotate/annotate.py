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
from goannotate.mappings import (SymbolResolver, collect_go_evidence, extract_accession, get_go_terms_for_hits,
                                 symbol_bases)
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


ORTHOLOG_VOTE_MARGIN = 0.75     # a distant ortholog symbol wins when its group scores this fraction of the best


def _symbol_group_result(idxs: list, hits_with_names: list, bad_name_filter: BadNameFilter, origin: str = None):
    """(protein_name, symbol, origin, indices) of one symbol group: the most specific name, the best hit's symbol."""
    best_idx = max(idxs, key=lambda i: hits_with_names[i][HIT].bitscore)
    name_idx = _best_name_index(idxs, hits_with_names, bad_name_filter)
    return (hits_with_names[name_idx][NAME], hits_with_names[best_idx][SYMBOL],
            origin or hits_with_names[best_idx][ORIGIN], idxs)


def select_name_and_symbol(winning_indices: list,
                           hits_with_names: list,
                           bad_name_filter: BadNameFilter,
                           resolver: SymbolResolver = None,
                           score_method: str = "mean",
                           prior: tuple = None) -> tuple:
    """
    Select the protein name and gene symbol from the winning cluster.

    Hits are grouped by gene symbol. The highest-scoring symbol group
    (see symbol_group_score) determines both the symbol and the protein
    name (from that group's best hit). A single-hit group can only win when
    no other group has 2+ hits, and a lone symbol representing <10% of a
    cluster of 10+ hits is rejected.

    prior: (symbol, close) from the orthology module. A close relatives'
    symbol decides: its group wins whenever the cluster has one (origin
    "relative"). A distant relatives' symbol is one more vote: its group
    wins when it scores at least ORTHOLOG_VOTE_MARGIN of the best group.

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
    # Reference readthrough symbols (IQCJ-SCHIP1) never name the gene, and a reference pseudogene symbol
    # (FAM90A27P) does not when only reference-species hits carry it. A pseudogene symbol carried by hits
    # of other species names their functional gene (rodent Uox, Gucy1b2, whose human orthologs are
    # pseudogenes) and stays. The hits stay in the cluster and its other symbol groups decide.
    if resolver is not None and symbol_groups and hasattr(resolver, "is_excluded"):
        ref_taxid = str(getattr(resolver, "reference_taxid", ""))
        for key in [k for k in symbol_groups if resolver.is_excluded(k)]:
            readthrough = key.lower() not in getattr(resolver, "reference_pseudogenes", ())
            other_species = any(str(hits_with_names[i][HIT].taxid) != ref_taxid for i in symbol_groups[key])
            if readthrough or not other_species:
                symbol_groups.pop(key)

    if prior and prior[0] and symbol_groups:
        prior_key = prior[0].lower()
        if prior_key in symbol_groups:
            idxs = symbol_groups[prior_key]
            if prior[1]:
                return _symbol_group_result(idxs, hits_with_names, bad_name_filter, "relative")
            scores = {s: symbol_group_score([hits_with_names[i][HIT].bitscore for i in g], score_method)
                      for s, g in symbol_groups.items()}
            if scores[prior_key] >= ORTHOLOG_VOTE_MARGIN * max(scores.values()):
                return _symbol_group_result(idxs, hits_with_names, bad_name_filter)

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

# Orthology evidence in the confidence score (orthology module decisions)
ORTHOLOG_CONFIDENCE = {
    "close_confirmed": ("add", 0.15),     # close relatives' one-to-one orthologs name the same gene
    "distant_confirmed": ("add", 0.05),
    "close_taken": ("range", 0.55, 0.69),  # symbol from close orthologs against the hits: medium, review
    "close_filled": ("floor", 0.55),       # no symbol from the hits; close orthologs name the gene
    "distant_filled": ("floor", 0.4),
    "conflict": ("add", -0.15),            # one close relative disagrees with a reciprocal best hit
}


# Confidence subtracted for each end missing from a translated transcript's protein
FRAGMENT_PENALTY = 0.05


def compute_confidence(annot: QueryAnnotation, ortholog: str = None) -> tuple:
    """
    Per-gene confidence in the assigned annotation.

    Score in [0, 1] combining:
      0.45 * name concordance across all merged hits
      0.25 * min(1, winning cluster size / 10)
      0.15 * symbol present
      0.15 * symbol concordance across all merged hits
    minus 0.1 when no hit of the gene passed the hit coverage filter and
    the annotation rests on partial alignments (coverage status fallback),
    and, for translated transcripts, FRAGMENT_PENALTY for each end missing
    from the protein (annot.completeness: 5_partial or 3_partial one end,
    partial both).

    ortholog: the orthology module's evidence for the gene's symbol
    (ORTHOLOG_CONFIDENCE; annot.ortholog_evidence when not given): close-relative confirmation adds 0.15, distant
    0.05; a symbol taken from close relatives against the hits is placed in
    0.55-0.69 (medium, for review); a filled symbol is at least 0.55 (close)
    or 0.4 (distant), also for genes without hits; a kept conflict
    subtracts 0.15.

    Labels: none (no hits and no ortholog evidence), low (<0.4), medium
    (<0.7), high (>=0.7). Returns (label, score).
    """
    if ortholog is None:
        ortholog = annot.ortholog_evidence or None
    if not annot.hits and ortholog not in ("close_filled", "distant_filled", "close_taken"):
        return "none", 0.0
    score = (0.45 * annot.name_concordance
             + 0.25 * min(1.0, annot.n_cluster_hits / 10.0)
             + 0.15 * (1.0 if annot.consensus_symbol else 0.0)
             + 0.15 * annot.symbol_concordance)
    if annot.coverage_status == "fallback":
        score = max(0.0, score - 0.1)
    missing_ends = {"5_partial": 1, "3_partial": 1, "partial": 2}.get(annot.completeness, 0)
    score = max(0.0, score - FRAGMENT_PENALTY * missing_ends)
    rule = ORTHOLOG_CONFIDENCE.get(ortholog)
    if rule:
        if rule[0] == "add":
            score = score + rule[1]
        elif rule[0] == "floor":
            score = max(score, rule[1])
        else:
            score = min(max(score, rule[1]), rule[2])
        score = min(1.0, max(0.0, score))
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
_worker_go_source_hits = 5
_worker_go_min_hits = 5
_worker_go_single_cov = 0.7
_worker_go_single_identity = 50.0
# Source hits counted at most for the majority rule of the GO support requirement
GO_MAX_AVAILABLE = 10
_worker_reference_taxid = ""


def _init_worker(go_mapping, go_hierarchy, bad_name_filter, resolver, word_weights=None,
                 forbidden_go=None, go_iea_min_hits=1, reference_taxid="", go_source_hits=10, go_min_hits=5,
                 go_single_cov=0.7, go_single_identity=50.0):
    """Initialize shared read-only data in each worker process."""
    global _worker_go_mapping, _worker_go_hierarchy, _worker_bad_name_filter, _worker_resolver
    global _worker_word_weights, _worker_forbidden_go, _worker_go_iea_min_hits, _worker_reference_taxid
    global _worker_go_source_hits, _worker_go_min_hits, _worker_go_single_cov, _worker_go_single_identity
    _worker_reference_taxid = str(reference_taxid or "")
    _worker_go_source_hits = max(1, int(go_source_hits or 10))
    _worker_go_min_hits = max(1, int(go_min_hits or 1))
    _worker_go_single_cov = float(go_single_cov)
    _worker_go_single_identity = float(go_single_identity)
    _worker_go_mapping = go_mapping
    _worker_go_hierarchy = go_hierarchy
    _worker_bad_name_filter = bad_name_filter
    _worker_resolver = resolver
    _worker_word_weights = word_weights
    _worker_forbidden_go = frozenset(forbidden_go or ())
    _worker_go_iea_min_hits = go_iea_min_hits


_ANCESTOR_CACHE = {}
_EXPERIMENTAL_CODES = frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                 "HTP", "HDA", "HMP", "HGI", "HEP"})


def _with_ancestors(terms: set, go_hierarchy) -> set:
    """A set of GO terms with all their ancestors (cached per term within a worker)."""
    out = set()
    for t in terms:
        norm = go_hierarchy.normalize_id(t) if go_hierarchy is not None else t
        if not norm:
            continue
        anc = _ANCESTOR_CACHE.get(norm)
        if anc is None:
            anc = go_hierarchy.get_ancestors(norm) if go_hierarchy is not None else set()
            _ANCESTOR_CACHE[norm] = anc
        out.add(norm)
        out |= anc
    return out


def filter_weak_go(annotations: dict, go_hierarchy, min_hits: int = 3, iea_min_hits: int = 2) -> dict:
    """
    Remove GO terms that rest on weak support, in place,
    once every source of GO terms has been added. A term confirmed by InterPro
    (go_interpro) is kept, and so is a term from the reference gene
    (go_reference), unless the gene's hits do not carry it, its only evidence
    for the reference gene is electronic (IEA), and no close relative's
    one-to-one ortholog carries it (go_ortholog_support). A term
    from the hits alone is kept when one of its annotations is experimental
    (EXP, IDA, IPI, IMP, IGI, IEP, HTP codes), or when a SwissProt entry of the
    gene itself (its symbol is the gene's symbol) carries it with non-electronic
    evidence (annot.go_curated); otherwise it needs at least
    min_hits supporting source hits, or a majority of the gene's source hits
    with GO annotations when that is smaller (never fewer than 2; annot.go_need,
    set per gene by the annotation), and a term whose only evidence is
    electronic (IEA) at least iea_min_hits. A gene whose only source hit covers
    both proteins well at high identity (annot.go_single_strong) keeps that
    hit's non-electronic terms. A hit supports a term when it is
    annotated with the term or a more specific one. Terms without recorded
    evidence codes (older databases) follow min_hits. A removed term is replaced
    by its most specific ancestors that pass (annot.go_fallback, computed from the
    source hits). Removed terms are kept in annot.go_weak_removed for the
    evidence report. Returns counts.
    """
    from goannotate.go import EVIDENCE_CLASSES
    experimental = EVIDENCE_CLASSES["experimental"]
    stats = Counter()
    for annot in annotations.values():
        if not annot.go_terms:
            continue
        removed = {}
        for t in list(annot.go_terms):
            if annot.go_interpro.get(t) or t in annot.go_curated:
                continue
            support = annot.go_support.get(t, 0) + annot.go_ortholog_support.get(t, 0)
            codes = {c for c in annot.go_evidence.get(t, "").split(",") if c}
            if annot.go_reference.get(t):
                # From the reference gene: kept, except a term the gene's hits do not carry whose only
                # evidence for the reference gene is electronic, unless a close relative's ortholog carries it
                if (not annot.go_support.get(t, 0) and codes and codes <= ELECTRONIC_CODES
                        and not annot.go_ortholog_support.get(t, 0)):
                    removed[t] = (0, ",".join(sorted(codes)), "reference")
                    stats["reference_electronic"] += 1
                continue
            if codes & experimental:
                continue
            main = annot.go_need if annot.go_need else min_hits
            need = max(main, iea_min_hits) if codes and codes <= ELECTRONIC_CODES else main
            if support < need:
                removed[t] = (support, ",".join(sorted(codes)), "hits")
        if not removed:
            continue
        specific_before = set(annot.specific_go_terms)
        annot.go_terms = set(annot.go_terms) - set(removed)
        # the supported part of a removed term: its most specific ancestors that pass
        for t in removed:
            for a, support in (annot.go_fallback.get(t) or {}).items():
                if a not in annot.go_terms:
                    annot.go_terms.add(a)
                    annot.go_support[a] = support + annot.go_ortholog_support.get(a, 0)
                    annot.go_evidence[a] = removed[t][1]
                    stats["ancestors_kept"] += 1
        annot.specific_go_terms = (go_hierarchy.filter_to_specific(annot.go_terms) if go_hierarchy is not None
                                   else set(annot.go_terms))
        annot.go_weak_removed = removed
        stats["genes"] += 1
        stats["terms"] += len(removed)
        stats["specific_removed"] += len(specific_before - annot.specific_go_terms)
        if not annot.go_terms:
            stats["genes_left_without_go"] += 1
    stats = dict(stats)
    stats.update({"min_hits": min_hits, "iea_min_hits": iea_min_hits})
    n_ref = stats.get("reference_electronic", 0)
    logger.info(f"GO terms with weak support removed: {stats.get('terms', 0) - n_ref:,} from the hits alone "
                f"(fewer than {min_hits} supporting hits, {iea_min_hits} for electronic-only terms; experimental "
                f"terms and terms confirmed by the reference gene, InterPro or the gene's SwissProt entries kept) "
                f"and {n_ref:,} electronic-only terms of the reference gene that neither the hits nor a close "
                f"relative's ortholog carry; {stats.get('terms', 0):,} removed from {stats.get('genes', 0):,} genes, "
                f"{stats.get('ancestors_kept', 0):,} supported ancestors kept in their place "
                f"({stats.get('genes_left_without_go', 0):,} left without GO terms)")
    return stats


def annotate_single_gene(gene_id: str, hits_with_names: list,
                         consensus_threshold: float, namespaces: list,
                         symbol_case: str = "upper", symbol_score: str = "mean",
                         prior: tuple = None) -> tuple:
    """
    Full annotation of a single gene (runs in a worker process).

    prior: (ortholog symbol, close) from the orthology module, or None. A
    close relatives' symbol also steers the cluster choice: when some
    cluster carries a hit with that symbol, the winning cluster is chosen
    among those clusters.

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
            'go_n_source': 0,
            'best_query_coverage': 0.0,
            'best_subject_coverage': 0.0,
            'hit_names': {},
            'hit_symbols': {},
            'hit_curated_names': {},
        }

        if not hits_with_names:
            return (gene_id, result)

        result['had_hits'] = True
        result['n_isoforms'] = len(set(tx for t in hits_with_names for tx in _isoforms(t)))
        result['hits'] = sorted([t[HIT] for t in hits_with_names],
                                key=lambda h: h.bitscore, reverse=True)
        result['hit_names'] = dict(Counter(normalize_gene_name(t[NAME]) for t in hits_with_names if t[NAME]))
        result['hit_symbols'] = dict(Counter(t[SYMBOL].upper() for t in hits_with_names if t[SYMBOL]))
        # curated (SwissProt) names given by other species than the reference: the names a reference
        # protein's name is checked against (assign.assign_symbols)
        curated = defaultdict(set)      # name -> symbols of the entries carrying it
        for t in hits_with_names:
            if t[NAME] and t[HIT].is_swissprot and t[HIT].taxid != _worker_reference_taxid:
                syms = curated[normalize_gene_name(t[NAME])]
                for sym in (t[SYMBOL], t[HIT].gene_symbol):
                    if sym:
                        syms.add(sym.upper())
        result['hit_curated_names'] = {name: sorted(syms) for name, syms in curated.items()}

        # Cluster and select the winning cluster
        clusters = cluster_hits_agglomerative(hits_with_names,
                                              similarity_threshold=consensus_threshold)
        if prior and prior[0] and prior[1]:
            key = prior[0].lower()
            with_prior = [c for c in clusters
                          if any(hits_with_names[i][SYMBOL] and hits_with_names[i][SYMBOL].lower() == key for i in c)]
            if with_prior:
                clusters = with_prior
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
            score_method=symbol_score, prior=prior)

        result['consensus_name'] = protein_name
        result['consensus_symbol'] = format_symbol(gene_symbol, symbol_case)

        if protein_name:
            result['had_protein_name'] = True
        if gene_symbol:
            result['had_symbol'] = True
            result['symbol_source'] = symbol_origin

        # Concordance of the voted name and symbol with the full merged hit pool (recomputed for the
        # final annotation by final_concordance)
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

        # GO terms from the winning symbol group (its top hits by bitscore), with
        # the evidence codes behind each term and the number of hits that
        # carry it
        go_source_indices = sorted(winning_symbol_indices,
                                   key=lambda i: hits_with_names[i][HIT].bitscore,
                                   reverse=True)[:_worker_go_source_hits]
        go_hits = [hits_with_names[i][HIT] for i in go_source_indices]
        term_codes, term_support = collect_go_evidence(go_hits, go_mapping)
        hit_term_sets = [set(collect_go_evidence([h], go_mapping)[0]) for h in go_hits]
        winning_sym_upper = gene_symbol.upper() if gene_symbol else ""
        # Curated (non-electronic) annotations of SwissProt entries of the gene itself (GN= symbol equal to
        # the winning symbol): taken without a support requirement
        go_curated = set()

        def add_curated(hit, symbol, codes_by_term):
            if not (hit.is_swissprot and winning_sym_upper and symbol and symbol.upper() == winning_sym_upper):
                return
            for t, codes in codes_by_term.items():
                if codes - ELECTRONIC_CODES:
                    go_curated.add(t)

        for i in go_source_indices:
            add_curated(hits_with_names[i][HIT], hits_with_names[i][SYMBOL],
                        collect_go_evidence([hits_with_names[i][HIT]], go_mapping)[0])
        go_source = f"cluster_top{_worker_go_source_hits}" if term_codes else ""
        source_with_go = [h for h in go_hits if any(acc in go_mapping for acc in extract_accession(h.subject_id))]
        n_source = len(source_with_go)

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
                    n_source += 1
                    source_with_go.append(hit)
                    hit_term_sets.append(set(sp_codes))
                    if sym_match:
                        add_curated(hit, symbol, sp_codes)
                    if any(t not in term_codes for t in sp_codes):
                        go_source = f"cluster_top{_worker_go_source_hits}+swissprot_supplement"
                    for t, codes in sp_codes.items():
                        term_codes.setdefault(t, set()).update(codes)
                        term_support[t] = term_support.get(t, 0) + sp_support[t]

        result['go_source'] = go_source
        result['go_n_source'] = n_source
        result['go_source_accs'] = {acc for h in source_with_go for acc in extract_accession(h.subject_id)}
        # Support of a term: the source hits annotated with it or with a more specific term (an
        # annotation to a term implies its ancestors)
        closures = [_with_ancestors(s, go_hierarchy) for s in hit_term_sets if s]
        term_support = {t: sum(1 for c in closures if t in c) for t in term_codes}
        result['go_curated'] = go_curated
        # Support required of a term from the hits alone: --go-min-hits, or a majority of the source
        # hits with GO annotations (at most GO_MAX_AVAILABLE counted) when that is smaller, never fewer
        # than 2; electronic-only terms --go-iea-min-hits. A gene whose only source hit is a near
        # full-length alignment at high identity keeps that hit's non-electronic terms.
        available = min(n_source, GO_MAX_AVAILABLE)
        need_main = max(2, min(_worker_go_min_hits, available // 2 + 1)) if _worker_go_min_hits > 1 else 1
        single_strong = (n_source == 1 and source_with_go[0].query_coverage >= _worker_go_single_cov
                         and source_with_go[0].subject_coverage >= _worker_go_single_cov
                         and source_with_go[0].pident >= _worker_go_single_identity)
        if single_strong:
            need_main = 1
        result['go_need'] = need_main
        result['go_single_strong'] = single_strong

        def needed(codes):
            return max(need_main, _worker_go_iea_min_hits) if codes and codes <= ELECTRONIC_CODES else need_main

        # For a term that may fail the support requirement, its most specific ancestors that pass it
        # (kept in its place when it is removed: the supported part of the claim)
        fallback = {}
        if go_hierarchy is not None and (need_main > 1 or _worker_go_iea_min_hits > 1):
            for t, codes in term_codes.items():
                if t in go_curated or codes & _EXPERIMENTAL_CODES:
                    continue
                need = needed(codes)
                if term_support.get(t, 0) >= need:
                    continue
                norm = go_hierarchy.normalize_id(t)
                if not norm:
                    continue
                passing = {a for a in _with_ancestors({norm}, go_hierarchy) - {norm}
                           if a not in GO_ROOT_TERMS and sum(1 for c in closures if a in c) >= need}
                if passing:
                    fallback[t] = {a: sum(1 for c in closures if a in c)
                                   for a in go_hierarchy.filter_to_specific(passing)}
        result['go_fallback'] = fallback

        go_terms = {t for t in term_codes if t not in GO_ROOT_TERMS}

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


def final_concordance(annotations: dict, consensus_threshold: float = 0.5) -> dict:
    """
    Concordance of each gene's final annotation with its own hits, computed
    after every step that changes a symbol or product (reference checks,
    orthology, InterProScan products), and the confidence recomputed from it.

    Name concordance: fraction of the gene's hits whose name matches the
    final product (name similarity at least consensus_threshold).
    Symbol concordance: fraction of the gene's hits that carry a symbol and
    name the gene of the final symbol: the symbol itself, a copy form of it
    or the gene it is a copy form of (one genome suffix: kdm2aa, KDM2A), the
    replaced voted symbol, or one of the gene's aliases (hit symbols of the
    same gene, reference-species synonyms). Hits without a symbol say nothing
    about the symbol and are not counted.
    Genes without hits get 0 and are left out of the statistics.
    Returns the statistics for the summary.
    """
    names, symbols = [], []
    for annot in annotations.values():
        annot.name_concordance = 0.0
        annot.symbol_concordance = 0.0
        if annot.hits and annot.consensus_name and annot.hit_names:
            final = normalize_gene_name(annot.consensus_name)
            total = sum(annot.hit_names.values())
            match = sum(n for name, n in annot.hit_names.items()
                        if name_similarity(final, name) >= consensus_threshold)
            annot.name_concordance = match / total if total else 0.0
            names.append(annot.name_concordance)
        if annot.hits and annot.consensus_symbol:
            final = annot.consensus_symbol.upper()
            same = {final} | symbol_bases(final)
            if annot.ortholog_evidence != "close_taken":
                # a symbol taken from close relatives names another gene than the vote:
                # the replaced symbol and the vote's aliases are not the same gene
                same |= {a.upper() for a in annot.symbol_aliases if a}
                if annot.voted_symbol:
                    same.add(annot.voted_symbol.upper())
            total = sum(annot.hit_symbols.values())
            match = sum(n for sym, n in annot.hit_symbols.items()
                        if sym in same or final in symbol_bases(sym))
            annot.symbol_concordance = match / total if total else 0.0
            symbols.append(annot.symbol_concordance)
        annot.confidence, annot.confidence_score = compute_confidence(annot)

    def describe(values):
        if not values:
            return {}
        return {"n": len(values), "mean": statistics.mean(values), "median": statistics.median(values),
                "min": min(values), "max": max(values)}
    stats = {"name": describe(names), "symbol": describe(symbols),
             "confidence": dict(Counter(a.confidence for a in annotations.values()))}
    for label, d in (("Name", stats["name"]), ("Symbol", stats["symbol"])):
        if d:
            logger.info(f"  {label} concordance with the final annotation ({d['n']:,} genes with hits): "
                        f"mean={d['mean']:.2f}, median={d['median']:.2f}")
    return stats


def propagate_symbols(annotations: dict, gene_hits: dict, symbol_case: str = "upper") -> tuple:
    """
    Symbol propagation: genes sharing an identical normalized protein name
    receive the most common symbol among them that their own hits carry, when
    they lack one. Returns (number propagated, propagated gene IDs).
    """
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
        ranked = [s for s, _ in Counter(existing_symbols).most_common()]
        for gid in genes_without:
            # Only a symbol the gene's own hits carry: identical product names
            # alone ("disks large 5") are not evidence that two genes are
            # copies of the same gene. The genes sharing the name may carry a
            # close relative's form (glula) that the hits give in the
            # reference form (GLUL); the hit's form is then given, since the
            # relatives' suffix designates one copy and this gene is another.
            hit_symbols = {}
            for t in gene_hits.get(gid, []):
                if t[SYMBOL]:
                    hit_symbols.setdefault(t[SYMBOL].lower(), t[SYMBOL])
            chosen = ""
            for s in ranked:
                if s.lower() in hit_symbols:
                    chosen = s
                    break
            if not chosen:
                for s in ranked:
                    bases = symbol_bases(s)
                    forms = [h for h in hit_symbols.values() if h.upper() in bases]
                    if forms:
                        chosen = format_symbol(forms[0], symbol_case)
                        break
            if not chosen:
                continue
            annotations[gid].consensus_symbol = chosen
            symbols_propagated += 1
            propagated_gene_ids.append(gid)
    return symbols_propagated, propagated_gene_ids


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
                     coverage_status: dict = None,
                     ortholog_prior: dict = None, reference_taxid: str = "",
                     go_source_hits: int = 10, go_min_hits: int = 5, go_single_cov: float = 0.7,
                     go_single_identity: float = 50.0) -> dict:
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
    ortholog_prior: gene -> (ortholog symbol, close) from the orthology
    module, steering the symbol choice (see select_name_and_symbol).

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

    ortholog_prior = ortholog_prior or {}
    if ortholog_prior:
        n_close = sum(1 for _, close in ortholog_prior.values() if close)
        logger.info(f"Ortholog symbols steer the vote for {len(ortholog_prior):,} genes "
                    f"({n_close:,} from close relatives decide, the rest count as one more vote)")
    tasks = [(gene_id, gene_hits.get(gene_id, []), consensus_threshold, namespaces, symbol_case, symbol_score,
              ortholog_prior.get(gene_id))
             for gene_id in all_gene_ids]

    start_time = time.time()

    # The main process uses the same globals (normalize_gene_name reads the
    # bad-name filter for comparison-only normalization during propagation).
    init_args = (go_mapping, go_hierarchy, bad_name_filter, resolver, word_weights, forbidden_go, go_iea_min_hits,
                 reference_taxid, go_source_hits, go_min_hits, go_single_cov, go_single_identity)
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
        'symbols_from_relative': 0,
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
        annot.hit_names = result_dict.get('hit_names', {})
        annot.hit_symbols = result_dict.get('hit_symbols', {})
        annot.hit_curated_names = result_dict.get('hit_curated_names', {})
        annot.go_n_source = result_dict.get('go_n_source', 0)
        annot.go_curated = result_dict.get('go_curated', set())
        annot.go_source_accs = result_dict.get('go_source_accs', set())
        annot.go_fallback = result_dict.get('go_fallback', {})
        annot.go_need = result_dict.get('go_need', 0)
        annot.go_single_strong = result_dict.get('go_single_strong', False)
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
            elif result_dict['symbol_source'] == 'relative':
                stats['symbols_from_relative'] += 1
        if result_dict['go_terms']:
            stats['genes_with_go_terms'] += 1
        if result_dict['cluster_fraction'] > 0:
            cluster_fractions.append(result_dict['cluster_fraction'])

    symbols_propagated, propagated_gene_ids = propagate_symbols(annotations, gene_hits, symbol_case)

    stats['symbols_from_propagation'] = symbols_propagated
    if symbols_propagated > 0:
        logger.info(f"  Symbol propagation: {symbols_propagated:,} genes received symbols "
                    f"from genes with identical protein names")

    for gid in propagated_gene_ids:
        annotations[gid].symbol_status = "propagated"

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
                     + stats['symbols_from_ncbi'] + stats['symbols_from_propagation'] + stats['symbols_from_relative'])
    logger.info(f"  Coding genes with gene symbol: {total_symbols:,} ({100*total_symbols/coding:.1f}%)")
    if stats['symbols_from_ortholog'] > 0:
        logger.info(f"    - From reference-species ortholog: {stats['symbols_from_ortholog']:,}")
    if stats['symbols_from_uniprot'] > 0:
        logger.info(f"    - From UniProt GN= field: {stats['symbols_from_uniprot']:,}")
    if stats['symbols_from_ncbi'] > 0:
        logger.info(f"    - From NCBI cross-reference: {stats['symbols_from_ncbi']:,}")
    if stats['symbols_from_relative'] > 0:
        logger.info(f"    - Decided by a close relative's ortholog (hits carrying its symbol): "
                    f"{stats['symbols_from_relative']:,}")
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

    return annotations
