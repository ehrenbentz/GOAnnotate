"""Parsing and per-transcript filtering of BLAST/Diamond results."""

import logging
from collections import defaultdict
from dataclasses import dataclass

from goannotate.filters import BadNameFilter, NameAudit
from goannotate.models import BlastHit, BlastStats

logger = logging.getLogger(__name__)


def parse_blast_results(blast_file: str, evalue_threshold: float = 1e-5) -> tuple:
    """
    Parse BLAST/Diamond output (format 6 with stitle) keyed by query ID.

    Collects all hits per query, applying only e-value filtering.
    No sorting or truncation is performed here.

    Two column layouts are read: the current layout written by search.py
    (... evalue bitscore qlen slen stitle) and the earlier one without
    qlen and slen (... evalue bitscore stitle). A line is taken to carry
    lengths when columns 13 and 14 are integers and a title follows.

    Returns:
        tuple of (dict of query_id -> list of BlastHit, BlastStats)
    """
    logger.info(f"Parsing BLAST results from {blast_file}")
    logger.info(f"E-value threshold: {evalue_threshold}")

    transcript_hits = defaultdict(list)
    all_queries_in_file = set()
    stats = BlastStats()
    with_lengths = without_lengths = 0

    with open(blast_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue

            stats.total_hits_in_file += 1
            query_id = parts[0]
            all_queries_in_file.add(query_id)
            evalue = float(parts[10])

            if evalue > evalue_threshold:
                stats.hits_filtered_by_evalue += 1
                continue

            if len(parts) >= 15 and parts[12].isdigit() and parts[13].isdigit():
                qlen, slen = int(parts[12]), int(parts[13])
                stitle = '\t'.join(parts[14:])
                with_lengths += 1
            else:
                qlen = slen = 0
                stitle = '\t'.join(parts[12:]) if len(parts) > 12 else ""
                without_lengths += 1

            hit = BlastHit(
                query_id=query_id,
                subject_id=parts[1],
                pident=float(parts[2]),
                length=int(parts[3]),
                mismatch=int(parts[4]),
                gapopen=int(parts[5]),
                qstart=int(parts[6]),
                qend=int(parts[7]),
                sstart=int(parts[8]),
                send=int(parts[9]),
                evalue=evalue,
                bitscore=float(parts[11]),
                stitle=stitle,
                qlen=qlen,
                slen=slen,
            )

            transcript_hits[query_id].append(hit)

            if hit.is_swissprot:
                stats.swissprot_hits += 1
            elif hit.is_trembl:
                stats.trembl_hits += 1
            else:
                stats.other_hits += 1

    stats.queries_in_file = len(all_queries_in_file)
    stats.queries_after_evalue_filter = len(transcript_hits)
    stats.hits_after_evalue_filter = stats.swissprot_hits + stats.trembl_hits + stats.other_hits
    stats.hits_with_lengths = with_lengths
    if without_lengths and with_lengths:
        logger.warning(f"  {without_lengths:,} hits lack the qlen/slen columns and {with_lengths:,} have them; "
                       f"the coverage filter can only test the latter")
    elif without_lengths:
        logger.info("  The results carry no qlen/slen columns (searched by an earlier version); "
                    "hit coverage cannot be tested")

    logger.info(f"Parsed hits for {len(transcript_hits):,} queries")
    logger.info(f"  Total hits in file: {stats.total_hits_in_file:,}")
    if stats.hits_filtered_by_evalue > 0:
        logger.info(f"  Filtered by e-value (>{evalue_threshold}): {stats.hits_filtered_by_evalue:,}")
    logger.info(f"  Hits passing e-value filter: {stats.hits_after_evalue_filter:,}")
    logger.info(f"  SwissProt (sp|): {stats.swissprot_hits:,}")
    logger.info(f"  TrEMBL (tr|): {stats.trembl_hits:,}")
    if stats.other_hits > 0:
        logger.info(f"  Other: {stats.other_hits:,}")

    return dict(transcript_hits), stats


@dataclass
class HitCoverageFilter:
    """
    Alignment coverage and identity requirements for a hit to take part in
    naming and symbol voting (after PANNZER2's hit filter).

    mode:
      strict      the alignment must cover min_query_cov of the query AND
                  min_subject_cov of the subject
      permissive  it must cover min_query_cov of the query OR min_subject_cov
                  of the subject, OR be at least min_align_length residues long
                  (fragments and multidomain proteins pass)
      off         no coverage or identity requirement
    In both active modes the identity must be at least min_identity percent.
    Hits without recorded lengths (results of searches run by earlier versions) are
    tested on identity only.

    fallback: when no hit of a transcript passes, its unfiltered hits are
    used and the transcript is marked (coverage status "fallback"); without
    fallback such transcripts keep no hits.
    """
    mode: str = "permissive"
    min_query_cov: float = 0.5
    min_subject_cov: float = 0.5
    min_identity: float = 40.0
    min_align_length: int = 100
    fallback: bool = True

    @property
    def active(self) -> bool:
        return self.mode != "off"

    def passes(self, hit: BlastHit) -> bool:
        if not self.active:
            return True
        if hit.pident < self.min_identity:
            return False
        if not hit.has_lengths:
            return True
        qcov, scov = hit.query_coverage, hit.subject_coverage
        if self.mode == "strict":
            return qcov >= self.min_query_cov and scov >= self.min_subject_cov
        return (qcov >= self.min_query_cov or scov >= self.min_subject_cov
                or hit.length >= self.min_align_length)

    def describe(self) -> str:
        if not self.active:
            return "off"
        if self.mode == "strict":
            s = (f"strict: query coverage >= {self.min_query_cov:.2f} and subject coverage >= "
                 f"{self.min_subject_cov:.2f}")
        else:
            s = (f"permissive: query coverage >= {self.min_query_cov:.2f} or subject coverage >= "
                 f"{self.min_subject_cov:.2f} or alignment >= {self.min_align_length} residues")
        s += f", identity >= {self.min_identity:g}%"
        s += "; transcripts without a passing hit fall back to all their hits" if self.fallback \
            else "; transcripts without a passing hit keep none"
        return s


def filter_and_select_hits(transcript_hits: dict,
                           bad_name_filter: BadNameFilter,
                           top_n: int = 100,
                           audit: NameAudit = None,
                           coverage: HitCoverageFilter = None) -> tuple:
    """
    For each transcript, clean names, filter bad names, apply the hit
    coverage filter, and select the top-N hits by bitscore. Top-N
    truncation happens after filtering so that uninformative hits do not
    consume slots.

    The GN= symbol is carried as found: symbols are tested against the
    symbol patterns only after ortholog resolution (resolve_symbols), so
    that a LOC number or locus tag can still reach the reference-species
    symbol of its gene.

    Returns:
        (dict of transcript_id -> list of (BlastHit, cleaned_name, gn_symbol) tuples,
         dict of transcript_id -> coverage status: "pass", "fallback" or "" when the filter is off)
    """
    logger.info(f"Filtering and selecting top-{top_n} hits per transcript")
    coverage = coverage or HitCoverageFilter(mode="off")
    if coverage.active:
        logger.info(f"  Hit coverage filter {coverage.describe()}")

    filtered = {}
    status = {}
    total_input = 0
    total_surviving = 0
    total_after_topn = 0
    cov_failed = 0
    n_fallback = n_dropped = 0

    for transcript_id, hits in transcript_hits.items():
        total_input += len(hits)

        surviving = []
        for hit in hits:
            cleaned_name = bad_name_filter.process_name(hit.gene_name, audit)
            if cleaned_name is None:
                continue
            surviving.append((hit, cleaned_name, hit.gene_symbol))

        total_surviving += len(surviving)
        if not surviving:
            continue

        if coverage.active:
            covered = [t for t in surviving if coverage.passes(t[0])]
            cov_failed += len(surviving) - len(covered)
            if covered:
                surviving = covered
                status[transcript_id] = "pass"
            elif coverage.fallback:
                status[transcript_id] = "fallback"
                n_fallback += 1
            else:
                n_dropped += 1
                continue
        else:
            status[transcript_id] = ""

        surviving.sort(key=lambda t: -t[0].bitscore)
        selected = surviving[:top_n]
        total_after_topn += len(selected)
        filtered[transcript_id] = selected

    logger.info(f"  Transcripts with surviving hits: {len(filtered):,}")
    logger.info(f"  Hits: {total_input:,} input -> "
                f"{total_surviving:,} after name filtering -> "
                f"{total_surviving - cov_failed:,} after the coverage filter -> "
                f"{total_after_topn:,} after top-{top_n}")
    if coverage.active:
        n_pass = sum(1 for s in status.values() if s == "pass")
        logger.info(f"  Coverage filter: {cov_failed:,} hits failed; {n_pass:,} transcripts have passing hits, "
                    f"{n_fallback:,} fell back to their unfiltered hits, {n_dropped:,} were left without hits")

    return filtered, status


def collect_symbol_pairs(filtered_transcript_hits: dict, reference_taxid: str) -> set:
    """(tax_id, lowercase GN symbol) pairs of filtered hits from species other than the reference."""
    pairs = set()
    for hits in filtered_transcript_hits.values():
        for hit, _name, symbol in hits:
            taxid = hit.taxid
            if symbol and taxid and taxid != reference_taxid:
                pairs.add((taxid, symbol.lower()))
    return pairs


def resolve_symbols(filtered_transcript_hits: dict, resolver, bad_name_filter: BadNameFilter,
                    audit: NameAudit = None) -> dict:
    """
    Assign the final symbol of every filtered hit through the SymbolResolver
    (reference-species ortholog symbol where available, else the hit's own
    GN= symbol). Symbols that match the symbol patterns (LOC numbers, locus
    tags, clone IDs) are blanked here, after resolution was attempted.

    Returns:
        dict of transcript_id -> list of (BlastHit, cleaned_name, symbol, symbol_origin) tuples
    """
    resolved = {}
    counts = {"ortholog": 0, "uniprot": 0, "": 0}
    raw_bad = rescued = blanked = 0
    for transcript_id, hits in filtered_transcript_hits.items():
        out = []
        for hit, cleaned_name, gn_symbol in hits:
            symbol, origin = resolver.resolve(hit, gn_symbol, bad_name_filter.is_bad_symbol)
            counts[origin] += 1
            if gn_symbol and bad_name_filter.is_bad_symbol(gn_symbol):
                raw_bad += 1
                if origin == "ortholog":
                    rescued += 1
                else:
                    blanked += 1
                    if audit is not None:
                        reason = bad_name_filter.explain_bad_symbol(gn_symbol)
                        audit.symbol_discards[reason] += 1
                        if len(audit.symbol_discard_examples[reason]) < 50:
                            audit.symbol_discard_examples[reason][gn_symbol] += 1
            if audit is not None and symbol:
                audit.kept_symbols[symbol] += 1
            out.append((hit, cleaned_name, symbol, origin))
        resolved[transcript_id] = out
    if audit is not None:
        audit.symbols_raw_bad = raw_bad
        audit.symbols_rescued = rescued
    total = sum(counts.values()) or 1
    if resolver.has_orthologs:
        logger.info(f"Symbol resolution: {counts['ortholog']:,} hits ({100*counts['ortholog']/total:.1f}%) "
                    f"use the reference-species ortholog symbol, {counts['uniprot']:,} keep their GN= symbol, "
                    f"{counts['']:,} have none")
    else:
        logger.info(f"Symbols: {counts['uniprot']:,} hits with a GN= symbol, {counts['']:,} without")
    if raw_bad:
        logger.info(f"  GN= symbols that are genome-specific identifiers: {raw_bad:,}; "
                    f"{rescued:,} resolved to a reference-species symbol, {blanked:,} blanked")
    return resolved


def merge_hits_by_gene(resolved_transcript_hits: dict, transcript_to_gene: dict) -> dict:
    """
    Merge filtered hits from all transcripts of each gene into a single pool.
    Each hit is tagged with its source transcript ID.

    Returns:
        dict of gene_id -> list of
        (BlastHit, cleaned_name, symbol, symbol_origin, source_transcript) tuples
    """
    gene_hits = defaultdict(list)
    for transcript_id, hits in resolved_transcript_hits.items():
        gene_id = transcript_to_gene.get(transcript_id, transcript_id)
        for hit, cleaned_name, symbol, origin in hits:
            gene_hits[gene_id].append((hit, cleaned_name, symbol, origin, transcript_id))
    return dict(gene_hits)
