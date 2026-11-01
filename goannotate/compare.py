"""Comparison of GOAnnotate results with functional annotation already present in the GFF3.

Nothing is resolved here. The comparison reports, for every gene, whether
the GOAnnotate symbol and product agree with the existing ones, and
summarizes agreement by any gene- or transcript-level attributes that
record how the gene model was produced or changed (for example GRAFT's
graft_status and graft_protein).
"""

import logging
import re
from collections import Counter, defaultdict
from pathlib import Path

from goannotate.annotate import name_similarity
from goannotate.output import product_for

logger = logging.getLogger(__name__)

PRODUCT_MATCH = 0.8
PRODUCT_PARTIAL = 0.5

_PRODUCT_NOISE = [
    re.compile(r",?\s*transcript variant\s+\S+", re.I),
    re.compile(r",?\s*isoform\s+\S+", re.I),
    re.compile(r"\s*\(fragment\)", re.I),
    re.compile(r"-like\b", re.I),
    re.compile(r"\bprecursor\b", re.I),
    re.compile(r"\bhomolog\b", re.I),
    re.compile(r"\bputative\b", re.I),
    re.compile(r"\bprobable\b", re.I),
]


def normalize_product(name: str) -> str:
    """Lowercase, remove variant/isoform qualifiers and punctuation, collapse whitespace."""
    if not name:
        return ""
    s = name
    for pat in _PRODUCT_NOISE:
        s = pat.sub("", s)
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def existing_symbol(gene, symbol_attrs: list) -> str:
    for attr in symbol_attrs:
        value = gene.attrs.get(attr)
        if value:
            return value
    for tx in gene.transcripts:
        for attr in symbol_attrs:
            value = tx.attrs.get(attr)
            if value:
                return value
    return ""


def existing_product(gene, product_attrs: list) -> str:
    """
    Existing product name of a gene. Attributes are tried in the order
    given; for each attribute the gene feature is checked first, then the
    coding transcripts, then any transcript. With the default order
    (product, description) a transcript product therefore takes
    precedence over a gene description.
    """
    for attr in product_attrs:
        value = gene.attrs.get(attr)
        if value:
            return value
        for coding_only in (True, False):
            for tx in gene.transcripts:
                if coding_only and not tx.coding:
                    continue
                value = tx.attrs.get(attr)
                if value:
                    return value
    return ""


def stratum_value(gene, attr: str, level: str) -> str:
    """
    Value of a stratification attribute for a gene. Transcript-level
    attributes are collapsed: one distinct value is reported as is,
    several as mixed(a|b).
    """
    if level == "gene":
        return gene.attrs.get(attr, "") or "NA"
    values = sorted({tx.attrs.get(attr, "") for tx in gene.transcripts} - {""})
    if not values:
        return "NA"
    if len(values) == 1:
        return values[0]
    return "mixed(" + "|".join(values) + ")"


def symbol_status(existing: str, new: str) -> str:
    if existing and new:
        return "match" if existing.lower() == new.lower() else "mismatch"
    if new:
        return "goannotate_only"
    if existing:
        return "existing_only"
    return "neither"


def product_status(existing: str, new: str) -> tuple:
    """Returns (status, similarity)."""
    if existing and new:
        sim = name_similarity(normalize_product(existing), normalize_product(new))
        if sim >= PRODUCT_MATCH:
            return "match", sim
        if sim >= PRODUCT_PARTIAL:
            return "partial", sim
        return "mismatch", sim
    if new:
        return "goannotate_only", 0.0
    if existing:
        return "existing_only", 0.0
    return "neither", 0.0


def compare_with_existing(annotations: dict, genes: list, output_dir: Path, prefix: str,
                          symbol_attrs: list, product_attrs: list,
                          gene_strata: list, transcript_strata: list, is_bad_symbol=None) -> dict:
    """
    Write {prefix}_comparison.tsv (one row per gene) and
    {prefix}_comparison_summary.txt. Returns summary counts.

    is_bad_symbol: existing symbols it recognizes as genome-specific
    identifiers (LOC numbers of an NCBI annotation, locus tags) are listed
    in the table but compared as absent, so that a GOAnnotate symbol on a
    LOC gene counts as goannotate_only rather than a mismatch.
    """
    output_dir = Path(output_dir)
    table_path = output_dir / f"{prefix}_comparison.tsv"
    summary_path = output_dir / f"{prefix}_comparison_summary.txt"
    review_path = output_dir / f"{prefix}_review.tsv"
    review_rows = []

    # Keep only strata attributes that occur at least once
    present_gene_strata = [a for a in gene_strata
                           if any(a in g.attrs for g in genes)]
    present_tx_strata = [a for a in transcript_strata
                         if any(a in t.attrs for g in genes for t in g.transcripts)]
    strata = [(a, "gene") for a in present_gene_strata] + [(a, "transcript") for a in present_tx_strata]

    n_existing_symbol = sum(1 for g in genes if existing_symbol(g, symbol_attrs))
    n_existing_identifier = sum(1 for g in genes if existing_symbol(g, symbol_attrs) and is_bad_symbol
                                and is_bad_symbol(existing_symbol(g, symbol_attrs)))
    n_existing_product = sum(1 for g in genes if existing_product(g, product_attrs))
    logger.info(f"Comparing with existing annotation: {n_existing_symbol:,} genes carry a symbol "
                f"({'/'.join(symbol_attrs)}), {n_existing_product:,} carry a product "
                f"({'/'.join(product_attrs)})")
    if strata:
        logger.info("  Stratifying by: " + ", ".join(f"{a} ({lvl})" for a, lvl in strata))

    symbol_counts = Counter()
    product_counts = Counter()
    by_stratum = {key: defaultdict(Counter) for key in strata}       # (attr, level) -> value -> Counter(symbol_status)
    by_stratum_product = {key: defaultdict(Counter) for key in strata}
    stratum_totals = {key: Counter() for key in strata}

    with open(table_path, "w") as fh:
        header = ["gene_id", "gene_type", "coding"]
        header += [a for a, _ in strata]
        header += ["existing_symbol", "goannotate_symbol", "symbol_status",
                   "existing_product", "goannotate_product", "product_similarity", "product_status",
                   "confidence", "n_go_terms", "assignment"]
        fh.write("\t".join(header) + "\n")

        for gene in genes:
            annot = annotations.get(gene.key)
            new_symbol = annot.consensus_symbol if annot else ""
            new_product = product_for(annot)
            old_symbol = existing_symbol(gene, symbol_attrs)
            old_product = existing_product(gene, product_attrs)

            compared = "" if (old_symbol and is_bad_symbol and is_bad_symbol(old_symbol)) else old_symbol
            s_status = symbol_status(compared, new_symbol)
            p_status, p_sim = product_status(old_product, new_product)
            symbol_counts[s_status] += 1
            product_counts[p_status] += 1

            row = [gene.key, gene.biotype, "yes" if gene.coding else "no"]
            for key in strata:
                value = stratum_value(gene, key[0], key[1])
                row.append(value)
                by_stratum[key][value][s_status] += 1
                by_stratum_product[key][value][p_status] += 1
                stratum_totals[key][value] += 1
            assignment = annot.symbol_status if annot else ""
            row += [old_symbol, new_symbol, s_status, old_product, new_product,
                    f"{p_sim:.2f}", p_status,
                    annot.confidence if annot else "", str(len(annot.specific_go_terms)) if annot else "0",
                    assignment]
            fh.write("\t".join(row) + "\n")
            if s_status == "mismatch":
                review_rows.append(row)

    # Review table: genes where both annotations name the gene and disagree.
    # The existing symbol is kept in the annotated GFF3 (fill mode); the
    # GOAnnotate symbol is recorded as goannotate_symbol.
    with open(review_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in review_rows:
            fh.write("\t".join(row) + "\n")

    symbol_order = ["match", "mismatch", "goannotate_only", "existing_only", "neither"]
    product_order = ["match", "partial", "mismatch", "goannotate_only", "existing_only", "neither"]

    def _agreement(counter: Counter, statuses: list) -> str:
        both = sum(counter.get(s, 0) for s in statuses)
        if both == 0:
            return "n/a"
        return f"{100 * counter.get('match', 0) / both:.1f}%"

    with open(summary_path, "w") as fh:
        fh.write("GOAnnotate comparison with existing annotation\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Genes: {len(genes):,}\n")
        fh.write(f"Genes with an existing symbol ({'/'.join(symbol_attrs)}): {n_existing_symbol:,}\n")
        fh.write(f"Genes with an existing product ({'/'.join(product_attrs)}): {n_existing_product:,}\n\n")

        fh.write("Gene symbol (case-insensitive)\n")
        fh.write("-" * 70 + "\n")
        for s in symbol_order:
            fh.write(f"  {s:<18}{symbol_counts.get(s, 0):>10,}\n")
        fh.write(f"  agreement where both present: {_agreement(symbol_counts, ['match', 'mismatch'])}\n\n")

        fh.write(f"Product name (normalized similarity; match >= {PRODUCT_MATCH}, partial >= {PRODUCT_PARTIAL})\n")
        fh.write("-" * 70 + "\n")
        for s in product_order:
            fh.write(f"  {s:<18}{product_counts.get(s, 0):>10,}\n")
        fh.write(f"  agreement where both present: "
                 f"{_agreement(product_counts, ['match', 'partial', 'mismatch'])}\n\n")

        for key in strata:
            attr, level = key
            fh.write(f"Symbol status by {attr} ({level} level)\n")
            fh.write("-" * 70 + "\n")
            fh.write(f"  {'value':<32}{'n':>8}" + "".join(f"{s[:12]:>14}" for s in symbol_order)
                     + f"{'agree':>10}\n")
            for value, _ in stratum_totals[key].most_common():
                c = by_stratum[key][value]
                fh.write(f"  {value[:32]:<32}{stratum_totals[key][value]:>8,}"
                         + "".join(f"{c.get(s, 0):>14,}" for s in symbol_order)
                         + f"{_agreement(c, ['match', 'mismatch']):>10}\n")
            fh.write("\n")
            fh.write(f"Product status by {attr} ({level} level)\n")
            fh.write("-" * 70 + "\n")
            fh.write(f"  {'value':<32}{'n':>8}" + "".join(f"{s[:12]:>14}" for s in product_order)
                     + f"{'agree':>10}\n")
            for value, _ in stratum_totals[key].most_common():
                c = by_stratum_product[key][value]
                fh.write(f"  {value[:32]:<32}{stratum_totals[key][value]:>8,}"
                         + "".join(f"{c.get(s, 0):>14,}" for s in product_order)
                         + f"{_agreement(c, ['match', 'partial', 'mismatch']):>10}\n")
            fh.write("\n")

    logger.info(f"  Symbol: " + ", ".join(f"{s}={symbol_counts.get(s, 0):,}" for s in symbol_order))
    logger.info(f"  Product: " + ", ".join(f"{s}={product_counts.get(s, 0):,}" for s in product_order))
    logger.info(f"  Wrote {table_path}, {summary_path} and {review_path} "
                f"({len(review_rows):,} genes to review)")

    return {"symbol": dict(symbol_counts), "product": dict(product_counts),
            "n_existing_symbol": n_existing_symbol, "n_existing_product": n_existing_product,
            "n_existing_identifier": n_existing_identifier}
