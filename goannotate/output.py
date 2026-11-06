"""Output writers: mapping files, evidence, annotated FASTA and GFF3, summary."""

import logging
import re
import statistics
import sys
import urllib.parse
from datetime import datetime
from pathlib import Path

from goannotate import __version__
from goannotate.gff3 import encode_attribute_value, is_gene_type
from goannotate.models import BlastStats, GeneUniverse

logger = logging.getLogger(__name__)


####################
# PRODUCT NAMES
####################

def strip_uniprot_prefix(protein_name: str) -> str:
    """Remove a leading sp|ACC|ENTRY or tr|ACC|ENTRY token from a protein name."""
    if not protein_name:
        return protein_name
    match = re.match(r'^(?:sp|tr)\|[A-Za-z0-9]+\|[A-Za-z0-9_]+\s+(.+)$', protein_name)
    if match:
        return match.group(1)
    return protein_name


_GREEK = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "ζ": "zeta",
    "η": "eta", "θ": "theta", "ι": "iota", "κ": "kappa", "λ": "lambda", "μ": "mu",
    "ν": "nu", "ξ": "xi", "π": "pi", "ρ": "rho", "σ": "sigma", "τ": "tau",
    "υ": "upsilon", "φ": "phi", "χ": "chi", "ψ": "psi", "ω": "omega",
    "Α": "alpha", "Β": "beta", "Γ": "gamma", "Δ": "delta", "Ω": "omega",
    "µ": "mu", "’": "'", "‘": "'", "“": '"', "”": '"', "–": "-", "—": "-",
}


_product_cleaner = None      # BadNameFilter whose product_regex rules are applied to products


def configure_products(bad_name_filter, cleanup: bool = True) -> None:
    """
    Set the product naming rules (product_regex lines of the pattern file)
    that sanitize_product applies, or none with cleanup False.
    """
    global _product_cleaner
    _product_cleaner = bad_name_filter if cleanup else None


def product_cleanup_audit(annotations: dict) -> dict:
    """Rule text -> number of products changed by it, over the final products of all annotations."""
    audit = {}
    if _product_cleaner is None:
        return audit
    for annot in annotations.values():
        if annot.consensus_name:
            _product_cleaner.clean_product(strip_uniprot_prefix(annot.consensus_name), audit)
    return audit


def sanitize_product(name: str, ncbi_case: bool = True) -> str:
    """
    Make a product name acceptable for GFF3 output and NCBI submission:
    database prefixes, bracketed qualifiers, evidence codes, unbalanced
    parentheses, non-ASCII characters, GFF3 separators and trailing
    punctuation are removed; whitespace is collapsed. The product_regex
    rules of the pattern file (configure_products) are applied next: EC
    numbers, molecular weights, species names, "<SYMBOL> protein" word
    order. With ncbi_case the first letter is lowercased unless the first
    word is an acronym.
    """
    if not name:
        return ""
    s = strip_uniprot_prefix(name)
    if _product_cleaner is not None:
        s = _product_cleaner.clean_product(s)
    s = re.sub(r'\s*\{ECO:[^}]*\}', '', s)
    s = re.sub(r'\s*\[[^\]]*\]', '', s)
    s = re.sub(r'\s*\{[^}]*\}', '', s)
    for ch, rep in _GREEK.items():
        s = s.replace(ch, rep)
    s = s.encode("ascii", "ignore").decode("ascii")
    if s.count("(") != s.count(")"):
        s = s.replace("(", " ").replace(")", " ")
    s = s.replace(";", ",").replace("=", " ").replace("\t", " ")
    s = re.sub(r'\s+', ' ', s).strip()
    s = re.sub(r'[\s.,;:\-]+$', '', s)
    # A trailing "-like" is the NCBI pipeline's paralog marker for one
    # genome and never carries over; internal "-like" ("Kinesin-like
    # protein KIF1A") is part of the name and stays.
    s = re.sub(r'[\s\-]+like$', '', s, flags=re.I)
    s = re.sub(r'[\s.,;:\-]+$', '', s)
    s = re.sub(r'^[\s.,;:\-]+', '', s)
    if ncbi_case and len(s) > 1 and s[0].isupper() and s[1].islower():
        first_word = s.split()[0]
        if not first_word.isupper():
            s = s[0].lower() + s[1:]
    return s


_missing_symbol_check = None    # existing symbols that count as missing in fill mode (LOC numbers, locus tags)
_missing_product_check = None   # existing products that count as missing ("uncharacterized protein LOC...")


def configure_existing(bad_name_filter, keep_identifiers: bool = False) -> None:
    """
    In fill mode an existing symbol or product is kept and GOAnnotate fills
    only what is missing. Existing symbols that are genome-specific
    identifiers (LOC numbers of an NCBI annotation, locus tags, clone IDs;
    the symbol patterns) and existing products that carry no information
    ("uncharacterized protein LOC122519524", "hypothetical protein") count
    as missing and are filled, unless keep_identifiers is set.
    """
    global _missing_symbol_check, _missing_product_check
    if keep_identifiers or bad_name_filter is None:
        _missing_symbol_check = _missing_product_check = None
        return
    _missing_symbol_check = bad_name_filter.is_bad_symbol
    _missing_product_check = lambda p: bad_name_filter.process_name(p) is None  # noqa: E731


def missing_symbol(value: str) -> bool:
    return not value or bool(_missing_symbol_check and _missing_symbol_check(value))


def missing_product(value: str) -> bool:
    return not value or bool(_missing_product_check and _missing_product_check(value))


def product_for(annot) -> str:
    return sanitize_product(annot.consensus_name) if annot and annot.consensus_name else ""


def top_hit_accession(annot) -> str:
    """UniProt accession of the best hit in the winning symbol group, if any."""
    if not annot or not annot.filtered_hits:
        return ""
    sid = annot.filtered_hits[0].subject_id
    parts = sid.split("|")
    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
        return parts[1]
    return ""


####################
# MAPPING FILES
####################

def write_mapping_files(annotations: dict, output_dir: Path, prefix: str,
                        universe: GeneUniverse = None, evidence_dir: Path = None) -> None:
    """gene2symbol, gene2name, gene2go, and the annotation evidence table (in evidence_dir)."""
    symbol_file = output_dir / f"{prefix}_gene2symbol.tsv"
    name_file = output_dir / f"{prefix}_gene2name.tsv"
    go_file = output_dir / f"{prefix}_gene2go.tsv"
    evidence_file = Path(evidence_dir or output_dir) / f"{prefix}_annotation_evidence.tsv"

    logger.info(f"Writing gene2symbol to {symbol_file}")
    with open(symbol_file, 'w') as f:
        f.write("gene_id\tgene_symbol\n")
        for gene_id, annot in sorted(annotations.items()):
            f.write(f"{gene_id}\t{annot.consensus_symbol}\n")

    logger.info(f"Writing gene2name to {name_file}")
    with open(name_file, 'w') as f:
        f.write("gene_id\tprotein_name\n")
        for gene_id, annot in sorted(annotations.items()):
            f.write(f"{gene_id}\t{product_for(annot)}\n")

    logger.info(f"Writing gene2go to {go_file}")
    with open(go_file, 'w') as f:
        f.write("gene_id\tgo_terms\n")
        for gene_id, annot in sorted(annotations.items()):
            go_str = ';'.join(sorted(annot.specific_go_terms)) if annot.specific_go_terms else ""
            f.write(f"{gene_id}\t{go_str}\n")

    logger.info(f"Writing annotation evidence to {evidence_file}")
    with open(evidence_file, 'w') as f:
        f.write("gene_id\tgene_symbol\tprotein_name\tn_isoforms\tn_total_hits\t"
                "n_cluster_hits\twinning_cluster_names\twinning_cluster_symbols\t"
                "top_hit_sseqid\ttop_hit_bitscore\ttop_hit_source\t"
                "mean_cluster_bitscore\tname_concordance\tsymbol_concordance\t"
                "go_source\tn_go_terms\tconfidence\tconfidence_score\tgene_type\tcoding\tsymbol_status\t"
                "voted_symbol\tcoverage_status\tbest_query_coverage\tbest_subject_coverage\t"
                "go_evidence_classes\tn_go_taxon_removed\tproduct_note\tsymbol_aliases\tsecondary_symbol\n")
        for gene_id, annot in sorted(annotations.items()):
            cluster_names = "|".join(annot.winning_cluster_names)
            cluster_symbols = "|".join(annot.winning_cluster_symbols)
            if annot.filtered_hits:
                top_hit = annot.filtered_hits[0]
                top_sseqid = top_hit.subject_id
                top_bitscore = f"{top_hit.bitscore:.2f}"
                if top_hit.subject_id.startswith('sp|'):
                    top_source = "sp"
                elif top_hit.subject_id.startswith('tr|'):
                    top_source = "tr"
                else:
                    top_source = ""
            else:
                top_sseqid = ""
                top_bitscore = "0.00"
                top_source = ""

            gene_type = universe.gene_types.get(gene_id, "") if universe else ""
            coding = universe.coding.get(gene_id) if universe else None
            coding_str = "" if coding is None else ("yes" if coding else "no")

            f.write(f"{gene_id}\t{annot.consensus_symbol}\t{product_for(annot)}\t"
                    f"{annot.n_isoforms}\t{len(annot.hits)}\t"
                    f"{annot.n_cluster_hits}\t{cluster_names}\t{cluster_symbols}\t"
                    f"{top_sseqid}\t{top_bitscore}\t{top_source}\t"
                    f"{annot.mean_cluster_bitscore:.2f}\t"
                    f"{annot.name_concordance:.2f}\t{annot.symbol_concordance:.2f}\t"
                    f"{annot.go_source}\t{len(annot.specific_go_terms)}\t"
                    f"{annot.confidence}\t{annot.confidence_score:.3f}\t{gene_type}\t{coding_str}\t"
                    f"{annot.symbol_status}\t{annot.voted_symbol}\t{annot.coverage_status}\t"
                    f"{annot.best_query_coverage:.2f}\t{annot.best_subject_coverage:.2f}\t"
                    f"{go_evidence_summary(annot)}\t{len(annot.go_taxon_removed)}\t{annot.product_note}\t"
                    f"{','.join(annot.symbol_aliases)}\t{annot.secondary_symbol}\n")


def go_evidence_summary(annot) -> str:
    """Counts of the gene's specific GO terms by strongest evidence class, e.g. experimental:12,electronic:30."""
    if not annot.specific_go_terms:
        return ""
    from goannotate.go import best_evidence_class
    counts = {}
    for t in annot.specific_go_terms:
        codes = annot.go_evidence.get(t, "")
        cls = best_evidence_class(set(codes.split(",")) if codes else set())
        counts[cls] = counts.get(cls, 0) + 1
    return ",".join(f"{k}:{v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1]))


def write_go_evidence(annotations: dict, go_hierarchy, path: Path, forbidden_go: dict = None) -> None:
    """
    One row per gene and GO term (the specific terms written to gene2go
    and the tables) with the term's name, namespace, the evidence codes
    behind it among the source hits, the strongest evidence class and the
    number of source hits supporting it; terms removed for weak support
    follow with status removed_support and their source (hits or reference), and terms removed
    by the taxon constraints with status removed_taxon and the constraint
    that removed them.
    """
    from goannotate.go import best_evidence_class
    forbidden_go = forbidden_go or {}
    logger.info(f"Writing GO evidence to {path}")
    n_rows = n_removed = 0
    with open(path, "w") as f:
        f.write("gene_id\tgo_id\tname\tnamespace\tstatus\tevidence_codes\tevidence_class\tn_hits\tsource\t"
                "interpro_entries\treference_proteins\tconstraint\n")
        for gene_id, annot in sorted(annotations.items()):
            for t in sorted(annot.specific_go_terms):
                codes = annot.go_evidence.get(t, "")
                cls = best_evidence_class(set(codes.split(",")) if codes else set())
                from_ipr = annot.go_interpro.get(t)
                from_ref = annot.go_reference.get(t)
                from_hits = annot.go_support.get(t, 0) > 0
                source = "+".join(s for s, ok in (("hits", from_hits), ("reference", bool(from_ref)),
                                                  ("interpro", bool(from_ipr))) if ok)
                f.write(f"{gene_id}\t{t}\t{go_hierarchy.names.get(t, '')}\t{go_hierarchy.namespace.get(t, '')}\t"
                        f"specific\t{codes}\t{cls}\t{annot.go_support.get(t, '') if from_hits else ''}\t{source}\t"
                        f"{','.join(from_ipr) if from_ipr else ''}\t{','.join(from_ref) if from_ref else ''}\t\n")
                n_rows += 1
            for t, (support, codes, source) in sorted(annot.go_weak_removed.items()):
                cls = best_evidence_class(set(codes.split(",")) if codes else set())
                refs = ",".join(annot.go_reference.get(t) or []) if source == "reference" else ""
                f.write(f"{gene_id}\t{t}\t{go_hierarchy.names.get(t, '')}\t{go_hierarchy.namespace.get(t, '')}\t"
                        f"removed_support\t{codes}\t{cls}\t{support}\t{source}\t\t{refs}\t\n")
                n_removed += 1
            for t in sorted(annot.go_taxon_removed):
                f.write(f"{gene_id}\t{t}\t{go_hierarchy.names.get(t, '')}\t{go_hierarchy.namespace.get(t, '')}\t"
                        f"removed_taxon\t\t\t\t\t\t\t{forbidden_go.get(t, '')}\n")
                n_removed += 1
    logger.info(f"  {n_rows:,} gene-term rows" + (f", {n_removed:,} removed (weak support, taxon constraints)" if n_removed else ""))


####################
# ANNOTATION TABLES
####################

_ISOFORM_SUFFIX = re.compile(r",?\s+(?:isoform|transcript variant)\s+\S+$", re.I)


def _choose(input_value: str, goa_value: str, mode: str, kind: str = "symbol") -> tuple:
    """
    (value, source) under the annotation mode: replace prefers GOAnnotate;
    fill and add keep the input unless it counts as missing (an identifier
    such as a LOC number, see configure_existing), in which case
    GOAnnotate's value is used when there is one.
    """
    if mode == "replace":
        if goa_value:
            return goa_value, "goannotate"
        return (input_value, "input") if input_value else ("", "")
    missing = missing_symbol if kind == "symbol" else missing_product
    if input_value and not (goa_value and missing(input_value)):
        return input_value, "input"
    return (goa_value, "goannotate") if goa_value else ("", "")


def _tx_input_product(tx) -> str:
    return tx.cds_attrs.get("product") or tx.attrs.get("product", "")


def _gene_input_product(gene) -> str:
    """Gene description, else the product of the first coding transcript without its isoform label."""
    desc = gene.attrs.get("description", "")
    if desc:
        return desc
    for coding_only in (True, False):
        for tx in gene.transcripts:
            if coding_only and not tx.coding:
                continue
            p = _tx_input_product(tx)
            if p:
                return _ISOFORM_SUFFIX.sub("", p)
    return ""


def _span(seqid, start, end, strand) -> str:
    return f"{seqid}:{start}-{end}({strand})" if seqid else ""


def write_annotation_tables(output_dir: Path, prefix: str, annotations: dict, universe: GeneUniverse,
                            genes: list = None, fasta_headers: str = None, mode: str = "fill",
                            protein_lengths: dict = None, modules: list = None,
                            extra_columns: tuple = None) -> tuple:
    """
    Two tab-delimited tables for everyday lookup, one row per gene and one
    row per transcript, with the identifiers in the first column(s) and the
    specific GO terms (as in gene2go) in the last column:

      {prefix}_gene_annotations.tsv
        gene_id, symbol, symbol_aliases, product, gene_type, coding, locus,
        n_transcripts, transcripts, confidence, symbol_source, product_source,
        uniprot_accession, go_terms
      {prefix}_transcript_annotations.tsv
        gene_id, transcript_id, symbol, symbol_aliases, product, transcript_type, coding,
        locus, cds_locus, protein_length, confidence, symbol_source,
        product_source, uniprot_accession, go_terms

    symbol and product follow --annotation-mode, as in the annotated GFF3
    and FASTA: with fill and add the input annotation is kept where it
    exists and GOAnnotate's result is used where it does not; with replace
    GOAnnotate's result is used wherever there is one. symbol_source and
    product_source say which ('input' or 'goannotate'). Transcript rows
    keep the input transcript's own product (with its isoform label) under
    fill and add. Loci are 1-based and inclusive; protein_length comes from
    the protein FASTA that was searched (protein_lengths).

    Without a GFF3 (transcript FASTA input), the rows come from the gene
    universe and the input values from the FASTA header fields (symbol=,
    product=, type=, locus=) where present.

    modules: optional modules that ran; their gene_columns() and
    transcript_columns() are appended after go_terms.
    extra_columns: (gene columns, transcript columns), each a list of
    (name, function of the ID), appended after go_terms before the
    modules' columns (the completeness of translated transcripts).
    """
    gene_extra, tx_extra = (list(extra_columns[0]), list(extra_columns[1])) if extra_columns else ([], [])
    gene_extra += [col for m in (modules or []) for col in m.gene_columns()]
    tx_extra += [col for m in (modules or []) for col in m.transcript_columns()]

    def extra(columns, key):
        return [fn(key) for _, fn in columns]
    gene_file = output_dir / f"{prefix}_gene_annotations.tsv"
    tx_file = output_dir / f"{prefix}_transcript_annotations.tsv"
    logger.info(f"Writing {gene_file.name} and {tx_file.name}")

    def aliases(a, symbol):
        """GOAnnotate's aliases, when the output symbol is GOAnnotate's symbol (any case)."""
        if a is None or not a.symbol_aliases or not symbol or symbol.lower() != a.consensus_symbol.lower():
            return ""
        return ",".join(a.symbol_aliases)

    def goa(gene_key):
        a = annotations.get(gene_key)
        if a is None:
            return None, "", "", "", ""
        go = ";".join(sorted(a.specific_go_terms)) if a.specific_go_terms else ""
        conf = a.confidence if a.hits else ""
        return a, a.consensus_symbol, product_for(a), conf, go

    gene_cols = ["gene_id", "symbol", "symbol_aliases", "secondary_symbol", "product", "gene_type", "coding", "locus", "n_transcripts", "transcripts",
                 "confidence", "symbol_source", "product_source", "uniprot_accession", "go_terms"]
    tx_cols = ["gene_id", "transcript_id", "symbol", "symbol_aliases", "secondary_symbol", "product", "transcript_type", "coding", "locus",
               "cds_locus", "protein_length", "confidence", "symbol_source", "product_source",
               "uniprot_accession", "go_terms"]
    gene_cols += [name for name, _ in gene_extra]
    tx_cols += [name for name, _ in tx_extra]
    n_genes = n_tx = 0

    with open(gene_file, "w") as fg, open(tx_file, "w") as ft:
        fg.write("\t".join(gene_cols) + "\n")
        ft.write("\t".join(tx_cols) + "\n")

        if genes is not None:
            from goannotate.sequences import gene_symbol_of, transcript_type_of
            for gene in genes:
                a, g_sym, g_prod, conf, go = goa(gene.key)
                acc = top_hit_accession(a)
                symbol, sym_src = _choose(gene_symbol_of(gene), g_sym, mode)
                product, prod_src = _choose(_gene_input_product(gene), g_prod, mode, "product")
                txs = gene.transcripts
                al = aliases(a, symbol)
                sec = a.secondary_symbol if a is not None else ""
                fg.write("\t".join([
                    gene.key, symbol, al, sec, product, gene.biotype, "yes" if gene.coding else "no",
                    _span(gene.seqid, gene.start, gene.end, gene.strand), str(len(txs)),
                    ",".join(t.key for t in txs), conf, sym_src, prod_src, acc, go]
                    + extra(gene_extra, gene.key)) + "\n")
                n_genes += 1
                for tx in txs:
                    t_prod, t_src = _choose(_tx_input_product(tx), g_prod, mode, "product")
                    if tx.exons:
                        start, end = min(s for s, _ in tx.exons), max(e for _, e in tx.exons)
                    else:
                        start, end = tx.start, tx.end
                    cds_locus = protein_length = ""
                    if tx.coding:
                        cds_locus = _span(tx.seqid, min(s for s, _e, _p in tx.cds),
                                          max(e for _s, e, _p in tx.cds), tx.strand)
                        if protein_lengths and tx.key in protein_lengths:
                            protein_length = str(protein_lengths[tx.key])
                    ft.write("\t".join([
                        gene.key, tx.key, symbol, al, sec, t_prod, transcript_type_of(tx),
                        "yes" if tx.coding else "no", _span(tx.seqid, start, end, tx.strand), cds_locus,
                        protein_length, conf, sym_src, t_src, acc, go]
                        + extra(tx_extra, tx.key)) + "\n")
                    n_tx += 1
        else:
            from goannotate.sequences import parse_header
            fields = {}
            if fasta_headers:
                with open(fasta_headers) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            seq_id, f = parse_header(line.rstrip("\n"))
                            fields[seq_id] = f
            for gene_id in universe.gene_ids:
                a, g_sym, g_prod, conf, go = goa(gene_id)
                acc = top_hit_accession(a)
                txs = universe.gene_to_transcripts.get(gene_id, [])
                first = fields.get(txs[0], {}) if txs else {}
                symbol, sym_src = _choose(first.get("symbol", ""), g_sym, mode)
                in_prod = _ISOFORM_SUFFIX.sub("", first.get("product", ""))
                product, prod_src = _choose(in_prod, g_prod, mode, "product")
                coding = universe.coding.get(gene_id) if universe.coding else None
                al = aliases(a, symbol)
                sec = a.secondary_symbol if a is not None else ""
                fg.write("\t".join([
                    gene_id, symbol, al, sec, product, universe.gene_types.get(gene_id, "") if universe.gene_types else "",
                    "" if coding is None else ("yes" if coding else "no"), "", str(len(txs)), ",".join(txs),
                    conf, sym_src, prod_src, acc, go] + extra(gene_extra, gene_id)) + "\n")
                n_genes += 1
                for t in txs:
                    f = fields.get(t, {})
                    t_prod, t_src = _choose(f.get("product", ""), g_prod, mode, "product")
                    ft.write("\t".join([
                        gene_id, t, symbol, al, sec, t_prod, f.get("type", ""), "", f.get("locus", ""), "", "",
                        conf, sym_src, t_src, acc, go] + extra(tx_extra, t)) + "\n")
                    n_tx += 1

    logger.info(f"  {n_genes:,} genes, {n_tx:,} transcripts")
    return gene_file, tx_file


####################
# ANNOTATED FASTA
####################

def write_annotated_fasta(input_fasta: str, annotations: dict, output_file: Path,
                          transcript_to_gene: dict = None, mode: str = "fill",
                          header_fields: dict = None) -> None:
    """
    Copy a FASTA, adding GOAnnotate results to the headers of annotated
    genes in the GRAFT header format (gene= symbol= product="..." type=
    locus=). mode follows --annotation-mode:
      add      goannotate_symbol= and goannotate_product="..." are appended
      fill     symbol= and product= are set where the header has none
      replace  symbol= and product= are set to the GOAnnotate values
    A legacy Name= field (GOAnnotate 1.x) is read as symbol=.
    header_fields (sequence ID -> fields) adds further fields, such as the
    ORF coordinates of translated transcripts, to the headers that have them.
    """
    from goannotate.sequences import format_header, parse_header
    logger.info(f"Writing annotated FASTA to {output_file} (mode: {mode})")

    annotated_count = 0
    total_count = 0

    with open(input_fasta) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if not line.startswith('>'):
                fout.write(line)
                continue
            total_count += 1
            seq_id, fields = parse_header(line.rstrip("\n"))
            extra = header_fields.get(seq_id) if header_fields else None
            if extra:
                fields.update(extra)
            gene_id = transcript_to_gene.get(seq_id, seq_id) if transcript_to_gene else seq_id
            annot = annotations.get(gene_id)
            if not annot or not (annot.consensus_symbol or annot.consensus_name):
                if extra:
                    fout.write(">" + format_header(seq_id, fields) + "\n")
                else:
                    fout.write(line)
                continue
            annotated_count += 1
            if "Name" in fields and "symbol" not in fields:
                fields["symbol"] = fields.pop("Name")
            symbol = annot.consensus_symbol
            product = product_for(annot)
            if mode == "add":
                fields["goannotate_symbol"] = symbol
                fields["goannotate_product"] = product
            else:
                overwrite = (mode == "replace")
                if symbol and (overwrite or missing_symbol(fields.get("symbol", ""))):
                    fields["symbol"] = symbol
                if product and (overwrite or missing_product(fields.get("product", ""))):
                    fields["product"] = product
            fout.write(">" + format_header(seq_id, fields) + "\n")

    logger.info(f"Annotated {annotated_count:,} of {total_count:,} sequences in FASTA")


####################
# ANNOTATED GFF3
####################

def _split_raw_attributes(text: str) -> list:
    """Split column 9 into (key, raw_value) pairs without decoding."""
    pairs = []
    for token in text.strip().strip(";").split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            pairs.append([k.strip(), v.strip()])
        else:
            pairs.append([token, None])
    return pairs


def _join_raw_attributes(pairs: list) -> str:
    return ";".join(k if v is None else f"{k}={v}" for k, v in pairs)


def _get(pairs: list, key: str):
    for k, v in pairs:
        if k == key:
            return v
    return None


def _set(pairs: list, key: str, value: str, overwrite: bool) -> None:
    """Add key=value (encoded) or replace the existing value when overwrite is set."""
    encoded = encode_attribute_value(value)
    for pair in pairs:
        if pair[0] == key:
            if overwrite:
                pair[1] = encoded
            return
    pairs.append([key, encoded])


def _merge_list(pairs: list, key: str, values: list) -> None:
    """Merge comma-separated values into an existing list attribute (Dbxref, Ontology_term)."""
    if not values:
        return
    encoded = [encode_attribute_value(v) for v in values]
    for pair in pairs:
        if pair[0] == key:
            existing = [x for x in (pair[1] or "").split(",") if x]
            for e in encoded:
                if e not in existing:
                    existing.append(e)
            pair[1] = ",".join(existing)
            return
    pairs.append([key, ",".join(encoded)])


# GFF3 attributes holding comma-separated lists: module values are merged into them, not overwritten
_LIST_ATTRIBUTES = {"Dbxref", "Ontology_term"}


def write_annotated_gff3(input_gff: str, output_file: Path, annotations: dict, genes: list,
                         mode: str = "fill", hypothetical: bool = True, modules: list = None) -> dict:
    """
    Stream the input GFF3 to output_file, passing every line through and
    adding GOAnnotate results to gene, transcript and CDS features.

    Always added where a gene has results: goannotate_symbol,
    goannotate_product and goannotate_confidence on the gene feature,
    Ontology_term (specific GO terms) and a UniProtKB Dbxref on the gene.

    mode controls the standard attributes:
      add      nothing else is changed
      fill     gene, Name and description on the gene, product on
               transcripts and CDS, written only where absent
      replace  the same attributes, overwritten with GOAnnotate values

    With hypothetical, coding transcripts and CDS that end up without any
    product receive product=hypothetical protein (NCBI convention) in the
    fill and replace modes.

    modules: optional modules that ran. Their gene_attributes() and
    transcript_attributes() are added to the gene and transcript features
    and their transcript_features() lines are written after each
    transcript's line.
    """
    logger.info(f"Writing annotated GFF3 to {output_file} (mode: {mode})")

    overwrite = (mode == "replace")
    standard = mode in ("fill", "replace")
    modules = modules or []
    tx_objects = {}        # transcript feature ID -> (Transcript, Gene), for module features

    # Decide once per gene and per transcript which standard attributes
    # receive GOAnnotate values, from the parsed (decoded) input attributes,
    # so that the decision does not depend on line order in the file.
    gene_key_by_feature = {}
    gene_decision = {}    # gene key -> dict
    tx_decision = {}      # transcript feature ID -> dict
    for g in genes:
        annot = annotations.get(g.key)
        symbol = annot.consensus_symbol if annot else ""
        product = product_for(annot)
        gene_decision[g.key] = {
            "symbol": symbol,
            "product": product,
            "write_gene": bool(standard and symbol and (overwrite or missing_symbol(g.attrs.get("gene", "")))),
            "write_name": bool(standard and symbol and (overwrite or missing_symbol(g.attrs.get("Name", "")))),
            "write_description": bool(standard and product
                                      and (overwrite or missing_product(g.attrs.get("description", "")))),
        }
        if not g.wrapper:
            gene_key_by_feature[g.id] = g.key
        for t in g.transcripts:
            if t.implicit:
                continue
            has_product = not missing_product(t.attrs.get("product", ""))
            tx_objects[t.id] = (t, g)
            tx_decision[t.id] = {
                "gene_key": g.key,
                "tx_key": t.key,
                "wrapper": g.wrapper,
                "write_product": bool(standard and product and (overwrite or not has_product)),
                "write_hypothetical": bool(standard and hypothetical and t.coding
                                           and not product and not has_product),
            }

    counts = {"genes_annotated": 0, "transcripts_annotated": 0, "hypothetical": 0}

    def gene_attrs(pairs, gene_key):
        annot = annotations.get(gene_key)
        d = gene_decision[gene_key]
        if d["symbol"]:
            _set(pairs, "goannotate_symbol", d["symbol"], True)
            if annot is not None and annot.symbol_aliases:
                _merge_list(pairs, "goannotate_symbol_aliases", annot.symbol_aliases)
        if annot is not None and annot.secondary_symbol:
            _set(pairs, "goannotate_secondary_symbol", annot.secondary_symbol, True)
        if d["product"]:
            _set(pairs, "goannotate_product", d["product"], True)
        if annot and annot.hits:
            _set(pairs, "goannotate_confidence", annot.confidence, True)
        if annot and annot.specific_go_terms:
            _merge_list(pairs, "Ontology_term", sorted(annot.specific_go_terms))
        acc = top_hit_accession(annot)
        if acc:
            _merge_list(pairs, "Dbxref", [f"UniProtKB:{acc}"])
        if d["write_gene"]:
            _set(pairs, "gene", d["symbol"], True)
        if d["write_name"]:
            _set(pairs, "Name", d["symbol"], True)
        if d["write_description"]:
            _set(pairs, "description", d["product"], True)
        if d["symbol"] or d["product"]:
            counts["genes_annotated"] += 1
        for m in modules:
            for k, v in m.gene_attributes(gene_key).items():
                if k in _LIST_ATTRIBUTES:
                    _merge_list(pairs, k, [x for x in v.split(",") if x])
                else:
                    _set(pairs, k, v, True)

    def transcript_attrs(pairs, tx_id, feature_type):
        t = tx_decision[tx_id]
        g = gene_decision[t["gene_key"]]
        if t["write_product"]:
            _set(pairs, "product", g["product"], True)
            if feature_type != "CDS":
                counts["transcripts_annotated"] += 1
        elif t["write_hypothetical"]:
            _set(pairs, "product", "hypothetical protein", False)
            if feature_type != "CDS":
                counts["hypothetical"] += 1
        if g["write_gene"]:
            _set(pairs, "gene", g["symbol"], True)
        if feature_type != "CDS":
            for m in modules:
                for k, v in m.transcript_attributes(t["tx_key"]).items():
                    if k in _LIST_ATTRIBUTES:
                        _merge_list(pairs, k, [x for x in v.split(",") if x])
                    else:
                        _set(pairs, k, v, True)

    def module_features(tx_id) -> list:
        if not modules or tx_id not in tx_objects:
            return []
        t, g = tx_objects[tx_id]
        lines = []
        for m in modules:
            lines.extend(m.transcript_features(t, g))
        return lines

    with open(input_gff) as fin, open(output_file, "w") as fout:
        in_fasta = False
        for line in fin:
            if in_fasta or line.startswith("#"):
                if line.startswith("##FASTA"):
                    in_fasta = True
                fout.write(line)
                continue
            stripped = line.rstrip("\n")
            if not stripped.strip():
                fout.write(line)
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                fout.write(line)
                continue

            ftype = fields[2]
            pairs = _split_raw_attributes(fields[8])
            fid = urllib.parse.unquote(_get(pairs, "ID") or "")
            parent = urllib.parse.unquote(_get(pairs, "Parent") or "")

            extra_lines = []
            if is_gene_type(ftype) and fid in gene_key_by_feature:
                gene_attrs(pairs, gene_key_by_feature[fid])
            elif fid in tx_decision:
                if tx_decision[fid]["wrapper"]:
                    gene_attrs(pairs, tx_decision[fid]["gene_key"])
                transcript_attrs(pairs, fid, ftype)
                extra_lines = module_features(fid)
            elif ftype == "CDS":
                for p in parent.split(","):
                    if p in tx_decision:
                        transcript_attrs(pairs, p, "CDS")
                        break
            else:
                fout.write(line)
                continue

            fields[8] = _join_raw_attributes(pairs)
            fout.write("\t".join(fields) + "\n")
            for extra_line in extra_lines:
                fout.write(extra_line + "\n")
                counts["module_features"] = counts.get("module_features", 0) + 1

    logger.info(f"  Genes with GOAnnotate symbol or product: {counts['genes_annotated']:,}; "
                f"transcripts given a product: {counts['transcripts_annotated']:,}; "
                f"coding transcripts set to hypothetical protein: {counts['hypothetical']:,}"
                + (f"; module features written: {counts['module_features']:,}" if counts.get("module_features") else ""))
    return counts


####################
# SUMMARY
####################

def write_summary(output_dir: Path, prefix: str, annotations: dict, args, blast_stats: BlastStats = None,
                  universe: GeneUniverse = None, extraction=None, gff_stats=None, translation=None,
                  comparison: dict = None,
                  search_results: str = None, assignment: dict = None, protein_check=None,
                  coverage=None, forbidden_go: dict = None, product_audit: dict = None,
                  word_weights=None, modules: list = None, reference_go: dict = None,
                  go_weak: dict = None) -> None:
    stats_file = output_dir / f"{prefix}_summary.txt"
    logger.info(f"Writing summary to {stats_file}")

    total = len(annotations)
    with_blast_hits = sum(1 for a in annotations.values() if a.hits)
    with_protein_name = sum(1 for a in annotations.values() if a.consensus_name)
    with_gene_symbol = sum(1 for a in annotations.values() if a.consensus_symbol)
    with_go = sum(1 for a in annotations.values() if a.go_terms)

    all_go = set()
    all_specific = set()
    for annot in annotations.values():
        all_go.update(annot.go_terms)
        all_specific.update(annot.specific_go_terms)

    confidence_counts = {}
    for annot in annotations.values():
        confidence_counts[annot.confidence] = confidence_counts.get(annot.confidence, 0) + 1

    # Coverage is reported against genes that can have hits (coding genes
    # of a GFF3; every gene of a transcript set)
    if universe is not None and universe.coding:
        n_coding = sum(1 for v in universe.coding.values() if v)
        pct_label = "of coding"
    else:
        n_coding = total
        pct_label = "of genes"

    def pct(n):
        return f"{100 * n / n_coding:.1f}% {pct_label}" if n_coding else "0.0%"

    with open(stats_file, 'w') as f:
        f.write(f"GOAnnotate {__version__} Summary\n")
        f.write("=" * 70 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("Command:\n")
        f.write("-" * 70 + "\n")
        f.write("  " + " ".join(sys.argv) + "\n\n")

        f.write("Inputs:\n")
        f.write("-" * 70 + "\n")
        if getattr(args, "gff3", None):
            f.write(f"  GFF3:           {args.gff3}\n")
        if getattr(args, "genome", None):
            f.write(f"  Genome:         {args.genome}\n")
        if getattr(args, "proteins", None):
            f.write(f"  Proteins:       {args.proteins}\n")
        if getattr(args, "transcripts", None):
            f.write(f"  Transcripts:    {args.transcripts}\n")
        if search_results:
            f.write(f"  BLAST results:  {search_results}\n")
        if getattr(args, "db", None):
            f.write(f"  SQLite DB:      {args.db}\n")
        if getattr(args, "go_mapping", None):
            f.write(f"  GO mapping:     {args.go_mapping}\n")
        f.write(f"  GO OBO:         {args.go_obo}\n")
        f.write(f"  Bad names:      {args.bad_names}\n\n")

        f.write("Parameters:\n")
        f.write("-" * 70 + "\n")
        f.write(f"  E-value threshold:        {args.evalue}\n")
        f.write(f"  Top N hits:               {args.top_n}\n")
        f.write(f"  Hit coverage filter:      {coverage.describe() if coverage else 'off'}\n")
        f.write(f"  Consensus threshold:      {args.consensus_threshold}\n")
        f.write(f"  Word weighting:           {word_weights.summary() if word_weights else 'off'}\n")
        f.write(f"  GO namespaces:            {', '.join(args.namespace)}\n")
        f.write(f"  GO evidence policy:       {getattr(args, 'go_evidence', 'all')}\n")
        f.write(f"  GO terms from hits:       top {getattr(args, 'go_source_hits', 10)} hits of the winning symbol "
                f"group plus SwissProt hits of the gene; a term from the hits alone needs "
                f"{getattr(args, 'go_min_hits', 5)} supporting hits, or a majority of fewer source hits (at least 2; "
                f"{getattr(args, 'go_iea_min_hits', 2)} if electronic only; one hit covering "
                f"{getattr(args, 'go_single_hit_coverage', 0.7):.0%} of both proteins at "
                f"{getattr(args, 'go_single_hit_identity', 50.0):.0f}% identity keeps its non-electronic terms) "
                f"unless experimental, curated in a SwissProt entry of the gene, or confirmed by the "
                f"reference gene or InterPro; a removed term is replaced by its supported ancestors\n")
        if getattr(args, "taxid", None):
            f.write(f"  Query taxon:              {args.taxid}; GO taxon constraints "
                    f"{'exclude ' + format(len(forbidden_go), ',') + ' terms' if forbidden_go else 'not applied'}\n")
        else:
            f.write("  GO taxon constraints:     not applied (no --taxid)\n")
        f.write(f"  Reference species:        {getattr(args, 'reference_taxid', '')} "
                f"(ortholog anchor {getattr(args, 'ortholog_anchor_taxid', '')})\n")
        f.write(f"  Secondary nomenclature:   {getattr(args, 'secondary_reference_taxid', None) or 'none'}\n")
        f.write(f"  Symbol policy:            {getattr(args, 'symbol_policy', 'reference')}\n")
        f.write(f"  Product naming rules:     {'off' if getattr(args, 'no_product_cleanup', False) else 'on'}\n")
        if getattr(args, "gff3", None):
            f.write(f"  GFF3 annotation mode:     {args.annotation_mode}\n")
        f.write(f"  Modules:                  {getattr(args, 'modules', None) or 'none'}\n")
        f.write("\n")

        if gff_stats is not None:
            f.write("GFF3 Input:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Genes:                        {gff_stats.n_genes:,} "
                    f"({gff_stats.n_coding_genes:,} coding)\n")
            f.write(f"  Transcripts:                  {gff_stats.n_transcripts:,} "
                    f"({gff_stats.n_coding_transcripts:,} coding)\n")
            if gff_stats.n_wrapper_genes:
                f.write(f"  Parentless transcripts:       {gff_stats.n_wrapper_genes:,}\n")
            if gff_stats.n_implicit_transcripts:
                f.write(f"  Implicit transcripts:         {gff_stats.n_implicit_transcripts:,}\n")
            f.write("\n")

        if extraction is not None:
            f.write("Sequence Extraction:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Transcripts:                  {extraction.n_transcripts:,}\n")
            f.write(f"  Coding transcripts:           {extraction.n_coding:,}\n")
            f.write(f"  Proteins written:             {extraction.n_proteins_written:,}\n")
            for flag, n in extraction.flag_counts.most_common():
                f.write(f"  {flag + ':':<30}{n:,}\n")
            f.write("\n")

        if translation is not None:
            f.write("ORF Translation (transcripts translated in the frame of their best blastx hit):\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Transcripts:                  {translation.n_transcripts:,}\n")
            f.write(f"  With blastx hits:             {translation.n_with_hits:,}\n")
            f.write(f"  Proteins written:             {translation.n_translated:,} "
                    f"({translation.n_short:,} shorter than 50 residues; transcripts without hits get no protein)\n")
            f.write("  Reading frames:               " + ", ".join(
                f"{status} {n:,}" for status, n in translation.completeness.most_common()) + "\n")
            if translation.notes:
                f.write("  Flagged (see the report):     " + ", ".join(
                    f"{k} {v:,}" for k, v in translation.notes.most_common()) + "\n")
            if getattr(translation, "n_blastx_fallback", 0):
                f.write(f"  Annotated from blastx hits:   {translation.n_blastx_fallback:,} "
                        f"(the protein found no blastp hit)\n")
            final = getattr(translation, "final_completeness", None)
            if final:
                f.write("  Completeness (reading frame and reference alignment; transcript tables, column "
                        "completeness):\n")
                for status in ("complete", "5_partial", "3_partial", "partial"):
                    f.write(f"    {status + ':':<28}{final.get(status, 0):,}\n")
            f.write(f"  Details:                      {translation.report}\n\n")

        if protein_check is not None:
            f.write("Protein Input Check:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Proteins in FASTA:            {protein_check.n_records:,} (all IDs match GFF3 transcripts)\n")
            f.write(f"  Coding transcripts without protein: {protein_check.n_coding_without_protein:,}\n")
            f.write(f"  Proteins for non-coding transcripts: {protein_check.n_noncoding_with_protein:,}\n")
            f.write(f"  Length differs from CDS:      {protein_check.n_length_mismatch:,}\n")
            if protein_check.n_compared:
                f.write(f"  Identical to genome extraction: {protein_check.n_identical:,} of "
                        f"{protein_check.n_compared:,}\n")
            f.write(f"  Details:                      {protein_check.report}\n\n")

        if blast_stats:
            f.write("BLAST Statistics:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Total hits in file:           {blast_stats.total_hits_in_file:,}\n")
            f.write(f"  Queries with hits (in file):  {blast_stats.queries_in_file:,}\n")
            if blast_stats.hits_filtered_by_evalue > 0:
                f.write(f"  Hits filtered by e-value:     {blast_stats.hits_filtered_by_evalue:,}\n")
            f.write(f"  Hits after e-value filter:    {blast_stats.hits_after_evalue_filter:,}\n")
            f.write(f"  Queries after e-value filter: {blast_stats.queries_after_evalue_filter:,}\n")
            queries_lost = blast_stats.queries_in_file - blast_stats.queries_after_evalue_filter
            if queries_lost > 0 and blast_stats.queries_in_file:
                f.write(f"  Queries lost to e-value:      {queries_lost:,} "
                        f"({100*queries_lost/blast_stats.queries_in_file:.1f}%)\n")
            f.write(f"  SwissProt hits:               {blast_stats.swissprot_hits:,}\n")
            f.write(f"  TrEMBL hits:                  {blast_stats.trembl_hits:,}\n")
            if blast_stats.other_hits > 0:
                f.write(f"  Other hits:                   {blast_stats.other_hits:,}\n")
            f.write("\n")

        f.write("Annotation Results:\n")
        f.write("-" * 70 + "\n")
        if universe is not None:
            f.write(f"  Total Sequences:              {universe.total_sequences:,}\n")
        f.write(f"  Total Genes:                  {total:,}\n")
        if universe is not None and universe.coding:
            f.write(f"  Coding genes:                 {n_coding:,}\n")
            f.write(f"  Non-coding genes:             {total - n_coding:,} "
                    f"(in the mapping files with empty annotations)\n")
        f.write(f"  Genes with BLAST hits:        {with_blast_hits:,} ({pct(with_blast_hits)})\n")
        f.write(f"  Genes with protein name:      {with_protein_name:,} ({pct(with_protein_name)})\n")
        f.write(f"  Genes with gene symbol:       {with_gene_symbol:,} ({pct(with_gene_symbol)})\n")
        f.write(f"  Genes with GO terms:          {with_go:,} ({pct(with_go)})\n")
        f.write("  Confidence:                   " + ", ".join(
            f"{k}={confidence_counts.get(k, 0):,}" for k in ("high", "medium", "low", "none")) + "\n")
        if coverage is not None and coverage.active:
            n_pass = sum(1 for a in annotations.values() if a.coverage_status == "pass")
            n_fb = sum(1 for a in annotations.values() if a.coverage_status == "fallback")
            f.write(f"  Hit coverage:                 {n_pass:,} genes from hits passing the filter, "
                    f"{n_fb:,} from unfiltered hits (partial alignments only)\n")
        f.write("\n")

        if assignment:
            n_distinct = len({a.consensus_symbol.upper() for a in annotations.values() if a.consensus_symbol})
            f.write("Symbol checks against the reference proteome:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Reciprocal best hits:             {assignment.get('rbh_total', 0):,}\n")
            f.write(f"  Vote confirmed by RBH:            {assignment.get('rbh', 0):,}"
                    + (f" ({assignment['rbh_same_gene_form']:,} voted in a copy form of the RBH symbol: adgrl1a for ADGRL1)"
                       if assignment.get('rbh_same_gene_form') else "") + "\n")
            f.write(f"  RBH replaced the vote:            {assignment.get('rbh_override', 0):,}\n")
            f.write(f"  RBH named an unnamed gene:        {assignment.get('rbh_added', 0):,}\n")
            f.write(f"  Vote supported by sequence:       {assignment.get('vote_supported', 0):,}\n")
            f.write(f"  Vote replaced by closest protein: {assignment.get('reference_best_hit', 0):,}\n")
            f.write(f"  Vote removed (partial homology):  {assignment.get('unsupported_removed', 0):,}\n")
            f.write(f"  Vote without reference hit, kept: {assignment.get('vote_no_reference_hit', 0):,}\n")
            if assignment.get("symbol_policy") == "any":
                f.write(f"  Non-reference symbol, kept:       {assignment.get('vote_foreign_symbol', 0):,}\n")
            elif assignment.get("symbol_policy") == "ortholog":
                f.write(f"  Non-reference symbol, replaced by the reference ortholog: "
                        f"{assignment.get('foreign_ortholog', 0):,} (suffix {assignment.get('foreign_ortholog_suffix', 0):,}, "
                        f"NCBI synonym {assignment.get('foreign_ortholog_renamed', 0):,}, "
                        f"RBH {assignment.get('foreign_ortholog_rbh', 0):,})\n")
                f.write(f"  Non-reference symbol, kept:       {assignment.get('vote_foreign_symbol', 0):,}\n")
                f.write(f"  Non-reference symbol, removed:    {assignment.get('foreign_removed', 0):,}\n")
            else:
                f.write(f"  Non-reference symbol, replaced:   {assignment.get('foreign_replaced', 0):,}\n")
                f.write(f"  Non-reference symbol, removed:    {assignment.get('foreign_removed', 0):,}\n")
            if "secondary_symbols" in assignment:
                f.write(f"  Secondary symbols ({assignment.get('secondary_reference', '')}): "
                        f"{assignment['secondary_symbols']:,} genes\n")
            if assignment.get("secondary_names"):
                f.write(f"  Products from the secondary reference: {assignment['secondary_names']:,} genes "
                        f"(the primary reference name is used by none of their hits)\n")
            if "reference_name_used" in assignment or "reference_name_not_supported" in assignment:
                f.write(f"  Reference names of confirmed symbols: {assignment.get('reference_name_used', 0):,} used, "
                        f"{assignment.get('reference_name_not_supported', 0):,} left to the hits' name (a reference "
                        f"species' own gene name that curated entries of other species do not use)\n")
            if assignment.get("established_preferred"):
                f.write(f"  Established gene preferred (PE=1 over weaker evidence): "
                        f"{assignment['established_preferred']:,}\n")
            if "designators_stripped" in assignment:
                f.write(f"  Paralog designators stripped:     {assignment['designators_stripped']:,} products\n")
            f.write(f"  Distinct symbols:                 {n_distinct:,} on {with_gene_symbol:,} genes\n\n")

        f.write("GO Term Statistics:\n")
        f.write("-" * 70 + "\n")
        f.write(f"  Total unique GO terms (all):      {len(all_go):,}\n")
        f.write(f"  Total unique GO terms (specific): {len(all_specific):,}\n")
        f.write(f"  GO term reduction: {len(all_go):,} -> {len(all_specific):,} "
                f"({100*(1-len(all_specific)/len(all_go)) if all_go else 0:.1f}% reduction)\n")
        class_counts = {}
        for annot in annotations.values():
            for cls_n in go_evidence_summary(annot).split(","):
                if ":" in cls_n:
                    k, v = cls_n.rsplit(":", 1)
                    class_counts[k] = class_counts.get(k, 0) + int(v)
        if class_counts:
            total_terms = sum(class_counts.values())
            f.write("  Specific term assignments by strongest evidence class:\n")
            for k, v in sorted(class_counts.items(), key=lambda kv: -kv[1]):
                f.write(f"    {k + ':':<16}{v:>10,}  ({100 * v / total_terms:.1f}%)\n")
        n_removed = sum(len(a.go_taxon_removed) for a in annotations.values())
        if forbidden_go:
            genes_removed = sum(1 for a in annotations.values() if a.go_taxon_removed)
            f.write(f"  Removed by taxon constraints:     {n_removed:,} assignments on {genes_removed:,} genes "
                    f"(reports/{prefix}_go_evidence.tsv)\n")
        if go_weak:
            n_ref = go_weak.get("reference_electronic", 0)
            f.write(f"  Removed for weak support:         {go_weak.get('terms', 0):,} terms on "
                    f"{go_weak.get('genes', 0):,} genes ({go_weak.get('terms', 0) - n_ref:,} from the hits alone; "
                    f"{n_ref:,} electronic-only terms of the reference gene that neither the hits nor a close "
                    f"relative's ortholog carry; {go_weak.get('specific_removed', 0):,} specific; "
                    f"{go_weak.get('ancestors_kept', 0):,} supported ancestors kept in their place; "
                    f"{go_weak.get('genes_left_without_go', 0):,} genes left without GO terms; "
                    f"reports/{prefix}_go_evidence.tsv, status removed_support)\n")
        if reference_go:
            if reference_go.get("mode") == "off":
                f.write("  From the reference gene:          off (--reference-go off)\n")
            else:
                f.write(f"  From the reference gene ({reference_go.get('mode')}):    "
                        f"{reference_go.get('terms_added', 0):,} specific terms added to "
                        f"{reference_go.get('genes_added', 0):,} genes ({reference_go.get('genes_only_reference', 0):,} "
                        f"had none from hits); {reference_go.get('genes_eligible', 0):,} genes carry a reference symbol\n")
        f.write("\n")

        if modules:
            f.write("Modules:\n")
            f.write("-" * 70 + "\n")
            for m in modules:
                if m.skipped_reason:
                    f.write(f"  {m.name}: not run ({m.skipped_reason})\n")
                    continue
                for line in m.summary_lines():
                    f.write(line + "\n")
            f.write("\n")

        if product_audit:
            f.write("Product naming rules applied:\n")
            f.write("-" * 70 + "\n")
            for rule, n in sorted(product_audit.items(), key=lambda kv: -kv[1]):
                f.write(f"  {n:>8,}  {rule}\n")
            f.write("\n")

        name_concs = [a.name_concordance for a in annotations.values() if a.consensus_name and a.hits]
        sym_concs = [a.symbol_concordance for a in annotations.values() if a.consensus_symbol and a.hits]
        if name_concs or sym_concs:
            f.write("Concordance of the final annotation with each gene's hits:\n")
            f.write("-" * 70 + "\n")
            f.write("  (name: hits whose name matches the product; symbol: hits with a symbol that name the same gene;\n"
                    "   genes without hits are not included)\n")
            if name_concs:
                f.write(f"  Name concordance:   {len(name_concs):,} genes, mean={statistics.mean(name_concs):.2f}, "
                        f"median={statistics.median(name_concs):.2f}, "
                        f"min={min(name_concs):.2f}, max={max(name_concs):.2f}\n")
            if sym_concs:
                f.write(f"  Symbol concordance: {len(sym_concs):,} genes, mean={statistics.mean(sym_concs):.2f}, "
                        f"median={statistics.median(sym_concs):.2f}, "
                        f"min={min(sym_concs):.2f}, max={max(sym_concs):.2f}\n")
            f.write("\n")

        if comparison:
            f.write("Comparison with existing annotation:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Genes with existing symbol:   {comparison['n_existing_symbol']:,}"
                    + (f" ({comparison['n_existing_identifier']:,} are identifiers such as LOC numbers, "
                       f"compared as absent)" if comparison.get('n_existing_identifier') else "") + "\n")
            f.write(f"  Genes with existing product:  {comparison['n_existing_product']:,}\n")
            f.write("  Symbol:  " + ", ".join(f"{k}={v:,}" for k, v in comparison["symbol"].items()) + "\n")
            f.write("  Product: " + ", ".join(f"{k}={v:,}" for k, v in comparison["product"].items()) + "\n")
            f.write(f"  Details: reports/{prefix}_comparison.tsv, reports/{prefix}_comparison_summary.txt\n\n")

    logger.info(f"Results written to {output_dir}")
