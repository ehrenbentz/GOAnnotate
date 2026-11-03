"""
GO terms from the reference gene of a gene's final symbol.

Once a gene carries a supported symbol of the reference species (or a
close relative's own form of one, through the replaced reference symbol; a
reciprocal best hit, a supported vote, an ortholog symbol from a close
relative, a fill), that symbol names one reference gene whose UniProt
entries carry curated GO annotations. Those terms are transferred to the
query gene (--reference-go add, the default): terms the hits already gave
are kept, terms new to the gene are added, and ancestors of terms already
present add nothing. With fill they go only to genes without GO terms from
hits; off leaves GO to the hits and InterProScan.

The reference gene's entries are the proteins of the reference proteome
(SwissProt and TrEMBL of the reference species, as built for the reciprocal
searches) that carry the symbol in GN=; their GO terms come from the
database's GO mapping under the run's --go-evidence policy, and the taxon
constraints and --namespace filters of the run apply.
"""

import logging
from collections import defaultdict

from goannotate.go import GO_ROOT_TERMS
from goannotate import mappings

logger = logging.getLogger(__name__)

# Symbols the reference gene's GO terms are transferred through: supported by the reference
# proteome, or taken from the relatives' orthologs
REFERENCE_GO_STATUSES = {"rbh", "rbh_override", "rbh_added", "vote_supported", "reference_best_hit",
                         "foreign_ortholog", "ortholog", "ortholog_fill"}


def add_reference_go(annotations: dict, ref_symbols: dict, db: str, go_hierarchy, forbidden_go: set,
                     namespaces: list, evidence_policy: str = "all", mode: str = "add") -> dict:
    """
    Add the reference gene's GO terms to each eligible gene in place and
    record them in annot.go_reference (GO ID -> accessions carrying it).

    ref_symbols: reference proteome accession -> GN= symbol (assign.read_reference_symbols).
    Returns the counts.
    """
    stats = {"mode": mode, "genes_eligible": 0, "genes_with_reference_go": 0, "genes_added": 0,
             "terms_added": 0, "genes_only_reference": 0}
    if mode == "off":
        return stats
    if not db:
        logger.info("Reference-gene GO terms need the SQLite database (--db); skipped")
        return stats
    by_symbol = defaultdict(set)
    for acc, sym in ref_symbols.items():
        by_symbol[sym.upper()].add(acc)
    eligible = {}
    for gene, annot in annotations.items():
        if not annot.consensus_symbol or (annot.symbol_status or "") not in REFERENCE_GO_STATUSES:
            continue
        accs = by_symbol.get(annot.consensus_symbol.upper())
        if not accs and (annot.symbol_status or "") in ("ortholog", "ortholog_fill"):
            # a species-specific form (kdm2aa) of a reference gene: the close relatives' symbol in the
            # reference nomenclature names it; for a refined symbol (the same gene as the core) also the
            # replaced core symbol and the aliases. A taken symbol names another gene than the replaced
            # one, whose terms are never used.
            alts = [annot.ortholog_reference_symbol]
            if annot.ortholog_evidence != "close_taken":
                alts += [annot.voted_symbol] + list(annot.symbol_aliases)
            for alt in alts:
                if alt and by_symbol.get(alt.upper()):
                    accs = by_symbol[alt.upper()]
                    break
        if accs:
            eligible[gene] = accs
    stats["genes_eligible"] = len(eligible)
    if not eligible:
        logger.info("Reference-gene GO terms: no gene carries a symbol of the reference proteome")
        return stats
    accessions = set().union(*eligible.values())
    go_map = mappings.load_go_mapping_sqlite(db, accession_filter=accessions, evidence_policy=evidence_policy)

    for gene, accs in eligible.items():
        annot = annotations[gene]
        if mode == "fill" and annot.go_terms:
            continue
        terms = {}          # GO ID -> (accessions, codes)
        for acc in accs:
            for t, codes in go_map.get(acc, {}).items():
                norm = go_hierarchy.normalize_id(t) if go_hierarchy is not None else t
                if not norm or norm in GO_ROOT_TERMS or norm in forbidden_go:
                    continue
                entry = terms.setdefault(norm, (set(), set()))
                entry[0].add(acc)
                if codes:
                    entry[1].update(codes.split(","))
        if namespaces and terms and go_hierarchy is not None:
            keep = go_hierarchy.filter_by_namespace(set(terms), namespaces)
            terms = {t: v for t, v in terms.items() if t in keep}
        if not terms:
            continue
        stats["genes_with_reference_go"] += 1
        had_go = bool(annot.go_terms)
        before = set(annot.go_terms)
        annot.go_terms = set(annot.go_terms) | set(terms)
        annot.specific_go_terms = (go_hierarchy.filter_to_specific(annot.go_terms) if go_hierarchy is not None
                                   else set(annot.go_terms))
        for t, (t_accs, codes) in terms.items():
            annot.go_reference[t] = sorted(t_accs)
            if t not in before:
                annot.go_evidence[t] = ",".join(sorted(codes)) if codes else ""
        added = [t for t in annot.specific_go_terms if t in terms and t not in before]
        if added:
            stats["genes_added"] += 1
            stats["terms_added"] += len(added)
            if not had_go:
                stats["genes_only_reference"] += 1
    logger.info(f"Reference-gene GO terms ({mode}): {stats['genes_eligible']:,} genes carry a symbol of the reference "
                f"proteome; {stats['genes_added']:,} genes received {stats['terms_added']:,} specific terms they did not "
                f"have ({stats['genes_only_reference']:,} had no GO terms from hits)")
    return stats
