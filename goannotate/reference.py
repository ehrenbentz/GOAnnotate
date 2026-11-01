"""
Choice of the reference species whose nomenclature gene symbols follow.

NCBI gene_orthologs links genes to a small set of anchor species. Symbols of
database hits can be converted to a reference species only through the
anchor whose rows cover the hit's species: vertebrates through human,
insects through Drosophila melanogaster, fungi through Saccharomyces
cerevisiae. With --reference-taxid auto (the default) the reference is
chosen from the lineage of the organism being annotated (--taxid):

    Vertebrata (7742)  -> Homo sapiens (9606)
    Arthropoda (6656)  -> Drosophila melanogaster (7227)
    Fungi (4751)       -> Saccharomyces cerevisiae S288C (559292)
    anything else      -> Homo sapiens (9606)

The secondary reference (--secondary-reference-taxid auto) is human when
the primary reference is another species, so that an insect gene carries
its FlyBase symbol and, in its own column, the symbol of its human
ortholog when there is an unambiguous one.
"""

import logging

from goannotate import GOAnnotateError

logger = logging.getLogger(__name__)

# (clade, reference taxon) in order of precedence
AUTO_REFERENCE = (
    (7742, 9606),     # Vertebrata -> human
    (6656, 7227),     # Arthropoda -> Drosophila melanogaster
    (4751, 559292),   # Fungi -> Saccharomyces cerevisiae S288C
)
DEFAULT_REFERENCE = 9606

# Anchor species of NCBI gene_orthologs (rows anchored on these taxa)
ORTHOLOG_ANCHORS = {9606, 7227, 7955, 559292, 330879, 660027, 7460, 7091, 347515, 36329, 237561, 7070,
                    1340429, 185431, 508771}

SPECIES_LABEL = {9606: "Homo sapiens", 7227: "Drosophila melanogaster", 559292: "Saccharomyces cerevisiae S288C",
                 10090: "Mus musculus", 7955: "Danio rerio", 6239: "Caenorhabditis elegans"}


def _parse(value, name):
    if value is None:
        return "auto"
    v = str(value).strip().lower()
    if v in ("auto", "none", "off"):
        return "none" if v == "off" else v
    if not v.isdigit():
        raise GOAnnotateError(f"--{name} must be a taxon ID, 'auto' or 'none' (got {value})")
    return int(v)


def lineage_of(args) -> list:
    """Lineage (taxon IDs, query first) of --taxid from --taxdump or the SQLite database; [] when unknown."""
    from goannotate import taxon
    taxid = getattr(args, "taxid", None)
    if not taxid:
        return []
    lineage = []
    if getattr(args, "taxdump", None):
        lineage = taxon.lineage_from_taxdump(args.taxdump, taxid)
    if not lineage and getattr(args, "db", None):
        lineage = taxon.lineage_sqlite(args.db, taxid)
    return lineage


def resolve_references(args) -> None:
    """
    Replace 'auto' in args.reference_taxid, args.ortholog_anchor_taxid and
    args.secondary_reference_taxid by taxon IDs (secondary may become None).
    """
    ref = _parse(getattr(args, "reference_taxid", "auto"), "reference-taxid")
    anchor = _parse(getattr(args, "ortholog_anchor_taxid", "auto"), "ortholog-anchor-taxid")
    secondary = _parse(getattr(args, "secondary_reference_taxid", "auto"), "secondary-reference-taxid")

    if ref == "auto":
        lineage = set(lineage_of(args))
        ref = DEFAULT_REFERENCE
        clade = None
        for c, r in AUTO_REFERENCE:
            if c in lineage:
                ref, clade = r, c
                break
        if not getattr(args, "taxid", None):
            logger.info(f"Reference species: {SPECIES_LABEL.get(ref, ref)} ({ref}) "
                        f"(default; give --taxid to choose it from the organism's lineage)")
        elif not lineage:
            logger.warning(f"Lineage of taxon {args.taxid} unknown (no taxonomy tables in the database and no "
                           f"--taxdump); reference species {SPECIES_LABEL.get(ref, ref)} ({ref})")
        else:
            how = f"lineage includes taxon {clade}" if clade else "no lineage-specific reference"
            logger.info(f"Reference species: {SPECIES_LABEL.get(ref, ref)} ({ref}), chosen automatically "
                        f"({how})")
    if ref == "none":
        raise GOAnnotateError("--reference-taxid cannot be 'none'")
    if anchor == "auto":
        anchor = ref if ref in ORTHOLOG_ANCHORS else DEFAULT_REFERENCE
    if anchor == "none":
        raise GOAnnotateError("--ortholog-anchor-taxid cannot be 'none'")
    if secondary == "auto":
        secondary = DEFAULT_REFERENCE if ref != DEFAULT_REFERENCE else None
    elif secondary == "none" or secondary == ref:
        secondary = None
    args.reference_taxid = ref
    args.ortholog_anchor_taxid = anchor
    args.secondary_reference_taxid = secondary
    if secondary:
        logger.info(f"Secondary nomenclature: {SPECIES_LABEL.get(secondary, secondary)} ({secondary}), "
                    f"written to the secondary_symbol column")


def label(taxid) -> str:
    return SPECIES_LABEL.get(int(taxid), str(taxid)) if taxid else ""
