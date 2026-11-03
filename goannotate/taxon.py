"""
GO taxon constraints and NCBI taxonomy lineages.

The Gene Ontology states for some terms the taxa they apply to:
  in_taxon (RO:0002162)        the term is only used in this taxon
  never_in_taxon (RO:0002161)  the term is never used in this taxon
Constraints pass to subtypes (is_a) and parts (part_of) of a term. Taxa are
NCBI taxonomy IDs or GO's own union classes (NCBITaxon_Union:0000004,
"Prokaryota"), whose members are listed in go-taxon-groupings.obo.

Sources (downloaded by build_databases.py into the SQLite database, or
given to GOAnnotate as files):
  http://current.geneontology.org/ontology/imports/go-computed-taxon-constraints.obo
  http://current.geneontology.org/ontology/imports/go-taxon-groupings.obo

A GO term is excluded for a query organism when a constraint on the term
or on one of its ancestors (is_a / part_of) is violated by the organism's
lineage: an in_taxon taxon that the lineage does not contain, or a
never_in_taxon taxon that it does.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

IN_TAXON = "in_taxon"
NEVER_IN_TAXON = "never_in_taxon"
_RELATION_IDS = {"RO:0002162": IN_TAXON, "RO:0002160": IN_TAXON, "RO:0002161": NEVER_IN_TAXON}
_TAXON = re.compile(r"(NCBITaxon(?:_Union)?:\d+)")


####################
# OBO PARSING
####################

def parse_taxon_constraints_obo(path) -> dict:
    """GO ID -> list of (relation, taxon) from go-computed-taxon-constraints.obo (or go.obo / go-plus)."""
    constraints = defaultdict(list)
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current = None
            elif line.startswith("id: GO:"):
                current = line.split()[1]
            elif current and (line.startswith("relationship: ") or line.startswith("property_value: ")):
                parts = line.split()
                if len(parts) < 3:
                    continue
                relation = _RELATION_IDS.get(parts[1]) or {"in_taxon": IN_TAXON, "only_in_taxon": IN_TAXON,
                                                           "never_in_taxon": NEVER_IN_TAXON}.get(parts[1])
                if not relation:
                    continue
                m = _TAXON.search(parts[2])
                if m:
                    constraints[current].append((relation, m.group(1)))
    n = sum(len(v) for v in constraints.values())
    logger.info(f"Read {n:,} taxon constraints on {len(constraints):,} GO terms from {Path(path).name}")
    return dict(constraints)


def parse_taxon_groupings_obo(path) -> dict:
    """Union class ID -> set of member taxa (NCBITaxon or NCBITaxon_Union IDs) from go-taxon-groupings.obo."""
    unions = defaultdict(set)
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current = None
            elif line.startswith("id: "):
                current = line[4:].strip()
            elif current and current.startswith("NCBITaxon_Union:") and line.startswith("union_of: "):
                m = _TAXON.search(line[10:])
                if m:
                    unions[current].add(m.group(1))
    logger.info(f"Read {len(unions):,} taxon union classes from {Path(path).name}")
    return dict(unions)


def expand_taxon(taxon: str, unions: dict) -> set:
    """NCBI taxon IDs (ints) that a taxon or union class stands for."""
    out = set()
    stack = [taxon]
    seen = set()
    while stack:
        t = stack.pop()
        if t in seen:
            continue
        seen.add(t)
        if t.startswith("NCBITaxon_Union:"):
            stack.extend(unions.get(t, ()))
        else:
            try:
                out.add(int(t.split(":")[1]))
            except (IndexError, ValueError):
                pass
    return out


####################
# SQLITE
####################

def load_taxon_constraints_sqlite(db_path: str) -> tuple:
    """(constraints, unions) from the go_taxon_constraints and go_taxon_unions tables, or (None, None)."""
    import sqlite3
    conn = sqlite3.connect(db_path)
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "go_taxon_constraints" not in tables:
        conn.close()
        return None, None
    constraints = defaultdict(list)
    for go_id, relation, taxon in conn.execute("SELECT go_id, relation, taxon FROM go_taxon_constraints"):
        constraints[go_id].append((relation, taxon))
    unions = defaultdict(set)
    if "go_taxon_unions" in tables:
        for union_id, member in conn.execute("SELECT union_id, member FROM go_taxon_unions"):
            unions[union_id].add(member)
    conn.close()
    n = sum(len(v) for v in constraints.values())
    logger.info(f"Loaded {n:,} taxon constraints on {len(constraints):,} GO terms and {len(unions):,} "
                f"union classes from the SQLite database")
    return dict(constraints), dict(unions)


def lineage_sqlite(db_path: str, taxid: int) -> list:
    """Taxon IDs from the query taxon up to the root, from the taxonomy_nodes table; [] when unavailable."""
    import sqlite3
    conn = sqlite3.connect(db_path)
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "taxonomy_nodes" not in tables:
        conn.close()
        return []
    if "taxonomy_merged" in tables:
        row = conn.execute("SELECT new_tax_id FROM taxonomy_merged WHERE old_tax_id = ?", (taxid,)).fetchone()
        if row:
            logger.info(f"Taxon {taxid} was merged into {row[0]}")
            taxid = row[0]
    lineage = []
    seen = set()
    t = taxid
    while t and t not in seen:
        seen.add(t)
        row = conn.execute("SELECT parent_tax_id FROM taxonomy_nodes WHERE tax_id = ?", (t,)).fetchone()
        if row is None:
            if t == taxid:
                conn.close()
                return []
            break
        lineage.append(t)
        if row[0] == t:
            break
        t = row[0]
    conn.close()
    return lineage


def lineage_nodes_sqlite(db_path: str, taxid: int) -> list:
    """[(tax_id, rank, name)] from the query taxon up to the root (taxonomy_nodes); [] when unavailable."""
    import sqlite3
    lineage = lineage_sqlite(db_path, taxid)
    if not lineage:
        return []
    conn = sqlite3.connect(db_path)
    columns = {r[1] for r in conn.execute("PRAGMA table_info(taxonomy_nodes)")}
    out = []
    for t in lineage:
        if "rank" in columns:
            row = conn.execute("SELECT rank, name FROM taxonomy_nodes WHERE tax_id = ?", (t,)).fetchone()
        else:
            row = conn.execute("SELECT 'no rank', name FROM taxonomy_nodes WHERE tax_id = ?", (t,)).fetchone()
        out.append((t, (row[0] if row else "no rank") or "no rank", (row[1] if row else "") or ""))
    conn.close()
    return out


# Ranks from the species upward, for the closeness of a relative (orthology module)
RANK_ORDER = ["subspecies", "species", "species group", "subgenus", "genus", "subtribe", "tribe", "subfamily",
              "family", "superfamily", "parvorder", "infraorder", "suborder", "order", "superorder", "infraclass",
              "subclass", "class", "superclass", "subphylum", "phylum", "kingdom", "superkingdom", "domain"]


def lca_rank(query_nodes: list, relative_nodes: list) -> tuple:
    """
    (rank, name) of the lowest common ancestor of two lineages given as
    [(tax_id, rank, name)]; an unranked ancestor takes the rank of the first
    ranked node above it. ('', '') without a common node.
    """
    relative_ids = {t for t, _, _ in relative_nodes}
    for i, (t, rank, name) in enumerate(query_nodes):
        if t in relative_ids:
            for _, r, n in query_nodes[i:]:
                if r in RANK_ORDER:
                    return r, name
            return rank, name
    return "", ""


def rank_within(rank: str, limit: str) -> bool:
    """True when rank is limit or a lower rank (genus is within family); unknown ranks are not."""
    return rank in RANK_ORDER and limit in RANK_ORDER and RANK_ORDER.index(rank) <= RANK_ORDER.index(limit)


def taxon_name_sqlite(db_path: str, taxid: int) -> str:
    import sqlite3
    conn = sqlite3.connect(db_path)
    try:
        row = conn.execute("SELECT name FROM taxonomy_nodes WHERE tax_id = ?", (taxid,)).fetchone()
    except Exception:  # noqa: BLE001
        row = None
    conn.close()
    return row[0] if row else ""


####################
# TAXDUMP FILES
####################

def lineage_from_taxdump(taxdump_dir, taxid: int) -> list:
    """Taxon IDs from the query taxon up to the root, from nodes.dmp (and merged.dmp) in taxdump_dir."""
    taxdump_dir = Path(taxdump_dir)
    nodes_path = taxdump_dir / "nodes.dmp"
    merged_path = taxdump_dir / "merged.dmp"
    if not nodes_path.exists():
        return []
    if merged_path.exists():
        with open(merged_path) as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 2 and parts[0].isdigit() and int(parts[0]) == taxid:
                    taxid = int(parts[1])
                    break
    parent = {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = line.split("|")
            parent[int(parts[0])] = int(parts[1])
    if taxid not in parent:
        return []
    lineage = []
    t = taxid
    while t in parent and t not in lineage:
        lineage.append(t)
        if parent[t] == t:
            break
        t = parent[t]
    return lineage


####################
# APPLICATION
####################

def forbidden_go_terms(constraints: dict, unions: dict, lineage, hierarchy) -> dict:
    """
    GO terms excluded for an organism with the given lineage (set or list of
    taxon IDs). Returns {GO ID: reason}, where reason names the constrained
    ancestor and the relation, for the constrained terms and all their
    descendants through is_a and part_of.
    """
    lineage = set(int(t) for t in lineage)
    forbidden = {}
    for go_id, rules in constraints.items():
        reason = None
        for relation, taxon in rules:
            taxa = expand_taxon(taxon, unions or {})
            if not taxa:
                continue
            if relation == IN_TAXON and not (taxa & lineage):
                reason = f"{go_id} only in {taxon}"
            elif relation == NEVER_IN_TAXON and (taxa & lineage):
                reason = f"{go_id} never in {taxon}"
            if reason:
                break
        if not reason:
            continue
        primary = hierarchy.normalize_id(go_id) or go_id
        forbidden.setdefault(primary, reason)
        for d in hierarchy.get_descendants(primary, include_part_of=True):
            forbidden.setdefault(d, reason)
    return forbidden
