"""Accession-to-GO mappings and NCBI cross-references (flat files and SQLite)."""

import logging
import re
import sys
from collections import defaultdict

logger = logging.getLogger(__name__)


####################
# ACCESSIONS
####################

def extract_accession(subject_id: str) -> list:
    """
    Extract possible accession formats from a BLAST subject ID.
    Returns multiple possibilities to maximize matching chances.
    """
    accessions = [subject_id]

    # UniProt format: sp|P12345|NAME or tr|A0A123|NAME
    if '|' in subject_id:
        parts = subject_id.split('|')
        if len(parts) >= 2:
            accessions.append(parts[1])  # The accession (P12345)
        if len(parts) >= 3:
            accessions.append(parts[2].split('_')[0])  # Gene name part

    # NCBI format: ref|NP_001234.1| or gi|12345|ref|NP_001234.1|
    for match in re.findall(r'[A-Z]{2}_\d+\.?\d*', subject_id):
        accessions.append(match)
        accessions.append(match.split('.')[0])  # Without version

    return accessions


def collect_accessions_from_blast(filtered_hits: dict) -> set:
    """
    Collect all possible accession formats from filtered BLAST results.
    Used to restrict GO mapping loading to relevant entries.

    filtered_hits: dict of transcript_id -> list of (BlastHit, cleaned_name, valid_symbol) tuples
    """
    accessions = set()
    for hits in filtered_hits.values():
        for hit_tuple in hits:
            accessions.update(extract_accession(hit_tuple[0].subject_id))
    return accessions


def get_go_terms_for_hits(hits: list, go_mapping: dict) -> set:
    """Get all GO terms associated with a list of BLAST hits."""
    go_terms = set()
    for hit in hits:
        for acc in extract_accession(hit.subject_id):
            if acc in go_mapping:
                go_terms.update(go_mapping[acc])
    return go_terms


def collect_go_evidence(hits: list, go_mapping: dict) -> tuple:
    """
    GO terms of a list of hits with their evidence: (term -> set of
    evidence codes, term -> number of hits carrying the term). Codes are
    empty for databases that recorded none.
    """
    term_codes = {}
    term_support = {}
    for hit in hits:
        seen = set()
        for acc in extract_accession(hit.subject_id):
            entry = go_mapping.get(acc)
            if not entry:
                continue
            for term, codes in entry.items():
                if term in seen:
                    continue
                seen.add(term)
                s = term_codes.setdefault(term, set())
                if codes:
                    s.update(codes.split(","))
                term_support[term] = term_support.get(term, 0) + 1
    return term_codes, term_support


_GO_WITH_CODES = re.compile(r'^(GO:\d{7})(?::([A-Za-z,]+))?$')


def parse_go_field(text: str) -> list:
    """
    (GO ID, codes) pairs from a GO field. Terms are separated by ';' (or
    ',' in files without evidence codes); a term may carry its evidence
    codes as GO:0005515:IDA,IPI (databases built by the current build_databases.py).
    """
    out = []
    for part in text.split(';'):
        part = part.strip()
        if not part:
            continue
        m = _GO_WITH_CODES.match(part)
        if m:
            out.append((m.group(1), (m.group(2) or "").upper()))
            continue
        for sub in part.split(','):
            sub = sub.strip()
            if sub.startswith('GO:'):
                out.append((sub, ""))
    return out


def _allowed_codes(policy: str):
    from goannotate.go import EVIDENCE_POLICIES
    if policy not in EVIDENCE_POLICIES:
        raise ValueError(f"Unknown GO evidence policy: {policy}")
    return EVIDENCE_POLICIES[policy]


def _add_terms(mapping: dict, accession: str, text: str, allowed, stats: dict) -> None:
    """Add the terms of one GO field to mapping[accession] under the evidence policy."""
    entry = mapping.get(accession)
    for term, codes in parse_go_field(text):
        if codes:
            stats["with_codes"] += 1
            code_set = set(codes.split(','))
            if "ND" in code_set:
                code_set.discard("ND")
                if not code_set:
                    stats["nd"] += 1
                    continue
            if allowed is not None:
                code_set &= allowed
                if not code_set:
                    stats["excluded"] += 1
                    continue
            codes = ",".join(sorted(code_set))
        else:
            stats["no_codes"] += 1
        if entry is None:
            entry = mapping[accession] = {}
        term = sys.intern(term)
        if term in entry and entry[term] and codes:
            merged = set(entry[term].split(',')) | set(codes.split(','))
            entry[term] = sys.intern(",".join(sorted(merged)))
        else:
            entry[term] = sys.intern(codes) if codes else entry.get(term, "")


def _log_evidence(stats: dict, policy: str, mapping: dict) -> None:
    n_terms = sum(len(v) for v in mapping.values())
    logger.info(f"Loaded GO mappings for {len(mapping):,} accessions ({n_terms:,} term assignments)")
    if stats["no_codes"] and not stats["with_codes"]:
        logger.info("  The GO mapping carries no evidence codes (built by an earlier build_databases.py); "
                    "the --go-evidence policy cannot be applied")
    else:
        logger.info(f"  GO evidence policy '{policy}': {stats['excluded']:,} assignments excluded, "
                    f"{stats['nd']:,} ND (no data) assignments dropped")


####################
# FLAT FILES
####################

def load_go_mapping(mapping_file: str, accession_filter: set = None, evidence_policy: str = "all") -> dict:
    """
    Load accession to GO term mapping.

    Expected format (tab-separated):
      accession<TAB>GO:0000001;GO:0000002;...
      accession<TAB>GO:0000001:IDA,IEA;GO:0000002:IEA   (with evidence codes)
    or UniProt ID mapping format:
      accession<TAB>GO:0000001<TAB>GO:0000002<TAB>...

    Returns: dict of accession -> {GO ID: evidence codes (comma-separated, may be "")}
    """
    logger.info(f"Loading GO mapping from {mapping_file}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter)} accessions from BLAST results")
    allowed = _allowed_codes(evidence_policy)

    mapping = {}
    stats = {"excluded": 0, "nd": 0, "no_codes": 0, "with_codes": 0}
    lines_scanned = 0
    lines_loaded = 0

    with open(mapping_file) as f:
        for line in f:
            lines_scanned += 1
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            accession = parts[0].strip()
            if accession_filter is not None and accession not in accession_filter:
                continue

            lines_loaded += 1
            for part in parts[1:]:
                _add_terms(mapping, accession, part, allowed, stats)

    _log_evidence(stats, evidence_policy, mapping)
    if accession_filter:
        logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return mapping


def load_ncbi_idmapping(idmapping_file: str, accession_filter: set = None) -> dict:
    """
    Load UniProt accession to NCBI GeneID mapping from idmapping_selected.tab.
    Column 1: UniProt accession; column 3: NCBI GeneID.

    Returns: dict of UniProt accession -> NCBI GeneID (as string)
    """
    logger.info(f"Loading UniProt->NCBI GeneID mapping from {idmapping_file}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter):,} accessions from BLAST results")

    mapping = {}
    lines_scanned = 0
    lines_loaded = 0

    with open(idmapping_file) as f:
        for line in f:
            lines_scanned += 1
            if lines_scanned % 10000000 == 0:
                logger.info(f"  Scanned {lines_scanned:,} lines...")

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue

            accession = parts[0].strip()
            gene_id = parts[2].strip()
            if not gene_id:
                continue
            if accession_filter is not None and accession not in accession_filter:
                continue

            lines_loaded += 1
            mapping[accession] = gene_id

    logger.info(f"Loaded {len(mapping):,} UniProt->GeneID mappings")
    logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return mapping


def load_ncbi_geneinfo(geneinfo_file: str, geneid_filter: set = None) -> dict:
    """
    Load NCBI GeneID to gene symbol mapping from gene_info.
    Column 2: GeneID; column 3: Symbol.

    Returns: dict of GeneID (as string) -> gene symbol
    """
    logger.info(f"Loading NCBI GeneID->Symbol mapping from {geneinfo_file}")
    if geneid_filter:
        logger.info(f"Filtering to {len(geneid_filter):,} GeneIDs")

    mapping = {}
    lines_scanned = 0
    lines_loaded = 0

    with open(geneinfo_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            lines_scanned += 1
            if lines_scanned % 5000000 == 0:
                logger.info(f"  Scanned {lines_scanned:,} lines...")

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue

            gene_id = parts[1].strip()
            symbol = parts[2].strip()
            if symbol == 'NEWENTRY' or symbol == '-':
                continue
            if geneid_filter is not None and gene_id not in geneid_filter:
                continue

            lines_loaded += 1
            mapping[gene_id] = symbol

    logger.info(f"Loaded {len(mapping):,} GeneID->Symbol mappings")
    logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return mapping


####################
# SQLITE
####################

def _filtered_cursor(conn, sql_all: str, sql_filtered: str, filter_values: set,
                     temp_table: str, temp_col: str):
    """Run a query, restricted to filter_values through a temporary table when given."""
    if filter_values:
        conn.execute(f"CREATE TEMP TABLE {temp_table} ({temp_col} TEXT PRIMARY KEY)")
        batch = [(v,) for v in filter_values]
        for i in range(0, len(batch), 100_000):
            conn.executemany(f"INSERT OR IGNORE INTO {temp_table} VALUES (?)",
                             batch[i:i + 100_000])
        return conn.execute(sql_filtered)
    return conn.execute(sql_all)


def load_go_mapping_sqlite(db_path: str, accession_filter: set = None, evidence_policy: str = "all") -> dict:
    """
    Load accession to GO term mapping from the GOAnnotate SQLite database.
    Returns: dict of accession -> {GO ID: evidence codes (comma-separated, "" when not recorded)}
    """
    import sqlite3

    logger.info(f"Loading GO mapping from SQLite: {db_path}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter):,} accessions from BLAST results")
    allowed = _allowed_codes(evidence_policy)

    conn = sqlite3.connect(db_path)
    mapping = {}
    stats = {"excluded": 0, "nd": 0, "no_codes": 0, "with_codes": 0}

    cursor = _filtered_cursor(
        conn,
        "SELECT accession, go_terms FROM go_mapping",
        "SELECT g.accession, g.go_terms FROM go_mapping g "
        "INNER JOIN _acc_filter f ON g.accession = f.acc",
        accession_filter, "_acc_filter", "acc")

    for accession, go_terms_raw in cursor:
        for part in go_terms_raw.split('\t'):
            _add_terms(mapping, accession, part, allowed, stats)

    conn.close()
    _log_evidence(stats, evidence_policy, mapping)
    return mapping


def load_ncbi_idmapping_sqlite(db_path: str, accession_filter: set = None) -> dict:
    """Load UniProt accession to NCBI GeneID mapping from the SQLite database."""
    import sqlite3

    logger.info(f"Loading UniProt->NCBI GeneID mapping from SQLite: {db_path}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter):,} accessions from BLAST results")

    conn = sqlite3.connect(db_path)
    cursor = _filtered_cursor(
        conn,
        "SELECT uniprot_acc, gene_id FROM ncbi_idmapping",
        "SELECT m.uniprot_acc, m.gene_id FROM ncbi_idmapping m "
        "INNER JOIN _acc_filter f ON m.uniprot_acc = f.acc",
        accession_filter, "_acc_filter", "acc")

    mapping = {}
    for uniprot_acc, gene_id in cursor:
        mapping[uniprot_acc] = gene_id

    conn.close()
    logger.info(f"Loaded {len(mapping):,} UniProt->GeneID mappings")
    return mapping


def load_ncbi_geneinfo_sqlite(db_path: str, geneid_filter: set = None) -> dict:
    """Load NCBI GeneID to gene symbol mapping from the SQLite database."""
    import sqlite3

    logger.info(f"Loading NCBI GeneID->Symbol mapping from SQLite: {db_path}")
    if geneid_filter:
        logger.info(f"Filtering to {len(geneid_filter):,} GeneIDs")

    conn = sqlite3.connect(db_path)
    cursor = _filtered_cursor(
        conn,
        "SELECT gene_id, symbol FROM ncbi_geneinfo",
        "SELECT g.gene_id, g.symbol FROM ncbi_geneinfo g "
        "INNER JOIN _gid_filter f ON g.gene_id = f.gid",
        geneid_filter, "_gid_filter", "gid")

    mapping = {}
    for gene_id, symbol in cursor:
        mapping[gene_id] = symbol

    conn.close()
    logger.info(f"Loaded {len(mapping):,} GeneID->Symbol mappings")
    return mapping


def sqlite_has_table(db_path: str, table: str) -> bool:
    import sqlite3
    conn = sqlite3.connect(db_path)
    row = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                       (table,)).fetchone()
    conn.close()
    return row is not None


def load_ncbi_orthologs_sqlite(db_path: str, geneid_filter: set, anchor_taxid: str,
                               reference_taxid: str) -> dict:
    """
    Map NCBI GeneIDs of any species to the GeneID of their ortholog in the
    reference species, using NCBI gene_orthologs (table ncbi_orthologs).

    gene_orthologs is anchored on a few species (human, tax 9606, for
    vertebrates): every row links an anchor-species gene to its ortholog
    in another species. GeneIDs are resolved to the anchor first; when the
    reference species differs from the anchor, the anchor gene is then
    mapped to its ortholog in the reference species.

    Returns: dict of GeneID -> reference-species GeneID
    """
    import sqlite3

    anchor_taxid = str(int(anchor_taxid))
    reference_taxid = str(int(reference_taxid))
    logger.info(f"Loading NCBI ortholog mapping from SQLite: {db_path} "
                f"(anchor taxon {anchor_taxid}, reference taxon {reference_taxid})")
    conn = sqlite3.connect(db_path)

    to_anchor = {}
    cursor = _filtered_cursor(
        conn,
        f"SELECT other_gene_id, gene_id FROM ncbi_orthologs "
        f"WHERE tax_id = '{anchor_taxid}' AND relationship = 'Ortholog'",
        f"SELECT o.other_gene_id, o.gene_id FROM ncbi_orthologs o "
        f"INNER JOIN _gid_filter f ON o.other_gene_id = f.gid "
        f"WHERE o.tax_id = '{anchor_taxid}' AND o.relationship = 'Ortholog'",
        geneid_filter, "_gid_filter", "gid")
    for other_gid, anchor_gid in cursor:
        to_anchor.setdefault(other_gid, anchor_gid)

    if reference_taxid == anchor_taxid:
        mapping = to_anchor
    else:
        anchor_ids = set(to_anchor.values()) | set(geneid_filter or [])
        anchor_to_ref = {}
        conn.execute("DROP TABLE IF EXISTS _gid_filter")
        cursor = _filtered_cursor(
            conn,
            f"SELECT gene_id, other_gene_id FROM ncbi_orthologs "
            f"WHERE tax_id = '{anchor_taxid}' AND other_tax_id = '{reference_taxid}' "
            f"AND relationship = 'Ortholog'",
            f"SELECT o.gene_id, o.other_gene_id FROM ncbi_orthologs o "
            f"INNER JOIN _gid_filter f ON o.gene_id = f.gid "
            f"WHERE o.tax_id = '{anchor_taxid}' AND o.other_tax_id = '{reference_taxid}' "
            f"AND o.relationship = 'Ortholog'",
            anchor_ids, "_gid_filter", "gid")
        for anchor_gid, ref_gid in cursor:
            anchor_to_ref.setdefault(anchor_gid, ref_gid)
        mapping = {}
        for other_gid, anchor_gid in to_anchor.items():
            ref_gid = anchor_to_ref.get(anchor_gid)
            if ref_gid:
                mapping[other_gid] = ref_gid
        # anchor-species genes themselves map to the reference species
        for anchor_gid, ref_gid in anchor_to_ref.items():
            mapping.setdefault(anchor_gid, ref_gid)

    conn.close()
    logger.info(f"Loaded {len(mapping):,} GeneID -> reference ortholog mappings")
    return mapping


def load_gene_symbol_ids_sqlite(db_path: str, pairs: set) -> dict:
    """
    Map (tax_id, lowercase symbol) pairs to NCBI GeneIDs through the
    ncbi_gene_symbols table.

    Databases built by the current build_databases.py hold three kinds of names per
    gene: the official symbol, the locus tag, and synonyms (column kind).
    A pair is resolved by the most authoritative kind that matches it
    (symbol, then locus_tag, then synonym); a pair that matches several
    GeneIDs within that kind is left out as ambiguous. Older databases
    without the kind column hold official symbols only.
    """
    import sqlite3

    logger.info(f"Looking up GeneIDs for {len(pairs):,} (species, symbol) pairs from SQLite")
    conn = sqlite3.connect(db_path)
    columns = {row[1] for row in conn.execute("PRAGMA table_info(ncbi_gene_symbols)")}
    conn.execute("CREATE TEMP TABLE _pair_filter (tax_id TEXT, symbol TEXT, PRIMARY KEY (tax_id, symbol))")
    batch = list(pairs)
    for i in range(0, len(batch), 100_000):
        conn.executemany("INSERT OR IGNORE INTO _pair_filter VALUES (?, ?)", batch[i:i + 100_000])
    if "kind" in columns:
        cursor = conn.execute(
            "SELECT g.tax_id, g.symbol, g.gene_id, g.kind FROM ncbi_gene_symbols g "
            "INNER JOIN _pair_filter f ON g.tax_id = f.tax_id AND g.symbol = f.symbol")
    else:
        cursor = conn.execute(
            "SELECT g.tax_id, g.symbol, g.gene_id, 'symbol' FROM ncbi_gene_symbols g "
            "INNER JOIN _pair_filter f ON g.tax_id = f.tax_id AND g.symbol = f.symbol")
    rank = {"symbol": 0, "locus_tag": 1, "synonym": 2}
    found = {}   # key -> {rank: set of GeneIDs}
    for tax_id, symbol, gene_id, kind in cursor:
        found.setdefault((tax_id, symbol), {}).setdefault(rank.get(kind, 3), set()).add(gene_id)
    mapping = {}
    ambiguous = 0
    by_kind = {0: 0, 1: 0, 2: 0, 3: 0}
    for key, ranks in found.items():
        best = min(ranks)
        ids = ranks[best]
        if len(ids) == 1:
            mapping[key] = next(iter(ids))
            by_kind[best] += 1
        else:
            ambiguous += 1
    conn.close()
    logger.info(f"  Resolved {len(mapping):,} pairs to GeneIDs ({by_kind[0]:,} by symbol, {by_kind[1]:,} by "
                f"locus tag, {by_kind[2]:,} by synonym; {ambiguous:,} ambiguous pairs skipped)")
    return mapping


class ReferenceNomenclature:
    """
    Official symbols and synonyms of the reference species from NCBI
    gene_info (ncbi_gene_symbols and ncbi_geneinfo). All lookups are
    case-insensitive.

      official:   lowercase symbol -> official symbol as NCBI writes it
      synonyms:   lowercase official symbol -> lowercase synonyms
      renamed:    lowercase synonym -> official symbol, for synonyms that
                  belong to one gene only and are not another gene's
                  official symbol (EIF2C3 -> AGO3, AARS -> AARS1)
    """

    def __init__(self):
        self.official = {}
        self.synonyms = {}
        self.renamed = {}
        self.pseudogenes = set()   # lowercase official symbols of genes NCBI types as pseudo

    def is_pseudogene(self, symbol: str) -> bool:
        return bool(symbol) and symbol.lower() in self.pseudogenes

    def is_readthrough(self, symbol: str) -> bool:
        """
        A readthrough locus of the reference species: an official symbol made
        of two official symbols joined by a hyphen (IQCJ-SCHIP1, KLRC4-KLRK1).
        HLA-A and similar are not, since their parts are not symbols.
        """
        return is_readthrough_symbol(symbol, self.official)

    def is_excluded(self, symbol: str) -> bool:
        """Pseudogenes and readthroughs never name a gene of the query species."""
        return self.is_pseudogene(symbol) or self.is_readthrough(symbol)

    def __bool__(self):
        return bool(self.official)

    def is_official(self, symbol: str) -> bool:
        return bool(symbol) and symbol.lower() in self.official

    def current(self, symbol: str) -> str:
        """The official symbol for a symbol or a unique synonym; '' when unknown."""
        if not symbol:
            return ""
        s = symbol.lower()
        if s in self.official:
            return self.official[s]
        return self.renamed.get(s, "")


# Genome-duplication and copy suffixes of gene symbols (.2, a, b, l, la, lb, l2)
_GENOME_SUFFIX = re.compile(r"(?:\.\d+|L?[AB]|L\d*)$", re.I)
# FlyBase transcript and protein identifiers carried as gene symbols by cDNA-clone entries (kay-RA, CG1647-PB)
_TRANSCRIPT_SUFFIX = re.compile(r"-[RP][A-Z]$")


def gene_symbol_of(gn: str) -> str:
    """A UniProt GN= value as a gene symbol: a FlyBase transcript or protein suffix (-RA, -PB) is dropped."""
    if gn and _TRANSCRIPT_SUFFIX.search(gn) and len(gn) > 3:
        return gn[:-3]
    return gn


def symbol_base(symbol: str) -> str:
    """Upper case without a genome suffix (.2, a, b, l, la, lb, l2): foxo3b -> FOXO3."""
    s = symbol.upper()
    stripped = _GENOME_SUFFIX.sub("", s)
    return stripped or s


_SUFFIX_FORMS = (re.compile(r"\.\d+$"), re.compile(r"[AB]$", re.I), re.compile(r"L[AB]$", re.I),
                 re.compile(r"L\d*$", re.I))


def symbol_bases(symbol: str) -> set:
    """
    Upper-case forms of a symbol without one genome suffix, each suffix read
    separately: glula -> {GLUL, GLU}, foxo3b -> {FOXO3}, krt18a.1 -> {KRT18A}.
    A gene named in one of these forms is the gene the copy belongs to.
    """
    s = symbol.upper()
    out = set()
    for pattern in _SUFFIX_FORMS:
        stripped = pattern.sub("", s)
        if stripped and stripped != s:
            out.add(stripped)
    return out


def is_readthrough_symbol(symbol: str, official) -> bool:
    """symbol is official and made of two or more official symbols joined by hyphens (official: lowercase keys)."""
    if not symbol or "-" not in symbol:
        return False
    parts = symbol.lower().split("-")
    return symbol.lower() in official and len(parts) >= 2 and all(p in official for p in parts)


_SYMBOL_LIKE = re.compile(r"^[a-z0-9][a-z0-9._@/-]*$")
# Disease-locus designations that NCBI gene_info lists among a gene's
# synonyms (CMT2N, DEE29, SCA3, SPG4, DFNA5, LGMD2A, MODY3, PARK2) and
# chromosome-band names of translocation partners (AF6Q21)
_DISEASE_LOCUS = re.compile(
    r"^(?:cmt|dee|eiee|hdls|ttd|sca|spg|mrd|mrx|mrxs|dfna|dfnb|lgmd|als|cmh|cmd|arvd|arvc|lqt|sqt|opa|lca|"
    r"cord|crd|mody|hsan|hmsn|nbia|dyt|epm|eca|atx|mcops|coxpd|mc\ddn|pch|pof|srtd|bpes|hpe|iddm|niddm|eofad|"
    r"park|cdg\d?[a-z]?|ddh|ectd|ehlvic|fhm|hlh|hsp|imd|ipex|jbts|mcph|mgc|mrt|mtdps|nphp|odp|pkd|rp|scar|sccd|"
    r"sgmd|smd|stgd|tnd|wws|xlid)\d+[a-z]?\d*$")
_BAND_NAME = re.compile(r"^af\d+q\d+$")


def _alias_synonym_ok(name: str, official_lower: str, gene_names: set) -> bool:
    """
    Whether an NCBI synonym (lowercase) is usable as an alias of the official
    symbol; gene_names holds the gene's official symbol and all its synonyms.
    """
    if not _SYMBOL_LIKE.match(name):
        return False
    if _DISEASE_LOCUS.match(name) or _BAND_NAME.match(name):
        return False
    # pseudogene named after the gene or one of its synonyms (FKHRL1P2 of FOXO3,
    # FKHRL1 being a FOXO3 synonym)
    m = re.match(r"^(.+)p\d+$", name)
    if m and m.group(1) in gene_names:
        return False
    return True


def load_reference_nomenclature(db_path: str, taxid: str) -> ReferenceNomenclature:
    """Official symbols and synonyms of one species from the SQLite database (empty without the tables)."""
    import sqlite3
    nom = ReferenceNomenclature()
    if not (sqlite_has_table(db_path, "ncbi_gene_symbols") and sqlite_has_table(db_path, "ncbi_geneinfo")):
        return nom
    conn = sqlite3.connect(db_path)
    columns = {row[1] for row in conn.execute("PRAGMA table_info(ncbi_gene_symbols)")}
    if "kind" not in columns:
        conn.close()
        return nom
    has_type = "type" in {row[1] for row in conn.execute("PRAGMA table_info(ncbi_geneinfo)")}
    rows = conn.execute(
        "SELECT s.symbol, s.gene_id, s.kind, g.symbol, " + ("g.type" if has_type else "''") + " "
        "FROM ncbi_gene_symbols s INNER JOIN ncbi_geneinfo g ON g.gene_id = s.gene_id WHERE s.tax_id = ?",
        (str(taxid),)).fetchall()
    conn.close()
    official_of_gene = {}
    synonym_genes = defaultdict(set)
    for name, gene_id, kind, official, gene_type in rows:
        if kind == "symbol":
            nom.official[name] = official
            official_of_gene[gene_id] = official
            if gene_type == "pseudo":
                nom.pseudogenes.add(official.lower())
        elif kind == "synonym":
            synonym_genes[name].add(gene_id)
    names_of_gene = defaultdict(set)
    for name, genes in synonym_genes.items():
        for gid in genes:
            names_of_gene[gid].add(name)
    for gid, off in official_of_gene.items():
        names_of_gene[gid].add(off.lower())
    for name, genes in synonym_genes.items():
        # As aliases only synonyms that read as a gene symbol, belong to one
        # gene and are not another gene's official symbol ("adam 28",
        # "emdc ii", a symbol shared by two genes are left out), and that are
        # not disease-locus or chromosome-band designations or pseudogene
        # names of the gene (see _alias_synonym_ok).
        if len(genes) == 1 and name not in nom.official:
            off = official_of_gene.get(next(iter(genes)))
            if off and _alias_synonym_ok(name, off.lower(), names_of_gene[next(iter(genes))]):
                nom.synonyms.setdefault(off.lower(), set()).add(name)
        if len(genes) == 1 and name not in nom.official:
            off = official_of_gene.get(next(iter(genes)))
            if off:
                nom.renamed[name] = off
    logger.info(f"Reference nomenclature (taxon {taxid}): {len(nom.official):,} official symbols, "
                f"{len(nom.renamed):,} unique synonyms, {len(nom.pseudogenes):,} pseudogenes"
                + ("" if has_type else " (gene types not in the database; rebuild with "
                   "build_databases.py --tables gene_info --update)"))
    return nom


class SymbolResolver:
    """
    Resolves gene symbols for Diamond hits.

    With an ortholog map, the symbol of a hit from another species is
    replaced by the symbol of its ortholog in the reference species, so
    that every gene is named in one nomenclature. Hits from the reference
    species keep their own symbol. Without an ortholog map, the GN= symbol
    is used as is and the NCBI cross-reference serves only as a fallback
    for hits without a GN= symbol.
    """

    def __init__(self, uniprot_to_geneid: dict = None, geneid_to_symbol: dict = None,
                 ortholog_map: dict = None, reference_taxid: str = "9606",
                 symbol_to_geneid: dict = None):
        self.uniprot_to_geneid = uniprot_to_geneid or {}
        self.geneid_to_symbol = geneid_to_symbol or {}
        self.ortholog_map = ortholog_map
        self.reference_taxid = reference_taxid
        self.symbol_to_geneid = symbol_to_geneid or {}   # (tax_id, lowercase symbol) -> GeneID
        self.reference_symbols = frozenset()             # lowercase official reference-species symbols
        self.reference_pseudogenes = frozenset()         # lowercase reference symbols of pseudogenes

    def is_reference_symbol(self, symbol: str) -> bool:
        """True when a symbol, in any case, is an official symbol of the reference species (il15 -> IL15)."""
        return bool(symbol) and symbol.lower() in self.reference_symbols

    def is_excluded(self, symbol: str) -> bool:
        """A reference pseudogene (CMAHP) or readthrough (IQCJ-SCHIP1) symbol: never names a query gene."""
        return bool(symbol) and (symbol.lower() in self.reference_pseudogenes
                                 or is_readthrough_symbol(symbol, self.reference_symbols))

    @property
    def has_ncbi(self) -> bool:
        return bool(self.uniprot_to_geneid) and bool(self.geneid_to_symbol)

    @property
    def has_orthologs(self) -> bool:
        return bool(self.ortholog_map)

    def gene_id_for(self, accession: str, taxid: str = "", gn_symbol: str = "") -> str:
        """NCBI GeneID of a hit: from the accession, else from its (species, symbol)."""
        gene_id = self.uniprot_to_geneid.get(accession)
        if gene_id:
            return gene_id
        if taxid and gn_symbol:
            return self.symbol_to_geneid.get((taxid, gn_symbol.lower()), "")
        return ""

    def reference_symbol(self, accession: str, taxid: str = "", gn_symbol: str = "") -> str:
        """Symbol of the reference-species ortholog of a hit, or empty."""
        if not self.ortholog_map:
            return ""
        gene_id = self.gene_id_for(accession, taxid, gn_symbol)
        if not gene_id:
            return ""
        ref_gid = self.ortholog_map.get(gene_id)
        if not ref_gid:
            return ""
        return self.geneid_to_symbol.get(ref_gid, "")

    def ncbi_symbol(self, accession: str) -> str:
        """Symbol from the NCBI cross-reference (reference ortholog first, then the gene itself)."""
        symbol = self.reference_symbol(accession)
        if symbol:
            return symbol
        gene_id = self.uniprot_to_geneid.get(accession)
        if gene_id:
            return self.geneid_to_symbol.get(gene_id, "")
        return ""

    def resolve(self, hit, gn_symbol: str, is_bad_symbol) -> tuple:
        """
        Symbol and origin for a hit: ('uniprot' for the hit's own GN= symbol,
        'ortholog' for a symbol taken from the reference-species ortholog,
        '' for none).
        """
        if self.ortholog_map and hit.taxid != self.reference_taxid:
            symbol = self.reference_symbol(hit.accession, hit.taxid, gn_symbol)
            # a reference pseudogene (CMAHP) does not name the functional gene
            # of another species; the hit keeps its own symbol (Cmah)
            if symbol and not is_bad_symbol(symbol) and symbol.lower() not in self.reference_pseudogenes:
                return symbol, "ortholog"
        # The hit's own symbol is used only when it is not a genome-specific
        # identifier (LOC number, locus tag, clone ID); such symbols were
        # still useful above for reaching the gene's reference ortholog.
        if gn_symbol and not is_bad_symbol(gn_symbol):
            return gn_symbol, "uniprot"
        return "", ""


####################
# CONVENIENCE
####################

def load_reference_data(db: str, go_mapping_file: str, ncbi_idmapping: str,
                        ncbi_geneinfo: str, accessions: set,
                        resolve_orthologs: bool = True, reference_taxid: str = "9606",
                        anchor_taxid: str = "9606", symbol_pairs: set = None,
                        go_evidence_policy: str = "all") -> tuple:
    """
    Load the GO mapping and NCBI cross-references, restricted to the
    accessions present in the filtered BLAST hits.

    symbol_pairs: (tax_id, lowercase GN symbol) pairs of hits from species
    other than the reference, used to reach GeneIDs for accessions that
    the UniProt ID mapping does not cover.
    go_evidence_policy: which GO evidence codes a term needs (all,
    non-electronic, experimental).

    Returns (go_mapping, SymbolResolver).
    """
    uniprot_to_geneid = None
    geneid_to_symbol = None
    ortholog_map = None
    symbol_to_geneid = None

    if db:
        go_mapping = load_go_mapping_sqlite(db, accession_filter=accessions, evidence_policy=go_evidence_policy)
        uniprot_to_geneid = load_ncbi_idmapping_sqlite(db, accession_filter=accessions)
        geneids = set(uniprot_to_geneid.values())
        if resolve_orthologs:
            if sqlite_has_table(db, "ncbi_orthologs"):
                if symbol_pairs and sqlite_has_table(db, "ncbi_gene_symbols"):
                    symbol_to_geneid = load_gene_symbol_ids_sqlite(db, symbol_pairs)
                    geneids = geneids | set(symbol_to_geneid.values())
                elif symbol_pairs:
                    logger.warning("The SQLite database has no ncbi_gene_symbols table; only accessions "
                                   "with a GeneID in the UniProt ID mapping can be resolved. "
                                   "Run build_databases.py --add-orthologs.")
                ortholog_map = load_ncbi_orthologs_sqlite(db, geneids, anchor_taxid, reference_taxid)
                geneids = geneids | set(ortholog_map.values())
            else:
                logger.warning("The SQLite database has no ncbi_orthologs table; symbols will not be "
                               "resolved to the reference species. Run build_databases.py --add-orthologs.")
        geneid_to_symbol = load_ncbi_geneinfo_sqlite(db, geneid_filter=geneids)
    else:
        go_mapping = load_go_mapping(go_mapping_file, accession_filter=accessions,
                                     evidence_policy=go_evidence_policy)
        if ncbi_idmapping and ncbi_geneinfo:
            uniprot_to_geneid = load_ncbi_idmapping(ncbi_idmapping, accession_filter=accessions)
            geneid_to_symbol = load_ncbi_geneinfo(
                ncbi_geneinfo, geneid_filter=set(uniprot_to_geneid.values()))
        if resolve_orthologs:
            logger.info("Ortholog resolution requires the SQLite database (--db); skipped")

    resolver = SymbolResolver(uniprot_to_geneid, geneid_to_symbol, ortholog_map, reference_taxid,
                              symbol_to_geneid)
    nomenclature = load_reference_nomenclature(db, reference_taxid) if db else ReferenceNomenclature()
    resolver.reference_symbols = frozenset(nomenclature.official)
    resolver.reference_pseudogenes = frozenset(nomenclature.pseudogenes)
    return go_mapping, resolver, nomenclature
