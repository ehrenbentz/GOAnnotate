#!/usr/bin/env python3
"""
GOAnnotate: Gene Ontology Annotation Pipeline using Blastx and UniProt
Written by E. J. Bentz (2025)

Required input files:
1) A FASTA file containing all the sequences to be annotated
2) Diamond/Blastx results obtained by blasting sequences to appropriate SwissProt and/or TrEMBL protein databases
3) A GO mapping file (created from the uniprot goa.gaf file)
4) The current go.obo hierarchy definitions file
5) A bad_names.txt file containing patterns that will NOT be used to annotate genes

Usage:
  ./GOAnnotate.py --transcripts transcripts.fasta --blast-results diamond_results.tsv \\
                  --go-mapping GO_mapping.tsv --bad-names bad_names.txt --go-obo go.obo -o output_dir

Note: The --transcripts file defines the complete gene set. Gene IDs are extracted from the first
field of each FASTA header (the portion before the first whitespace). All genes from this file
will appear in the output, including those without BLAST hits. Statistics are calculated based
on the number of unique gene IDs, not just those with BLAST hits.

"""
import sys
import re
import argparse
import logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class BlastStats:
    """Statistics from BLAST results parsing."""
    total_hits_in_file: int = 0
    hits_filtered_by_evalue: int = 0
    hits_after_evalue_filter: int = 0
    queries_in_file: int = 0
    queries_after_evalue_filter: int = 0
    swissprot_hits: int = 0
    trembl_hits: int = 0
    other_hits: int = 0


# =============================================================================
# BAD NAMES FILTER
# =============================================================================

class BadNameFilter:
    """
    Filter for removing uninformative or spurious gene/protein names.

    Supports three types of patterns:
      1. Exact matches (case-insensitive) - patterns <=6 chars
      2. Substring matches (case-insensitive) - patterns >6 chars
      3. Regex patterns (prefixed with 'regex:')
    """

    def __init__(self, patterns_file: str):
        self.exact_matches = set()
        self.substrings = []
        self.regex_patterns = []

        if not Path(patterns_file).exists():
            raise FileNotFoundError(f"Bad names pattern file not found: {patterns_file}")

        with open(patterns_file) as f:
            self._parse_patterns(f.read())

        logger.info(
            f"BadNameFilter: loaded {len(self.exact_matches)} exact, "
            f"{len(self.substrings)} substring, {len(self.regex_patterns)} regex patterns "
            f"from '{patterns_file}'"
        )

    def _parse_patterns(self, text: str):
        """Parse pattern text into exact, substring, and regex categories."""
        for line in text.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('regex:'):
                pattern = line[6:].strip()
                try:
                    self.regex_patterns.append(re.compile(pattern, re.IGNORECASE))
                except re.error as e:
                    logger.warning(f"Invalid regex pattern '{pattern}': {e}")
            else:
                # Short patterns (<=6 chars) are exact matches, longer ones are substrings
                lower = line.lower()
                if len(line) <= 6:
                    self.exact_matches.add(lower)
                else:
                    self.substrings.append(lower)

    def is_bad_name(self, name: str) -> bool:
        """Check if a name matches any bad pattern."""
        if not name:
            return True

        name_lower = name.lower().strip()

        # Check exact matches
        if name_lower in self.exact_matches:
            return True

        # Check substrings
        for substring in self.substrings:
            if substring in name_lower:
                return True

        # Check regex patterns
        for pattern in self.regex_patterns:
            if pattern.search(name):
                return True

        return False

    def filter_names(self, names: list) -> list:
        """Return only names that are NOT bad."""
        return [n for n in names if not self.is_bad_name(n)]


# =============================================================================
# BLAST HIT DATA STRUCTURE
# =============================================================================

@dataclass
class BlastHit:
    """Represents a single BLAST hit."""
    query_id: str
    subject_id: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    stitle: str = ""

    @property
    def gene_name(self) -> str:
        """Extract gene/protein name from subject title."""
        if not self.stitle:
            return self.subject_id

        # Common patterns in BLAST titles:
        # >sp|P12345|GENE_SPECIES Description OS=Species ...
        # >tr|A0A123|A0A123_SPECIES Description OS=Species ...
        # >gi|123456|ref|NP_001234.1| Description [Species]
        # >ENSP00000123456 Description

        title = self.stitle

        # Remove OS= and everything after (UniProt format)
        if ' OS=' in title:
            title = title.split(' OS=')[0]

        # Remove [Species] suffix
        title = re.sub(r'\s*\[.*?\]\s*$', '', title)

        # Remove GN= gene name tag if present (we want the description)
        title = re.sub(r'\s*GN=\S+', '', title)

        # Remove PE= and SV= tags
        title = re.sub(r'\s*PE=\d+', '', title)
        title = re.sub(r'\s*SV=\d+', '', title)

        # Remove ECO evidence codes from TrEMBL entries {ECO:...}
        title = re.sub(r'\s*\{ECO:[^}]+\}', '', title)

        return title.strip()

    @property
    def gene_symbol(self) -> str:
        """Extract gene symbol from GN= field in UniProt title."""
        if not self.stitle:
            return ""

        # Look for GN=SYMBOL pattern
        match = re.search(r'GN=(\S+)', self.stitle)
        if match:
            return match.group(1)

        return ""

    @property
    def is_swissprot(self) -> bool:
        """Check if this hit is from Swiss-Prot (reviewed) vs TrEMBL (unreviewed)."""
        return self.subject_id.startswith('sp|')

    @property
    def is_trembl(self) -> bool:
        """Check if this hit is from TrEMBL (unreviewed)."""
        return self.subject_id.startswith('tr|')

    def is_valid_gene_symbol(self, symbol: str) -> bool:
        """Check if a gene symbol is informative (not a placeholder or internal ID)."""
        if not symbol:
            return False

        symbol_upper = symbol.upper()

        # Filter out LOC numbers (NCBI automatic naming)
        if re.match(r'^LOC\d+$', symbol_upper):
            return False

        # Filter out zebrafish clone IDs (SI:*, si:dkey-*, zgc:*, wu:*, zmp:*, etc.)
        if re.match(r'^SI:', symbol, re.IGNORECASE):
            return False
        if re.match(r'^zgc:\d+$', symbol, re.IGNORECASE):
            return False
        if re.match(r'^wu:', symbol, re.IGNORECASE):
            return False
        if re.match(r'^im:', symbol, re.IGNORECASE):
            return False
        if re.match(r'^zmp:', symbol, re.IGNORECASE):
            return False

        # Filter out internal lab/genome annotation naming patterns
        # Pattern: alphanumeric prefix + underscore + 5+ digits
        # Examples: D9C73_027734, KOW79_015307, JOB18_020029, F2P81_025073
        if re.match(r'^[A-Z0-9]+_\d{5,}$', symbol_upper):
            return False

        # Pattern: PREFIX_LETTER+DIGITS (4+ digits)
        # Examples: XNOV1_A014589, JOQ06_022217, PECUL_23A059567
        if re.match(r'^[A-Z0-9]+_[A-Z]*\d{4,}$', symbol_upper):
            return False

        # Genome annotation style: letter + digits + underscore + digits + G + digits
        # Examples: D4764_01G0004220, D4764_11G0002770
        if re.match(r'^[A-Z]\d+_\d+G\d+$', symbol_upper):
            return False

        # Fish genome annotation patterns
        # D5F01_LYC10690 style: PREFIX_PREFIXNUMBER
        if re.match(r'^[A-Z]\d+[A-Z]+_[A-Z]+\d+$', symbol_upper):
            return False

        # Nfu_g_1_016631 style (Nile tilapia/Nothobranchius furzeri genome)
        if re.match(r'^[A-Z][A-Z]+_[A-Z]_\d+_\d+$', symbol_upper):
            return False

        # EXN66_Car014809 style: PREFIX_PrefixNumber
        if re.match(r'^[A-Z0-9]+_[A-Z][A-Z]+\d{4,}$', symbol_upper):
            return False

        # DR999_PMT18446, E1301_Tti021464 style
        if re.match(r'^[A-Z]\d+_[A-Z]+\d+$', symbol_upper):
            return False

        # LOCUS patterns: MMEN_LOCUS7967, PLEPLA_LOCUS46505
        if re.match(r'^[A-Z]+_LOCUS\d+$', symbol_upper):
            return False

        # Medaka/Oryzias latipes: OLA.12830, OLA.95
        if re.match(r'^OLA\.\d+$', symbol_upper):
            return False

        # Xenopus: XELAEV_18009728mg
        if re.match(r'^XELAEV_\d+', symbol_upper):
            return False

        # Long alphanumeric genome IDs: GSTENG00028376001, GSONMT00077921001
        if re.match(r'^[A-Z]{3,}\d{8,}$', symbol_upper):
            return False

        # PPUP9740 style (short prefix + many digits)
        if re.match(r'^[A-Z]{2,4}\d{4,}$', symbol_upper):
            return False

        # Accession numbers with version (e.g., CU459095.1, AB123456.2, CABZ01015475.1)
        if re.match(r'^[A-Z]{1,4}\d{5,}\.\d+$', symbol_upper):
            return False

        # Filter out automatic paralog/copy numbering patterns with underscores
        # Examples: LIN1_28, Pol_27, PO21_2, CFDP2_13, G2E3_0, Ppm1l_0, YME1L1_1
        # Pattern: alphanumeric name (2-10 chars) + underscore + small number (1-2 digits)
        if re.match(r'^[A-Z0-9]{2,10}_\d{1,2}$', symbol_upper):
            return False

        # PECUL_23A059567 style: PREFIX_MIXED ending in many digits
        if re.match(r'^[A-Z]+_[A-Z0-9]*\d{5,}$', symbol_upper):
            return False

        # Filter out transposon/retrotransposon naming
#        if re.match(r'^POL_?\d*$', symbol_upper):
#            return False
#        if re.match(r'^PO\d+$', symbol_upper):
#            return False
#        if re.match(r'^TPASE$', symbol_upper):
#            return False
#        if re.match(r'^TY3B-', symbol_upper):
#            return False

        # Filter out KIAA protein project names (uninformative)
        if re.match(r'^KIAA\d+$', symbol_upper):
            return False

        # Filter out Ensembl-style IDs
        if re.match(r'^ENS[A-Z]+\d+$', symbol_upper):
            return False

        # Filter out patterns like GRMZM, Glyma, etc. (plant genome IDs)
        if re.match(r'^GRMZM\d+G\d+$', symbol_upper):
            return False
        if re.match(r'^Glyma\.\d+G\d+$', symbol, re.IGNORECASE):
            return False

        # Filter out chromosome ORF naming (C1orf43, C2orf74, CXorf38, etc.)
        if re.match(r'^C[0-9X]+ORF\d+$', symbol_upper):
            return False
        # Filter out cross-species chromosome ORF naming (e.g., C1H1ORF43 - chicken homolog)
        if re.match(r'^C\d+H\d+ORF\d+$', symbol_upper):
            return False
        # Filter out CUNH#ORF# pattern (zebrafish/other fish referencing human ORFs)
        if re.match(r'^CUNH\d+ORF\d+$', symbol_upper):
            return False
        # Filter out LG#H#ORF# pattern (linkage group + human ORF reference)
        if re.match(r'^LG\d+H\d+ORF\d+$', symbol_upper):
            return False
        # Filter out gene names with C#ORF# suffix (e.g., ZHX1-C8ORF76)
        if re.search(r'-C\d+ORF\d+$', symbol_upper):
            return False

        # Filter out generic ORF labels (just "ORF" with optional number)
        if re.match(r'^ORF\d*$', symbol_upper):
            return False

        # Filter out RIKEN clone IDs (mouse cDNA project naming)
        # Pattern: digits + letters + digits + "Rik" (e.g., 4933434E20Rik, 1700010I14Rik)
        if re.match(r'^\d+[A-Z]+\d+RIK$', symbol_upper):
            return False

        # Filter out very short symbols that are likely not real (1-2 chars)
        if len(symbol) < 2:
            return False

        # Filter out symbols that are mostly numbers (e.g., "123456")
        if re.match(r'^\d+$', symbol):
            return False

        return True


@dataclass
class QueryAnnotation:
    """Stores annotation results for a single query sequence."""
    query_id: str
    hits: list = field(default_factory=list)
    filtered_hits: list = field(default_factory=list)
    consensus_name: str = ""
    consensus_symbol: str = ""
    go_terms: set = field(default_factory=set)
    specific_go_terms: set = field(default_factory=set)


# =============================================================================
# GO HIERARCHY (for filtering to specific terms)
# =============================================================================

class GOHierarchy:
    """
    Parse GO OBO file and track parent-child relationships.
    Used to filter GO terms to the most specific (leaf-like) terms.
    """

    def __init__(self, obo_file: str):
        self.names = {}              # GO ID -> name
        self.namespace = {}          # GO ID -> namespace (BP/MF/CC)
        self.parents = defaultdict(set)   # GO ID -> set of parent IDs
        self.children = defaultdict(set)  # GO ID -> set of child IDs
        self.alt_ids = {}            # alt_id -> primary GO ID

        self._parse_obo(obo_file)
        logger.info(
            f"Loaded GO hierarchy: {len(self.names)} terms "
            f"({sum(1 for ns in self.namespace.values() if ns == 'biological_process')} BP, "
            f"{sum(1 for ns in self.namespace.values() if ns == 'molecular_function')} MF, "
            f"{sum(1 for ns in self.namespace.values() if ns == 'cellular_component')} CC)"
        )

    def _parse_obo(self, obo_file: str):
        """Parse the OBO file."""
        current_id = None
        current_name = None
        current_ns = None
        current_alt_ids = []
        is_obsolete = False

        with open(obo_file) as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    # Save previous term
                    if current_id and current_name and current_ns and not is_obsolete:
                        self.names[current_id] = current_name
                        self.namespace[current_id] = current_ns
                        for alt_id in current_alt_ids:
                            self.alt_ids[alt_id] = current_id

                    # Reset for new term
                    current_id = None
                    current_name = None
                    current_ns = None
                    current_alt_ids = []
                    is_obsolete = False
                    continue

                if line.startswith("[") and line != "[Term]":
                    # Non-term stanza, save previous term if valid
                    if current_id and current_name and current_ns and not is_obsolete:
                        self.names[current_id] = current_name
                        self.namespace[current_id] = current_ns
                        for alt_id in current_alt_ids:
                            self.alt_ids[alt_id] = current_id
                    current_id = None
                    continue

                if line.startswith("id: GO:"):
                    current_id = line.split()[1]
                elif line.startswith("alt_id: GO:"):
                    current_alt_ids.append(line.split()[1])
                elif line.startswith("name:"):
                    current_name = line[5:].strip()
                elif line.startswith("namespace:"):
                    current_ns = line.split()[1]
                elif line.startswith("is_a: GO:"):
                    parent_id = line.split()[1]
                    if current_id:
                        self.parents[current_id].add(parent_id)
                        self.children[parent_id].add(current_id)
                elif line == "is_obsolete: true":
                    is_obsolete = True

            # Don't forget the last term
            if current_id and current_name and current_ns and not is_obsolete:
                self.names[current_id] = current_name
                self.namespace[current_id] = current_ns
                for alt_id in current_alt_ids:
                    self.alt_ids[alt_id] = current_id

    def normalize_id(self, go_id: str) -> Optional[str]:
        """Normalize a GO ID, resolving alt_ids to primary IDs."""
        if go_id in self.names:
            return go_id
        if go_id in self.alt_ids:
            return self.alt_ids[go_id]
        return None

    def get_ancestors(self, go_id: str) -> set:
        """Get all ancestors of a GO term."""
        ancestors = set()
        stack = list(self.parents.get(go_id, []))
        while stack:
            parent = stack.pop()
            if parent not in ancestors:
                ancestors.add(parent)
                stack.extend(self.parents.get(parent, []))
        return ancestors

    def get_descendants(self, go_id: str) -> set:
        """Get all descendants of a GO term."""
        descendants = set()
        stack = list(self.children.get(go_id, []))
        while stack:
            child = stack.pop()
            if child not in descendants:
                descendants.add(child)
                stack.extend(self.children.get(child, []))
        return descendants

    def filter_to_specific(self, go_terms: set) -> set:
        """
        Filter a set of GO terms to retain only the most specific terms.
        Removes any term that is an ancestor of another term in the set.
        """
        if not go_terms:
            return set()

        # Normalize all IDs first
        normalized = set()
        for term in go_terms:
            norm = self.normalize_id(term)
            if norm:
                normalized.add(norm)

        if not normalized:
            return set()

        # For each term, check if any other term in the set is its descendant
        # If so, this term is not specific (it's an ancestor of something more specific)
        specific = set()
        for term in normalized:
            descendants = self.get_descendants(term)
            # If none of our terms are descendants of this term, it's specific
            if not descendants.intersection(normalized):
                specific.add(term)

        return specific

    def filter_by_namespace(self, go_terms: set, namespaces: list) -> set:
        """Filter GO terms to only include those in specified namespaces."""
        ns_map = {
            'BP': 'biological_process',
            'MF': 'molecular_function',
            'CC': 'cellular_component',
            'biological_process': 'biological_process',
            'molecular_function': 'molecular_function',
            'cellular_component': 'cellular_component'
        }
        target_ns = {ns_map.get(ns, ns) for ns in namespaces}
        return {t for t in go_terms if self.namespace.get(t) in target_ns}


# =============================================================================
# FASTA HEADER PARSING
# =============================================================================

def parse_fasta_headers(fasta_file: str) -> tuple:
    """
    Parse a FASTA file and extract gene IDs from headers.

    Extracts the first whitespace-delimited field from each header line
    (the portion after '>' up to the first space/tab).

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        Tuple of (total_sequences, list of unique gene IDs in order)
    """
    logger.info(f"Parsing gene IDs from {fasta_file}")

    gene_ids = []
    seen = set()
    total_sequences = 0

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                total_sequences += 1
                # Extract first field (gene ID) from header
                header = line[1:].strip()
                gene_id = header.split()[0] if header else ""

                if gene_id and gene_id not in seen:
                    gene_ids.append(gene_id)
                    seen.add(gene_id)

    logger.info(f"Found {total_sequences:,} sequences, {len(gene_ids):,} unique gene IDs")
    return total_sequences, gene_ids


# =============================================================================
# BLAST RESULTS PARSING
# =============================================================================

def parse_blast_results(blast_file: str, top_n: int = 25,
                        evalue_threshold: float = 1e-5,
                        prefer_swissprot: bool = True) -> tuple:
    """
    Parse BLAST/Diamond output (format 6 with stitle) and return top N hits per query.

    When prefer_swissprot is True, SwissProt (sp|) hits are prioritized over TrEMBL (tr|) hits.
    This ensures that well-annotated SwissProt entries are used when available, even if
    TrEMBL hits have slightly higher bitscores. Hits are sorted by: SwissProt first, then
    by bitscore within each category.

    Args:
        blast_file: Path to BLAST output file
        top_n: Maximum hits to keep per query
        evalue_threshold: E-value cutoff for filtering hits (default: 1e-5)
        prefer_swissprot: Prioritize SwissProt hits over TrEMBL (default: True)

    Returns: tuple of (dict of query_id -> list of BlastHit objects, BlastStats)
    """
    logger.info(f"Parsing BLAST results from {blast_file}")
    logger.info(f"E-value threshold: {evalue_threshold}")

    if prefer_swissprot:
        logger.info("SwissProt hits will be prioritized over TrEMBL hits")

    # Collect ALL hits first (no top_n limit during parsing)
    # This allows proper sorting before truncation
    results = defaultdict(list)
    all_queries_in_file = set()  # Track all unique queries before filtering
    stats = BlastStats()

    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            stats.total_hits_in_file += 1
            query_id = parts[0]
            all_queries_in_file.add(query_id)
            evalue = float(parts[10])

            # Apply e-value filter
            if evalue > evalue_threshold:
                stats.hits_filtered_by_evalue += 1
                continue

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
                stitle=parts[12] if len(parts) > 12 else ""
            )

            results[hit.query_id].append(hit)

            # Track database source counts
            if hit.is_swissprot:
                stats.swissprot_hits += 1
            elif hit.is_trembl:
                stats.trembl_hits += 1
            else:
                stats.other_hits += 1

    # Update stats
    stats.queries_in_file = len(all_queries_in_file)
    stats.queries_after_evalue_filter = len(results)
    stats.hits_after_evalue_filter = stats.swissprot_hits + stats.trembl_hits + stats.other_hits

    # Sort and truncate hits for each query
    # When prefer_swissprot is True: SwissProt first, then by bitscore
    # Otherwise: just by bitscore
    for query_id in results:
        if prefer_swissprot:
            # Sort key: (not is_swissprot, -bitscore)
            # This puts SwissProt first (False < True), then sorts by bitscore descending
            results[query_id].sort(key=lambda h: (not h.is_swissprot, -h.bitscore))
        else:
            # Sort by bitscore only
            results[query_id].sort(key=lambda h: -h.bitscore)
        # Truncate to top_n
        results[query_id] = results[query_id][:top_n]

    logger.info(f"Parsed hits for {len(results)} queries")
    logger.info(f"  Total hits in file: {stats.total_hits_in_file}")
    if stats.hits_filtered_by_evalue > 0:
        logger.info(f"  Filtered by e-value (>{evalue_threshold}): {stats.hits_filtered_by_evalue}")
    logger.info(f"  Hits passing e-value filter: {stats.hits_after_evalue_filter}")
    logger.info(f"  SwissProt (sp|): {stats.swissprot_hits}")
    logger.info(f"  TrEMBL (tr|): {stats.trembl_hits}")
    if stats.other_hits > 0:
        logger.info(f"  Other: {stats.other_hits}")

    # Count retained hits after truncation
    retained_sp = sum(1 for hits in results.values() for h in hits if h.is_swissprot)
    retained_tr = sum(1 for hits in results.values() for h in hits if h.is_trembl)
    logger.info(f"  Retained after top-{top_n} selection: {retained_sp} SwissProt, {retained_tr} TrEMBL")

    return dict(results), stats


# =============================================================================
# CONSENSUS NAME FINDING
# =============================================================================

def normalize_gene_name(name: str) -> str:
    """
    Normalize a gene name for comparison purposes.
    """
    if not name:
        return ""

    # Lowercase
    name = name.lower()

    # Remove common suffixes/prefixes
    name = re.sub(r'\s*-like\s*', ' ', name)
    name = re.sub(r'\s*homolog\s*', ' ', name)
    name = re.sub(r'\s*isoform\s*\w*', '', name)
    name = re.sub(r'\s*variant\s*\w*', '', name)
    name = re.sub(r'\s*precursor\s*', '', name)
    name = re.sub(r'\s*fragment\s*', '', name)

    # Remove numbers at end (isoform numbers)
    name = re.sub(r'\s+\d+$', '', name)

    # Remove extra whitespace
    name = ' '.join(name.split())

    return name.strip()


def tokenize_name(name: str) -> set:
    """Convert a name to a set of tokens for similarity comparison."""
    # Remove punctuation and split
    tokens = re.findall(r'\b\w+\b', name.lower())
    # Remove very short tokens and common words
    stopwords = {'the', 'a', 'an', 'of', 'and', 'or', 'in', 'to', 'for', 'with', 'by'}
    return {t for t in tokens if len(t) > 2 and t not in stopwords}


def name_similarity(name1: str, name2: str) -> float:
    """
    Calculate similarity between two gene names using Jaccard index of tokens.
    """
    tokens1 = tokenize_name(name1)
    tokens2 = tokenize_name(name2)

    if not tokens1 or not tokens2:
        return 0.0

    intersection = len(tokens1 & tokens2)
    union = len(tokens1 | tokens2)

    return intersection / union if union > 0 else 0.0


def find_consensus_name(names: list, similarity_threshold: float = 0.5,
                        min_cluster_fraction: float = 0.4,
                        bitscores: list = None) -> tuple:
    """
    Find the consensus gene name from a list of names.

    Args:
        names: List of gene names from BLAST hits
        similarity_threshold: Minimum similarity to join a cluster
        min_cluster_fraction: Minimum fraction of hits in consensus cluster
        bitscores: Optional list of bitscores (same order as names) for weighting

    Returns: (consensus_name, indices_of_matching_hits)
    """
    if not names:
        return "", []

    if len(names) == 1:
        return names[0], [0]

    # Normalize names
    normalized = [normalize_gene_name(n) for n in names]

    # Build similarity clusters using simple greedy clustering
    clusters = []  # list of (representative_name, [indices])

    for i, name in enumerate(normalized):
        if not name:
            continue

        # Find best matching cluster
        best_cluster = None
        best_sim = 0

        for cluster_idx, (rep_name, indices) in enumerate(clusters):
            sim = name_similarity(name, rep_name)
            if sim > best_sim and sim >= similarity_threshold:
                best_sim = sim
                best_cluster = cluster_idx

        if best_cluster is not None:
            clusters[best_cluster][1].append(i)
        else:
            # Start new cluster
            clusters.append((name, [i]))

    if not clusters:
        return "", []

    # Score clusters - use bitscore weighting if available
    if bitscores and len(bitscores) == len(names):
        # Weight by sum of bitscores in cluster (better hits = more weight)
        def cluster_score(cluster):
            rep_name, indices = cluster
            total_score = sum(bitscores[i] for i in indices)
            # Also factor in cluster size to avoid single high-scoring outlier dominating
            size_factor = len(indices) ** 0.5  # sqrt to balance size vs score
            return total_score * size_factor

        clusters.sort(key=cluster_score, reverse=True)
    else:
        # Fall back to largest cluster
        clusters.sort(key=lambda x: len(x[1]), reverse=True)

    largest_cluster = clusters[0]

    # Check if it meets the minimum fraction requirement
    cluster_fraction = len(largest_cluster[1]) / len(names)

    if cluster_fraction >= min_cluster_fraction:
        # Return the original (non-normalized) name from the best hit in cluster
        # (first one, which has highest bitscore if sorted)
        best_idx = largest_cluster[1][0]
        return names[best_idx], largest_cluster[1]

    # No clear consensus
    return "", []


# =============================================================================
# GO TERM MAPPING
# =============================================================================

def collect_accessions_from_blast(blast_results: dict) -> set:
    """
    Collect all possible accession formats from parsed BLAST results.
    Used to filter GO mapping loading to only relevant entries.
    """
    accessions = set()
    for hits in blast_results.values():
        for hit in hits:
            accessions.update(extract_accession(hit.subject_id))
    return accessions


def load_go_mapping(mapping_file: str, accession_filter: set = None) -> dict:
    """
    Load accession to GO term mapping.

    Expected format (tab-separated):
      accession<TAB>GO:0000001;GO:0000002;...

    Or UniProt ID mapping format:
      accession<TAB>GO:0000001<TAB>GO:0000002<TAB>...

    Args:
        mapping_file: Path to the mapping file
        accession_filter: If provided, only load entries for these accessions.
                          This dramatically reduces memory for large mapping files.

    Returns: dict of accession -> set of GO IDs
    """
    logger.info(f"Loading GO mapping from {mapping_file}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter)} accessions from BLAST results")

    mapping = defaultdict(set)
    lines_scanned = 0
    lines_loaded = 0

    with open(mapping_file) as f:
        for line in f:
            lines_scanned += 1
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            accession = parts[0].strip()

            # Skip if not in filter set (when filter is provided)
            if accession_filter is not None and accession not in accession_filter:
                continue

            lines_loaded += 1

            # Handle both formats; use sys.intern to deduplicate GO term strings
            for part in parts[1:]:
                for go_term in re.split(r'[;,]', part):
                    go_term = go_term.strip()
                    if go_term.startswith('GO:'):
                        # Intern strings to avoid storing duplicate GO term objects
                        mapping[accession].add(sys.intern(go_term))

    logger.info(f"Loaded GO mappings for {len(mapping)} accessions")
    if accession_filter:
        logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return dict(mapping)


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
    # Just add the full ID and any recognizable parts
    for match in re.findall(r'[A-Z]{2}_\d+\.?\d*', subject_id):
        accessions.append(match)
        accessions.append(match.split('.')[0])  # Without version

    return accessions


def get_go_terms_for_hits(hits: list, go_mapping: dict) -> set:
    """
    Get all GO terms associated with a list of BLAST hits.
    """
    go_terms = set()

    for hit in hits:
        accessions = extract_accession(hit.subject_id)
        for acc in accessions:
            if acc in go_mapping:
                go_terms.update(go_mapping[acc])

    return go_terms


# =============================================================================
# NCBI CROSS-REFERENCE LOADING
# =============================================================================

def load_ncbi_idmapping(idmapping_file: str, accession_filter: set = None) -> dict:
    """
    Load UniProt accession to NCBI GeneID mapping from idmapping_selected.tab.

    File format (tab-separated):
      Column 1: UniProt accession (e.g., P31946)
      Column 3: NCBI GeneID (e.g., 7529)

    Args:
        idmapping_file: Path to idmapping_selected.tab
        accession_filter: If provided, only load entries for these accessions

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

            # Skip if no GeneID
            if not gene_id:
                continue

            # Skip if not in filter set (when filter is provided)
            if accession_filter is not None and accession not in accession_filter:
                continue

            lines_loaded += 1
            mapping[accession] = gene_id

    logger.info(f"Loaded {len(mapping):,} UniProt->GeneID mappings")
    logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return mapping


def load_ncbi_geneinfo(geneinfo_file: str, geneid_filter: set = None) -> dict:
    """
    Load NCBI GeneID to gene symbol mapping from gene_info file.

    File format (tab-separated):
      Column 2: GeneID
      Column 3: Symbol

    Args:
        geneinfo_file: Path to gene_info file
        geneid_filter: If provided, only load entries for these GeneIDs

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
            # Skip header
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

            # Skip placeholder entries
            if symbol == 'NEWENTRY' or symbol == '-':
                continue

            # Skip if not in filter set (when filter is provided)
            if geneid_filter is not None and gene_id not in geneid_filter:
                continue

            lines_loaded += 1
            mapping[gene_id] = symbol

    logger.info(f"Loaded {len(mapping):,} GeneID->Symbol mappings")
    logger.info(f"  Scanned {lines_scanned:,} lines, loaded {lines_loaded:,} entries")

    return mapping


def get_ncbi_symbol(accession: str, uniprot_to_geneid: dict, geneid_to_symbol: dict) -> str:
    """
    Look up gene symbol via NCBI cross-reference.

    Args:
        accession: UniProt accession (e.g., P31946)
        uniprot_to_geneid: Dict mapping UniProt accession -> GeneID
        geneid_to_symbol: Dict mapping GeneID -> Symbol

    Returns: Gene symbol if found, empty string otherwise
    """
    gene_id = uniprot_to_geneid.get(accession)
    if gene_id:
        return geneid_to_symbol.get(gene_id, "")
    return ""


# =============================================================================
# MAIN ANNOTATION PIPELINE
# =============================================================================

def annotate_queries(blast_results: dict,
                     go_mapping: dict,
                     go_hierarchy: GOHierarchy,
                     bad_name_filter: BadNameFilter,
                     all_gene_ids: list,
                     top_n: int = 25,
                     consensus_threshold: float = 0.5,
                     min_consensus_fraction: float = 0.4,
                     namespaces: list = None,
                     uniprot_to_geneid: dict = None,
                     geneid_to_symbol: dict = None) -> dict:
    """
    Main annotation pipeline for parsed BLAST results.

    Args:
        blast_results: dict of query_id -> list of BlastHit objects
        go_mapping: dict of accession -> set of GO IDs
        go_hierarchy: GOHierarchy object
        bad_name_filter: BadNameFilter object
        all_gene_ids: List of all gene IDs from the input FASTA file.
                      Annotations will be created for all IDs (including those without BLAST hits).
        top_n: Number of top hits to consider
        consensus_threshold: Similarity threshold for name clustering
        min_consensus_fraction: Minimum fraction of hits for consensus
        namespaces: GO namespaces to include
        uniprot_to_geneid: Optional dict mapping UniProt accession -> NCBI GeneID
        geneid_to_symbol: Optional dict mapping NCBI GeneID -> gene symbol

    Returns: dict of query_id -> QueryAnnotation
    """
    logger.info(f"Annotating {len(all_gene_ids):,} genes ({len(blast_results):,} have BLAST hits)...")

    annotations = {}

    stats = {
        'total_genes': len(all_gene_ids),
        'genes_with_blast_hits': 0,
        'genes_with_good_hits': 0,
        'genes_with_consensus': 0,
        'genes_with_go_terms': 0,
        'total_hits': 0,
        'hits_after_bad_name_filter': 0,
        'hits_after_consensus_filter': 0,
        'symbols_from_uniprot': 0,
        'symbols_from_ncbi': 0,
    }

    # Check if NCBI cross-reference is available
    ncbi_available = uniprot_to_geneid is not None and geneid_to_symbol is not None

    for query_id in all_gene_ids:
        # Get BLAST hits for this query (may be empty)
        hits = blast_results.get(query_id, [])

        # Sort hits: Swiss-Prot first, then by bitscore within each category
        sorted_hits = sorted(hits, key=lambda h: (not h.is_swissprot, -h.bitscore))
        annot = QueryAnnotation(query_id=query_id, hits=sorted_hits[:top_n])

        # Track genes with BLAST hits
        if hits:
            stats['genes_with_blast_hits'] += 1
        stats['total_hits'] += len(annot.hits)

        # Step 1: Filter out bad names
        good_hits = []
        for hit in annot.hits:
            if not bad_name_filter.is_bad_name(hit.gene_name):
                good_hits.append(hit)

        stats['hits_after_bad_name_filter'] += len(good_hits)

        if not good_hits:
            annotations[query_id] = annot
            continue

        stats['genes_with_good_hits'] += 1

        # Step 2: Find consensus name and filter
        names = [hit.gene_name for hit in good_hits]
        bitscores = [hit.bitscore for hit in good_hits]
        consensus_name, matching_indices = find_consensus_name(
            names,
            similarity_threshold=consensus_threshold,
            min_cluster_fraction=min_consensus_fraction,
            bitscores=bitscores
        )

        if consensus_name and matching_indices:
            annot.consensus_name = consensus_name
            annot.filtered_hits = [good_hits[i] for i in matching_indices]
            stats['genes_with_consensus'] += 1
        else:
            # No clear consensus - keep all good hits but no consensus name
            annot.filtered_hits = good_hits
            # Use the top hit's name as a fallback
            annot.consensus_name = good_hits[0].gene_name if good_hits else ""

        stats['hits_after_consensus_filter'] += len(annot.filtered_hits)

        # Step 2.5: Find consensus gene symbol from filtered hits
        # Priority: use symbol from the same hit that provided the consensus name
        # Fallback: use most common valid symbol from the cluster

        # First, try to get symbol from the best hit (the one that provided consensus_name)
        best_hit_symbol = None
        if annot.filtered_hits:
            best_hit = annot.filtered_hits[0]
            sym = best_hit.gene_symbol
            if sym and best_hit.is_valid_gene_symbol(sym):
                best_hit_symbol = sym

        if best_hit_symbol:
            # Use symbol from the same hit as the protein name (keeps them coupled)
            annot.consensus_symbol = best_hit_symbol.lower()
            stats['symbols_from_uniprot'] += 1
        else:
            # Fallback 1: best hit lacks GN= field, use most common symbol from cluster
            valid_symbols = []
            for hit in annot.filtered_hits:
                sym = hit.gene_symbol
                if sym and hit.is_valid_gene_symbol(sym):
                    valid_symbols.append(sym)

            if valid_symbols:
                # Find most common symbol (assign even if only 1 exists)
                symbol_counts = Counter(valid_symbols)
                most_common_symbol, count = symbol_counts.most_common(1)[0]
                annot.consensus_symbol = most_common_symbol.lower()
                stats['symbols_from_uniprot'] += 1

            # Fallback 2: No UniProt symbol found, try NCBI cross-reference
            elif ncbi_available and annot.filtered_hits:
                # Try to get symbol from NCBI for the best hit
                best_hit = annot.filtered_hits[0]
                # Extract UniProt accession from subject_id
                accessions = extract_accession(best_hit.subject_id)
                for acc in accessions:
                    ncbi_symbol = get_ncbi_symbol(acc, uniprot_to_geneid, geneid_to_symbol)
                    # Validate NCBI symbol using same filtering as UniProt symbols
                    if ncbi_symbol and best_hit.is_valid_gene_symbol(ncbi_symbol):
                        annot.consensus_symbol = ncbi_symbol.lower()
                        stats['symbols_from_ncbi'] += 1
                        break

        # Step 3: Get GO terms from top 5 highest bitscore hits in consensus cluster
        top_hits_for_go = sorted(annot.filtered_hits, key=lambda h: h.bitscore, reverse=True)[:5]
        go_terms = get_go_terms_for_hits(top_hits_for_go, go_mapping)

        # Filter by namespace if specified
        if namespaces and go_terms:
            go_terms = go_hierarchy.filter_by_namespace(go_terms, namespaces)

        annot.go_terms = go_terms

        if go_terms:
            stats['genes_with_go_terms'] += 1

            # Step 4: Filter to most specific terms
            annot.specific_go_terms = go_hierarchy.filter_to_specific(go_terms)

        annotations[query_id] = annot

    # Log statistics
    logger.info("Annotation statistics:")
    logger.info(f"  Total genes: {stats['total_genes']:,}")
    logger.info(f"  Genes with BLAST hits: {stats['genes_with_blast_hits']:,} "
                f"({100*stats['genes_with_blast_hits']/stats['total_genes']:.1f}%)")
    logger.info(f"  Genes with good hits (after bad name filter): {stats['genes_with_good_hits']:,} "
                f"({100*stats['genes_with_good_hits']/stats['total_genes']:.1f}%)")
    logger.info(f"  Genes with consensus name: {stats['genes_with_consensus']:,} "
                f"({100*stats['genes_with_consensus']/stats['total_genes']:.1f}%)")
    total_symbols = stats['symbols_from_uniprot'] + stats['symbols_from_ncbi']
    logger.info(f"  Genes with gene symbol: {total_symbols:,} "
                f"({100*total_symbols/stats['total_genes']:.1f}%)")
    if stats['symbols_from_uniprot'] > 0:
        logger.info(f"    - From UniProt GN= field: {stats['symbols_from_uniprot']:,}")
    if stats['symbols_from_ncbi'] > 0:
        logger.info(f"    - From NCBI cross-reference: {stats['symbols_from_ncbi']:,}")
    logger.info(f"  Genes with GO terms: {stats['genes_with_go_terms']:,} "
                f"({100*stats['genes_with_go_terms']/stats['total_genes']:.1f}%)")
    logger.info(f"  Total hits: {stats['total_hits']:,}")
    logger.info(f"  Hits after bad name filter: {stats['hits_after_bad_name_filter']:,}")
    logger.info(f"  Hits after consensus filter: {stats['hits_after_consensus_filter']:,}")

    return annotations


# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

def strip_uniprot_prefix(protein_name: str) -> str:
    """
    Strip UniProt ID prefix from protein name if present.

    UniProt format: sp|ACCESSION|ENTRY_NAME Description
                    tr|ACCESSION|ENTRY_NAME Description

    Examples:
        sp|Q02084|A33_PLEWA Zinc-binding protein A33 -> Zinc-binding protein A33
        tr|A0A123|A0A123_HUMAN Some protein -> Some protein

    Args:
        protein_name: The protein name potentially containing UniProt prefix

    Returns:
        The protein description without the UniProt ID prefix
    """
    if not protein_name:
        return protein_name

    # Match UniProt prefix pattern: sp|XXX|XXX_XXX or tr|XXX|XXX_XXX followed by space
    match = re.match(r'^(?:sp|tr)\|[A-Za-z0-9]+\|[A-Za-z0-9_]+\s+(.+)$', protein_name)
    if match:
        return match.group(1)

    return protein_name


def write_annotated_fasta(input_fasta: str, annotations: dict, output_file: Path):
    """
    Write an annotated FASTA file with Name= and product= fields added to headers.
    Replaces any existing Name= or product= fields in the header.

    Args:
        input_fasta: Path to the original FASTA file
        annotations: dict of query_id -> QueryAnnotation
        output_file: Path to write the annotated FASTA
    """
    logger.info(f"Writing annotated FASTA to {output_file}")

    annotated_count = 0
    total_count = 0

    with open(input_fasta) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                total_count += 1
                header = line[1:].rstrip()
                # Extract gene ID (first field)
                gene_id = header.split()[0] if header else ""

                # Look up annotation
                annot = annotations.get(gene_id)

                if annot and (annot.consensus_symbol or annot.consensus_name):
                    annotated_count += 1

                    # Remove any existing Name= or product= fields from header
                    header_cleaned = re.sub(r'\s+Name=\S+', '', header)
                    header_cleaned = re.sub(r'\s+product=[^\s]+(?:\s+[^\s=]+)*(?=\s+\w+=|$)', '', header_cleaned)
                    # Simpler approach: remove Name=value and product=value patterns
                    header_cleaned = re.sub(r'\bName=\S+\s*', '', header_cleaned)
                    header_cleaned = re.sub(r'\bproduct=\S+\s*', '', header_cleaned)
                    header_cleaned = header_cleaned.strip()

                    # Build annotation fields
                    fields = []
                    if annot.consensus_symbol:
                        fields.append(f"Name={annot.consensus_symbol}")
                    if annot.consensus_name:
                        # Strip UniProt ID prefix from protein name
                        product_name = strip_uniprot_prefix(annot.consensus_name)
                        fields.append(f"product={product_name}")

                    # Append to cleaned header
                    new_header = f">{header_cleaned} {' '.join(fields)}\n"
                    fout.write(new_header)
                else:
                    fout.write(line)
            else:
                fout.write(line)

    logger.info(f"Annotated {annotated_count:,} of {total_count:,} sequences in FASTA")


def write_annotated_gff3(input_gff: str, annotations: dict, output_file: Path):
    """
    Write an annotated GFF3 file with Name= and product= attributes added.

    Args:
        input_gff: Path to the original GFF3 file
        annotations: dict of query_id -> QueryAnnotation
        output_file: Path to write the annotated GFF3
    """
    logger.info(f"Writing annotated GFF3 to {output_file}")

    annotated_count = 0
    feature_count = 0

    with open(input_gff) as fin, open(output_file, 'w') as fout:
        for line in fin:
            # Pass through comments and directives unchanged
            if line.startswith('#'):
                fout.write(line)
                continue

            line = line.rstrip()
            if not line:
                fout.write('\n')
                continue

            parts = line.split('\t')
            if len(parts) != 9:
                # Not a valid GFF3 feature line, pass through
                fout.write(line + '\n')
                continue

            feature_count += 1
            attributes = parts[8]

            # Extract ID from attributes
            gene_id = None
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match:
                gene_id = id_match.group(1)

            # Look up annotation
            annot = annotations.get(gene_id) if gene_id else None

            if annot and (annot.consensus_symbol or annot.consensus_name):
                annotated_count += 1
                # Build new attributes
                new_attrs = []

                # Parse existing attributes
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    # Skip existing Name and product attributes (we'll replace them)
                    if attr.startswith('Name=') or attr.startswith('product='):
                        continue
                    new_attrs.append(attr)

                # Add our annotation attributes
                if annot.consensus_symbol:
                    new_attrs.append(f"Name={annot.consensus_symbol}")
                if annot.consensus_name:
                    # Strip UniProt ID prefix and replace spaces with underscores for GFF3
                    product = strip_uniprot_prefix(annot.consensus_name)
                    product = product.replace(' ', '_')
                    new_attrs.append(f"product={product}")

                parts[8] = ';'.join(new_attrs)

            fout.write('\t'.join(parts) + '\n')

    logger.info(f"Annotated {annotated_count:,} of {feature_count:,} features in GFF3")


def write_annotation_results(annotations: dict,
                             go_hierarchy: GOHierarchy,
                             output_dir: Path,
                             prefix: str,
                             input_fasta: str = None,
                             input_gff: str = None,
                             args: argparse.Namespace = None,
                             blast_stats: BlastStats = None,
                             total_sequences: int = None):
    """
    Write annotation results to various output files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Main annotation table
    main_file = output_dir / f"{prefix}_annotations.tsv"
    logger.info(f"Writing main annotations to {main_file}")

    with open(main_file, 'w') as f:
        f.write("query_id\tgene_symbol\tconsensus_name\tnum_hits\tnum_filtered_hits\t"
                "num_go_terms\tnum_specific_go_terms\ttop_hit_evalue\ttop_hit_pident\ttop_hit_bitscore\n")

        for query_id, annot in sorted(annotations.items()):
            top_evalue = annot.filtered_hits[0].evalue if annot.filtered_hits else "NA"
            top_pident = annot.filtered_hits[0].pident if annot.filtered_hits else "NA"
            top_bitscore = annot.filtered_hits[0].bitscore if annot.filtered_hits else "NA"

            f.write(f"{query_id}\t{annot.consensus_symbol}\t{annot.consensus_name}\t{len(annot.hits)}\t"
                    f"{len(annot.filtered_hits)}\t{len(annot.go_terms)}\t"
                    f"{len(annot.specific_go_terms)}\t{top_evalue}\t{top_pident}\t{top_bitscore}\n")

    # 2. GO term assignments (all GO terms)
    go_all_file = output_dir / f"{prefix}_GO_all.tsv"
    logger.info(f"Writing all GO terms to {go_all_file}")

    with open(go_all_file, 'w') as f:
        f.write("query_id\tgene_symbol\tprotein_name\tgo_terms\n")
        for query_id, annot in sorted(annotations.items()):
            gene_symbol = annot.consensus_symbol or ""
            protein_name = annot.consensus_name or ""
            go_terms_str = ';'.join(sorted(annot.go_terms)) if annot.go_terms else ""
            f.write(f"{query_id}\t{gene_symbol}\t{protein_name}\t{go_terms_str}\n")

    # 3. GO term assignments (specific/leaf terms only)
    go_specific_file = output_dir / f"{prefix}_GO_specific.tsv"
    logger.info(f"Writing specific GO terms to {go_specific_file}")

    with open(go_specific_file, 'w') as f:
        f.write("query_id\tgene_symbol\tprotein_name\tgo_terms\n")
        for query_id, annot in sorted(annotations.items()):
            gene_symbol = annot.consensus_symbol or ""
            protein_name = annot.consensus_name or ""
            go_terms_str = ';'.join(sorted(annot.specific_go_terms)) if annot.specific_go_terms else ""
            f.write(f"{query_id}\t{gene_symbol}\t{protein_name}\t{go_terms_str}\n")

    # 4. Detailed GO term table with names
    go_detailed_file = output_dir / f"{prefix}_GO_detailed.tsv"
    logger.info(f"Writing detailed GO terms to {go_detailed_file}")

    with open(go_detailed_file, 'w') as f:
        f.write("query_id\tgo_term\tgo_name\tnamespace\tis_specific\n")
        for query_id, annot in sorted(annotations.items()):
            if annot.go_terms:
                for go_term in sorted(annot.go_terms):
                    go_name = go_hierarchy.names.get(go_term, "unknown")
                    namespace = go_hierarchy.namespace.get(go_term, "unknown")
                    is_specific = "yes" if go_term in annot.specific_go_terms else "no"
                    f.write(f"{query_id}\t{go_term}\t{go_name}\t{namespace}\t{is_specific}\n")
            else:
                # Include genes without GO terms with empty values
                f.write(f"{query_id}\t\t\t\t\n")

    # 5. Protein names table
    names_file = output_dir / f"{prefix}_protein_names.tsv"
    logger.info(f"Writing protein names to {names_file}")

    with open(names_file, 'w') as f:
        f.write("query_id\tgene_symbol\tprotein_name\tall_symbols\tall_hit_names\n")
        for query_id, annot in sorted(annotations.items()):
            all_symbols = "|".join(h.gene_symbol for h in annot.filtered_hits[:5] if h.gene_symbol)
            all_names = "|".join(h.gene_name for h in annot.filtered_hits[:5])
            f.write(f"{query_id}\t{annot.consensus_symbol}\t{annot.consensus_name}\t{all_symbols}\t{all_names}\n")

    # 6. Annotated FASTA file
    if input_fasta:
        fasta_out = output_dir / f"{prefix}_annotated.fasta"
        write_annotated_fasta(input_fasta, annotations, fasta_out)

    # 7. Annotated GFF3 file (if GFF input was provided)
    if input_gff:
        gff_out = output_dir / f"{prefix}_annotated.gff3"
        write_annotated_gff3(input_gff, annotations, gff_out)

    # 8. Summary statistics
    stats_file = output_dir / f"{prefix}_summary.txt"
    logger.info(f"Writing summary to {stats_file}")

    total = len(annotations)
    with_blast_hits = sum(1 for a in annotations.values() if a.hits)
    with_protein_name = sum(1 for a in annotations.values() if a.consensus_name)
    with_gene_symbol = sum(1 for a in annotations.values() if a.consensus_symbol)
    with_go = sum(1 for a in annotations.values() if a.go_terms)
    with_specific = sum(1 for a in annotations.values() if a.specific_go_terms)

    all_go = set()
    all_specific = set()
    for annot in annotations.values():
        all_go.update(annot.go_terms)
        all_specific.update(annot.specific_go_terms)

    with open(stats_file, 'w') as f:
        f.write("GOAnnotate Summary\n")
        f.write("=" * 70 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Command section
        if args:
            f.write("Command:\n")
            f.write("-" * 70 + "\n")
            cmd_parts = ["GOAnnotate.py"]
            cmd_parts.extend(["--transcripts", str(args.transcripts)])
            cmd_parts.extend(["--blast-results", str(args.blast_results)])
            cmd_parts.extend(["--go-mapping", str(args.go_mapping)])
            cmd_parts.extend(["--go-obo", str(args.go_obo)])
            cmd_parts.extend(["--bad-names", str(args.bad_names)])
            cmd_parts.extend(["--output", str(args.output)])
            cmd_parts.extend(["--prefix", str(args.prefix)])
            cmd_parts.extend(["--evalue", str(args.evalue)])
            cmd_parts.extend(["--top-n", str(args.top_n)])
            cmd_parts.extend(["--consensus-threshold", str(args.consensus_threshold)])
            cmd_parts.extend(["--min-consensus-fraction", str(args.min_consensus_fraction)])
            cmd_parts.extend(["--namespace"] + args.namespace)
            if args.no_prefer_swissprot:
                cmd_parts.append("--no-prefer-swissprot")
            if args.gff:
                cmd_parts.extend(["--gff", str(args.gff)])
            f.write(" \\\n    ".join(cmd_parts) + "\n\n")

        # Input files section
        if args:
            f.write("Input Files:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Transcripts:    {args.transcripts}\n")
            f.write(f"  BLAST results:  {args.blast_results}\n")
            f.write(f"  GO mapping:     {args.go_mapping}\n")
            f.write(f"  GO OBO:         {args.go_obo}\n")
            f.write(f"  Bad names:      {args.bad_names}\n")
            if args.gff:
                f.write(f"  GFF3:           {args.gff}\n")
            f.write("\n")

        # Parameters section
        if args:
            f.write("Parameters:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  E-value threshold:        {args.evalue}\n")
            f.write(f"  Top N hits:               {args.top_n}\n")
            f.write(f"  Consensus threshold:      {args.consensus_threshold}\n")
            f.write(f"  Min consensus fraction:   {args.min_consensus_fraction}\n")
            f.write(f"  Prefer SwissProt:         {not args.no_prefer_swissprot}\n")
            f.write(f"  GO namespaces:            {', '.join(args.namespace)}\n\n")

        # BLAST statistics section
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
            if queries_lost > 0:
                f.write(f"  Queries lost to e-value:      {queries_lost:,} "
                        f"({100*queries_lost/blast_stats.queries_in_file:.1f}%)\n")
            f.write(f"  SwissProt hits:               {blast_stats.swissprot_hits:,}\n")
            f.write(f"  TrEMBL hits:                  {blast_stats.trembl_hits:,}\n")
            if blast_stats.other_hits > 0:
                f.write(f"  Other hits:                   {blast_stats.other_hits:,}\n")
            f.write("\n")

        # Annotation coverage section
        f.write("Annotation Results:\n")
        f.write("-" * 70 + "\n")
        if total_sequences is not None:
            f.write(f"  Total Sequences:              {total_sequences:,}\n")
        f.write(f"  Total Genes:                  {total:,}\n")
        f.write(f"  Genes with BLAST hits:        {with_blast_hits:,} ({100*with_blast_hits/total:.1f}%)\n")
        f.write(f"  Genes with protein name:      {with_protein_name:,} ({100*with_protein_name/total:.1f}%)\n")
        f.write(f"  Genes with gene symbol:       {with_gene_symbol:,} ({100*with_gene_symbol/total:.1f}%)\n")
        f.write(f"  Genes with GO terms:          {with_go:,} ({100*with_go/total:.1f}%)\n\n")

        # GO term statistics
        f.write("GO Term Statistics:\n")
        f.write("-" * 70 + "\n")
        f.write(f"  Total unique GO terms (all):      {len(all_go):,}\n")
        f.write(f"  Total unique GO terms (specific): {len(all_specific):,}\n")
        f.write(f"  GO term reduction: {len(all_go):,} -> {len(all_specific):,} "
                f"({100*(1-len(all_specific)/len(all_go)) if all_go else 0:.1f}% reduction)\n")

    logger.info(f"Results written to {output_dir}")


# =============================================================================
# CLI
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="GOAnnotate: Gene Ontology annotation pipeline using UniProt BLAST/Diamond results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input options
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument(
        "--transcripts", required=True, metavar="FILE",
        help="Input FASTA file containing all transcripts (used to define the complete gene set)"
    )
    input_group.add_argument(
        "--blast-results", "--blast", required=True, metavar="FILE",
        help="BLAST/Diamond results file (format 6 with stitle)"
    )
    input_group.add_argument(
        "--gff", metavar="FILE",
        help="Optional GFF3 file to annotate with gene symbols and product names"
    )

    # GO mapping
    go_group = parser.add_argument_group("GO options")
    go_group.add_argument(
        "--go-mapping", required=True, metavar="FILE",
        help="Accession to GO term mapping file (accession<TAB>GO:xxxx;GO:yyyy)"
    )
    go_group.add_argument(
        "--go-obo", required=True, metavar="FILE",
        help="GO OBO hierarchy file"
    )
    go_group.add_argument(
        "--namespace", nargs='+', default=["BP", "MF", "CC"],
        choices=["BP", "MF", "CC"],
        help="GO namespaces to include"
    )

    # NCBI cross-reference options
    ncbi_group = parser.add_argument_group("NCBI cross-reference options (for improved gene symbol coverage)")
    ncbi_group.add_argument(
        "--ncbi-idmapping", metavar="FILE",
        help="UniProt ID mapping file (idmapping_selected.tab) for UniProt->NCBI GeneID lookup"
    )
    ncbi_group.add_argument(
        "--ncbi-geneinfo", metavar="FILE",
        help="NCBI gene_info file for GeneID->Symbol lookup"
    )

    # Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold for filtering hits (use to apply stricter filtering than BLAST/Diamond)"
    )
    filter_group.add_argument(
        "--top-n", type=int, default=20,
        help="Number of top BLAST hits to consider per query"
    )
    filter_group.add_argument(
        "--bad-names", required=True, metavar="FILE",
        help="Bad names pattern file for filtering uninformative protein names"
    )
    filter_group.add_argument(
        "--consensus-threshold", type=float, default=0.5,
        help="Similarity threshold for name clustering (0-1)"
    )
    filter_group.add_argument(
        "--min-consensus-fraction", type=float, default=0.4,
        help="Minimum fraction of hits that must agree for consensus"
    )
    filter_group.add_argument(
        "--no-prefer-swissprot", action="store_true",
        help="Don't prioritize SwissProt hits over TrEMBL (by default, SwissProt is preferred)"
    )

    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output", "-o", default="./GOAnnotate_results", metavar="DIR",
        help="Output directory"
    )
    output_group.add_argument(
        "--prefix", default="annotation",
        help="Output file prefix"
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Initialize components
    logger.info("Initializing GOAnnotate pipeline...")

    # Validate input files exist
    if not Path(args.transcripts).exists():
        logger.error(f"Transcripts FASTA file not found: {args.transcripts}")
        sys.exit(1)

    if not Path(args.blast_results).exists():
        logger.error(f"BLAST results file not found: {args.blast_results}")
        sys.exit(1)

    if not Path(args.go_obo).exists():
        logger.error(f"GO OBO file not found: {args.go_obo}")
        sys.exit(1)

    if not Path(args.go_mapping).exists():
        logger.error(f"GO mapping file not found: {args.go_mapping}")
        sys.exit(1)

    if args.gff and not Path(args.gff).exists():
        logger.error(f"GFF3 file not found: {args.gff}")
        sys.exit(1)

    # Validate NCBI files if provided (both must be provided together)
    if args.ncbi_idmapping and not args.ncbi_geneinfo:
        logger.error("--ncbi-geneinfo is required when using --ncbi-idmapping")
        sys.exit(1)
    if args.ncbi_geneinfo and not args.ncbi_idmapping:
        logger.error("--ncbi-idmapping is required when using --ncbi-geneinfo")
        sys.exit(1)
    if args.ncbi_idmapping and not Path(args.ncbi_idmapping).exists():
        logger.error(f"NCBI ID mapping file not found: {args.ncbi_idmapping}")
        sys.exit(1)
    if args.ncbi_geneinfo and not Path(args.ncbi_geneinfo).exists():
        logger.error(f"NCBI gene_info file not found: {args.ncbi_geneinfo}")
        sys.exit(1)

    # Step 0: Parse transcripts file to get complete gene list
    total_sequences, all_gene_ids = parse_fasta_headers(args.transcripts)

    # Step 1: Parse BLAST results first (needed to filter GO mapping loading)
    blast_results, blast_stats = parse_blast_results(
        args.blast_results,
        top_n=args.top_n,
        evalue_threshold=args.evalue,
        prefer_swissprot=not args.no_prefer_swissprot
    )

    # Step 2: Collect accessions from BLAST hits for filtering GO mapping
    blast_accessions = collect_accessions_from_blast(blast_results)

    # Step 3: Load GO mapping (filtered to only accessions in BLAST results)
    go_mapping = load_go_mapping(args.go_mapping, accession_filter=blast_accessions)

    # Load GO hierarchy
    go_hierarchy = GOHierarchy(args.go_obo)

    # Initialize bad name filter
    bad_name_filter = BadNameFilter(args.bad_names)

    # Step 3.5: Load NCBI cross-reference mappings if provided
    uniprot_to_geneid = None
    geneid_to_symbol = None
    if args.ncbi_idmapping and args.ncbi_geneinfo:
        # Load UniProt -> GeneID mapping (filtered to BLAST accessions)
        uniprot_to_geneid = load_ncbi_idmapping(
            args.ncbi_idmapping,
            accession_filter=blast_accessions
        )
        # Collect GeneIDs that were found
        found_geneids = set(uniprot_to_geneid.values())
        # Load GeneID -> Symbol mapping (filtered to found GeneIDs)
        geneid_to_symbol = load_ncbi_geneinfo(
            args.ncbi_geneinfo,
            geneid_filter=found_geneids
        )

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 4: Annotate
    annotations = annotate_queries(
        blast_results=blast_results,
        go_mapping=go_mapping,
        go_hierarchy=go_hierarchy,
        bad_name_filter=bad_name_filter,
        all_gene_ids=all_gene_ids,
        top_n=args.top_n,
        consensus_threshold=args.consensus_threshold,
        min_consensus_fraction=args.min_consensus_fraction,
        namespaces=args.namespace,
        uniprot_to_geneid=uniprot_to_geneid,
        geneid_to_symbol=geneid_to_symbol
    )

    # Step 5: Write output
    write_annotation_results(
        annotations=annotations,
        go_hierarchy=go_hierarchy,
        output_dir=output_dir,
        prefix=args.prefix,
        input_fasta=args.transcripts,
        input_gff=args.gff,
        args=args,
        blast_stats=blast_stats,
        total_sequences=total_sequences
    )

    logger.info("GOAnnotate pipeline completed successfully")


if __name__ == "__main__":
    main()
