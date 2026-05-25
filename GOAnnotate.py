#!/usr/bin/env python3
"""
GOAnnotate: Gene Ontology Annotation Pipeline using BLAST/Diamond and UniProt
Written by E. J. Bentz

Annotation-only mode: annotates pre-computed BLAST results against a
transcript-level CDS FASTA with transcript-to-gene mapping.

Required: --transcripts, --blast-results, --go-mapping (or --db), --bad-names

The --transcripts CDS FASTA defines the complete gene set. Transcript IDs
are extracted from the first field of each header; gene IDs from the gene=
field (if absent, transcript ID is used as gene ID). All genes appear in
the output, including those without BLAST hits.

The GO OBO hierarchy file (--go-obo) is automatically downloaded if not
provided or if the existing file is more than 30 days old.

Usage:
  ./GOAnnotate.py --transcripts CDS.fasta --blast-results diamond_results.tsv \\
                  --go-mapping GO_mapping.tsv --bad-names bad_names.txt --go-obo go.obo -o output_dir

"""
import sys
import os
import re
import argparse
import logging
import math
import multiprocessing
import time
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional
from difflib import SequenceMatcher
import statistics
import shutil
import urllib.parse

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


####################
# BAD NAMES FILTER
####################
class BadNameFilter:
    """
    Filter for cleaning and removing uninformative gene/protein names and symbols.

    Supports five types of patterns:
      1. Strip patterns (prefixed with 'strip_regex:') - artifacts stripped from names
      2. Exact matches (case-insensitive) - plain text patterns <=6 chars
      3. Substring matches (case-insensitive) - plain text patterns >6 chars
      4. Regex patterns (prefixed with 'regex:') - matched against protein names
      5. Symbol regex patterns (prefixed with 'symbol_regex:') - matched against gene symbols

    Processing order: clean_name() applies strip patterns first, then
    is_bad_name() checks the cleaned result against discard patterns.
    """

    def __init__(self, patterns_file: str):
        self.exact_matches = set()
        self.substrings = []
        self.regex_patterns = []
        self.symbol_regex_patterns = []
        self.strip_patterns = []

        if not Path(patterns_file).exists():
            raise FileNotFoundError(f"Bad names pattern file not found: {patterns_file}")

        with open(patterns_file) as f:
            self._parse_patterns(f.read())

        logger.info(
            f"BadNameFilter: loaded {len(self.strip_patterns)} strip, "
            f"{len(self.exact_matches)} exact, "
            f"{len(self.substrings)} substring, {len(self.regex_patterns)} regex, "
            f"{len(self.symbol_regex_patterns)} symbol_regex patterns "
            f"from '{patterns_file}'"
        )

    def _parse_patterns(self, text: str):
        """Parse pattern text into strip, exact, substring, regex, and symbol_regex categories."""
        for line in text.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('strip_regex:'):
                pattern = line[12:].strip()
                try:
                    self.strip_patterns.append(re.compile(pattern, re.IGNORECASE))
                except re.error as e:
                    logger.warning(f"Invalid strip_regex pattern '{pattern}': {e}")
            elif line.startswith('symbol_regex:'):
                pattern = line[13:].strip()
                try:
                    self.symbol_regex_patterns.append(re.compile(pattern, re.IGNORECASE))
                except re.error as e:
                    logger.warning(f"Invalid symbol_regex pattern '{pattern}': {e}")
            elif line.startswith('regex:'):
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
        """Check if a protein name matches any bad pattern."""
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

    def clean_name(self, name: str) -> Optional[str]:
        """
        Strip database-specific artifacts from a protein name.

        Applies all strip_regex patterns loaded from the patterns file,
        then cleans up residual punctuation and whitespace.

        This method runs BEFORE the bad_names discard filter so that
        informative names obscured by database prefixes, suffixes, or
        qualifiers are preserved rather than discarded.

        Args:
            name: Raw protein name string from a BLAST stitle.

        Returns:
            Cleaned protein name, or None if the name is empty after cleaning.
        """
        if not name:
            return None

        cleaned = name
        for pattern in self.strip_patterns:
            cleaned = pattern.sub('', cleaned)

        # Remove trailing commas left after stripping
        cleaned = re.sub(r',\s*$', '', cleaned)
        # Collapse multiple spaces into one and strip leading/trailing whitespace
        cleaned = re.sub(r'\s+', ' ', cleaned).strip()

        if not cleaned:
            return None

        return cleaned

    def is_bad_symbol(self, symbol: str) -> bool:
        """Check if a gene symbol matches any bad symbol pattern."""
        if not symbol:
            return True

        for pattern in self.symbol_regex_patterns:
            if pattern.search(symbol):
                return True

        return False

    def filter_names(self, names: list) -> list:
        """Return only names that are NOT bad."""
        return [n for n in names if not self.is_bad_name(n)]


###########################
# BLAST HIT DATA STRUCTURE
###########################

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
    name_concordance: float = 0.0
    symbol_concordance: float = 0.0
    n_isoforms: int = 0
    n_cluster_hits: int = 0
    winning_cluster_names: list = field(default_factory=list)
    winning_cluster_symbols: list = field(default_factory=list)
    mean_cluster_bitscore: float = 0.0
    go_source: str = ""


#########################
# GO HIERARCHY STRUCTURE
#########################

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


########################
# FASTA HEADER PARSING
########################

def parse_fasta_headers(fasta_file: str) -> tuple:
    """
    Parse CDS FASTA headers to build transcript-to-gene mapping.

    Expects headers where transcript/mRNA ID is the first field and
    gene ID is in a gene= field. If no gene= field exists, the
    transcript ID is used as the gene ID (1:1 mapping).

    Returns:
        Tuple of (total_sequences,
                  all_gene_ids: list of unique gene IDs in first-seen order,
                  transcript_to_gene: dict of transcript_id -> gene_id,
                  gene_to_transcripts: dict of gene_id -> list of transcript_ids)
    """
    logger.info(f"Parsing CDS FASTA headers from {fasta_file}")

    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(list)
    gene_ids_ordered = []
    gene_ids_seen = set()
    total_sequences = 0

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                total_sequences += 1
                header = line[1:].strip()
                fields = header.split()
                transcript_id = fields[0] if fields else ""

                # Scan remaining fields for gene=XXXX
                gene_id = None
                for fld in fields[1:]:
                    if fld.startswith('gene='):
                        gene_id = fld[5:]
                        break

                # If no gene= field, use transcript ID as gene ID
                if gene_id is None:
                    gene_id = transcript_id

                transcript_to_gene[transcript_id] = gene_id
                gene_to_transcripts[gene_id].append(transcript_id)

                if gene_id not in gene_ids_seen:
                    gene_ids_ordered.append(gene_id)
                    gene_ids_seen.add(gene_id)

    logger.info(f"Found {total_sequences:,} sequences, "
                f"{len(gene_ids_ordered):,} unique genes, "
                f"{len(transcript_to_gene):,} transcripts")

    n_multi = sum(1 for txs in gene_to_transcripts.values() if len(txs) > 1)
    if n_multi > 0:
        logger.info(f"  {n_multi:,} genes have multiple transcripts")

    return (total_sequences, gene_ids_ordered,
            transcript_to_gene, dict(gene_to_transcripts))


####################
# GO OBO DOWNLOAD
####################

GO_OBO_URL = "https://purl.obolibrary.org/obo/go.obo"
GO_OBO_MAX_AGE_DAYS = 30


def download_go_obo(output_path: str):
    """
    Download the current GO OBO hierarchy file from the Gene Ontology.

    Args:
        output_path: Path to write the downloaded go.obo file.
    """
    import urllib.request
    import urllib.error

    logger.info(f"Downloading GO OBO file from {GO_OBO_URL} ...")

    try:
        req = urllib.request.Request(GO_OBO_URL, headers={"User-Agent": "GOAnnotate/2.0"})
        with urllib.request.urlopen(req) as response, open(output_path, 'wb') as f:
            f.write(response.read())
    except urllib.error.URLError as e:
        logger.error(f"Failed to download GO OBO file: {e}")
        sys.exit(1)

    logger.info(f"GO OBO file saved to: {output_path}")


def resolve_go_obo(go_obo_path: str) -> str:
    """
    Ensure a valid, up-to-date GO OBO file is available.

    If a path is provided and the file exists and is less than 30 days old,
    it is used as-is. If the file is older than 30 days, a fresh copy is
    re-downloaded to the same path. If no path is provided, defaults to
    ./go.obo in the current working directory.

    Args:
        go_obo_path: User-supplied path to a go.obo file, or None.

    Returns:
        Path to the go.obo file to use.
    """
    # Default to ./go.obo when no path provided
    if not go_obo_path:
        go_obo_path = "go.obo"

    obo = Path(go_obo_path)

    if obo.exists():
        age_days = (datetime.now() - datetime.fromtimestamp(
            obo.stat().st_mtime)).days
        if age_days <= GO_OBO_MAX_AGE_DAYS:
            logger.info(f"Using existing GO OBO file: {go_obo_path} ({age_days} days old)")
            return go_obo_path
        else:
            logger.warning(f"GO OBO file is {age_days} days old (>{GO_OBO_MAX_AGE_DAYS}); "
                           f"re-downloading to {go_obo_path}")
    else:
        logger.info(f"GO OBO file not found at {go_obo_path}; downloading")

    download_go_obo(str(obo))
    return str(obo)


#########################
# BLAST RESULTS PARSING
#########################

def parse_blast_results(blast_file: str,
                        evalue_threshold: float = 1e-5) -> tuple:
    """
    Parse BLAST/Diamond output (format 6 with stitle) keyed by transcript ID.

    Collects all hits per transcript, applying only e-value filtering.
    No sorting or truncation is performed here.

    Args:
        blast_file: Path to BLAST output file
        evalue_threshold: E-value cutoff for filtering hits (default: 1e-5)

    Returns:
        tuple of (dict of transcript_id -> list of BlastHit, BlastStats)
    """
    logger.info(f"Parsing BLAST results from {blast_file}")
    logger.info(f"E-value threshold: {evalue_threshold}")

    transcript_hits = defaultdict(list)
    all_queries_in_file = set()
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
                stitle='\t'.join(parts[12:]) if len(parts) > 12 else ""
            )

            transcript_hits[query_id].append(hit)

            # Track database source counts
            if hit.is_swissprot:
                stats.swissprot_hits += 1
            elif hit.is_trembl:
                stats.trembl_hits += 1
            else:
                stats.other_hits += 1

    # Update stats
    stats.queries_in_file = len(all_queries_in_file)
    stats.queries_after_evalue_filter = len(transcript_hits)
    stats.hits_after_evalue_filter = stats.swissprot_hits + stats.trembl_hits + stats.other_hits

    logger.info(f"Parsed hits for {len(transcript_hits):,} transcripts")
    logger.info(f"  Total hits in file: {stats.total_hits_in_file:,}")
    if stats.hits_filtered_by_evalue > 0:
        logger.info(f"  Filtered by e-value (>{evalue_threshold}): {stats.hits_filtered_by_evalue:,}")
    logger.info(f"  Hits passing e-value filter: {stats.hits_after_evalue_filter:,}")
    logger.info(f"  SwissProt (sp|): {stats.swissprot_hits:,}")
    logger.info(f"  TrEMBL (tr|): {stats.trembl_hits:,}")
    if stats.other_hits > 0:
        logger.info(f"  Other: {stats.other_hits:,}")

    return dict(transcript_hits), stats


##########################################
# PER-TRANSCRIPT FILTERING AND SELECTION 
##########################################

def filter_and_select_hits(transcript_hits: dict,
                           bad_name_filter: BadNameFilter,
                           top_n: int = 30) -> dict:
    """
    For each transcript, clean names, filter bad names, apply conditional
    SwissProt preference, and select top-N hits.

    Returns:
        dict of transcript_id -> list of (BlastHit, cleaned_name, valid_symbol) tuples
    """
    logger.info(f"Filtering and selecting top-{top_n} hits per transcript")

    filtered = {}
    total_input = 0
    total_surviving = 0
    total_after_topn = 0

    for transcript_id, hits in transcript_hits.items():
        total_input += len(hits)

        # Step 2a: Clean and filter every hit
        surviving = []
        for hit in hits:
            raw_name = hit.gene_name
            cleaned_name = bad_name_filter.clean_name(raw_name)
            if cleaned_name is None:
                continue
            if bad_name_filter.is_bad_name(cleaned_name):
                continue

            symbol = hit.gene_symbol
            if bad_name_filter.is_bad_symbol(symbol):
                valid_symbol = ""
            else:
                valid_symbol = symbol

            surviving.append((hit, cleaned_name, valid_symbol))

        total_surviving += len(surviving)

        if not surviving:
            continue

        # Step 2b: Sort by bitscore descending
        surviving.sort(key=lambda t: -t[0].bitscore)

        # Step 2c: Top-N truncation
        selected = surviving[:top_n]
        total_after_topn += len(selected)
        filtered[transcript_id] = selected

    logger.info(f"  Transcripts with surviving hits: {len(filtered):,}")
    logger.info(f"  Hits: {total_input:,} input -> "
                f"{total_surviving:,} after filtering -> "
                f"{total_after_topn:,} after top-{top_n}")

    return filtered


###################
# ISOFORM MERGING
###################

def merge_hits_by_gene(filtered_transcript_hits: dict,
                       transcript_to_gene: dict) -> dict:
    """
    Merge filtered hits from all transcripts of each gene into a single pool.

    Each hit is tagged with its source transcript ID.

    Returns:
        dict of gene_id -> list of (BlastHit, cleaned_name, valid_symbol, source_transcript) tuples
    """
    gene_hits = defaultdict(list)

    for transcript_id, hits in filtered_transcript_hits.items():
        gene_id = transcript_to_gene.get(transcript_id, transcript_id)
        for hit, cleaned_name, valid_symbol in hits:
            gene_hits[gene_id].append((hit, cleaned_name, valid_symbol, transcript_id))

    return dict(gene_hits)


############################
# CLUSTERING AND SIMILARITY 
############################

def normalize_gene_name(name: str) -> str:
    """
    Normalize a gene name for comparison purposes.
    Performs only lowercase and whitespace collapsing.
    """
    if not name:
        return ""
    name = name.lower()
    name = ' '.join(name.split())
    return name


def name_similarity(norm1: str, norm2: str) -> float:
    """
    Compute name-only similarity between two normalized protein names.
    Returns max(Jaccard token similarity, SequenceMatcher ratio).
    Both inputs must already be normalized (lowercase, whitespace-collapsed).
    """
    tokens1 = set(re.findall(r'\w+', norm1))
    tokens2 = set(re.findall(r'\w+', norm2))
    if tokens1 and tokens2:
        jaccard = len(tokens1 & tokens2) / len(tokens1 | tokens2)
    else:
        jaccard = 0.0

    seq_ratio = SequenceMatcher(None, norm1, norm2).ratio()

    return max(jaccard, seq_ratio)


def pair_similarity(name1: str, symbol1: str, name2: str, symbol2: str) -> float:
    """
    Compute similarity between two hits using both protein name and gene symbol.

    Uses max of Jaccard token similarity and SequenceMatcher ratio on
    normalized names. When both hits have gene symbols, symbol agreement
    boosts similarity and symbol disagreement caps it.
    """
    norm1 = normalize_gene_name(name1)
    norm2 = normalize_gene_name(name2)

    name_sim = name_similarity(norm1, norm2)

    # Symbol modulation
    if symbol1 and symbol2:
        sym_sim = SequenceMatcher(None, symbol1.lower(), symbol2.lower()).ratio()
        if sym_sim == 1.0:
            return min(1.0, name_sim + 0.2)
        if sym_sim < 0.8:
            return min(name_sim, 0.3)
        return name_sim

    return name_sim


def cluster_hits_agglomerative(hits_with_names: list,
                                similarity_threshold: float = 0.5) -> list:
    """
    Perform average-linkage agglomerative clustering on hits.

    Args:
        hits_with_names: list of (BlastHit, cleaned_name, valid_symbol, source_transcript) tuples
        similarity_threshold: minimum similarity to merge clusters

    Returns:
        list of clusters, where each cluster is a list of indices into hits_with_names
    """
    n = len(hits_with_names)
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    # Compute pairwise similarity matrix
    sim = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            s = pair_similarity(
                hits_with_names[i][1], hits_with_names[i][2],
                hits_with_names[j][1], hits_with_names[j][2]
            )
            sim[i][j] = s
            sim[j][i] = s

    # Initialize each hit as its own cluster
    # Use a dict so we can delete merged clusters efficiently
    clusters = {i: [i] for i in range(n)}

    while len(clusters) > 1:
        # Find the pair of clusters with highest average inter-cluster similarity
        best_sim = -1.0
        best_pair = None
        cluster_ids = list(clusters.keys())

        for idx_a in range(len(cluster_ids)):
            for idx_b in range(idx_a + 1, len(cluster_ids)):
                ca_id = cluster_ids[idx_a]
                cb_id = cluster_ids[idx_b]
                ca = clusters[ca_id]
                cb = clusters[cb_id]

                # Average linkage: mean of all pairwise similarities
                total = sum(sim[i][j] for i in ca for j in cb)
                avg = total / (len(ca) * len(cb))

                if avg > best_sim:
                    best_sim = avg
                    best_pair = (ca_id, cb_id)

        if best_sim < similarity_threshold:
            break

        # Merge the two closest clusters
        a_id, b_id = best_pair
        clusters[a_id] = clusters[a_id] + clusters[b_id]
        del clusters[b_id]

    return list(clusters.values())


################################
# CLUSTER SCORING AND SELECTION 
################################

def select_winning_cluster(clusters: list,
                           hits_with_names: list) -> tuple:
    """
    Score clusters and select the winner. ALWAYS returns the best-scoring
    cluster. There is no fallback to all hits.

    Score = sum(bitscores) * sqrt(cluster_size) * (1 + 0.1 * (n_isoforms - 1))

    Returns:
        (winning_cluster_indices: list, cluster_fraction: float)

    cluster_fraction is the fraction of total hits in the winning cluster.
    This is metadata for logging/output only. It does NOT change which
    indices are returned.
    """
    if not clusters:
        return [], 0.0

    total_hits = len(hits_with_names)
    best_score = -1.0
    best_cluster = None

    for cluster in clusters:
        bitscores_sum = sum(hits_with_names[i][0].bitscore for i in cluster)
        size = len(cluster)
        # Count distinct source transcripts in this cluster
        n_isoforms = len(set(hits_with_names[i][3] for i in cluster))
        score = bitscores_sum * math.sqrt(size) * (1 + 0.1 * (n_isoforms - 1))

        if score > best_score:
            best_score = score
            best_cluster = cluster

    cluster_fraction = len(best_cluster) / total_hits if total_hits > 0 else 0.0

    return best_cluster, cluster_fraction


###############################
# PAIRED NAME/SYMBOL SELECTION 
###############################

def select_name_and_symbol(winning_indices: list,
                           hits_with_names: list,
                           bad_name_filter: BadNameFilter,
                           uniprot_to_geneid: dict = None,
                           geneid_to_symbol: dict = None) -> tuple:
    """
    Select the best protein name and gene symbol from the winning cluster.

    Groups hits by gene symbol. Highest-scoring symbol group determines
    both the gene symbol and the protein name (from that group's best hit).

    Symbol group scoring: mean(bitscores) = sum(bitscores) / count.

    Single-hit exception: a symbol group with only 1 hit can only win
    if no other symbol group has 2+ hits.

    Returns:
        (protein_name: str, gene_symbol: str, winning_symbol_indices: list)

    winning_symbol_indices are the indices (into hits_with_names) of the
    hits in the winning symbol group. These are used downstream for GO
    term collection and for populating filtered_hits.
    """
    if not winning_indices:
        return "", "", []

    # Group winning_indices by valid_symbol (lowercased), tracking original indices
    symbol_groups = defaultdict(list)  # sym_key -> list of indices into hits_with_names
    for idx in winning_indices:
        hit, cleaned_name, valid_symbol, source_transcript = hits_with_names[idx]
        key = valid_symbol.lower() if valid_symbol else ""
        symbol_groups[key].append(idx)

    # Separate no-symbol group
    no_symbol_indices = symbol_groups.pop("", [])

    if symbol_groups:
        # Separate into multi-hit (2+ hits) and single-hit (1 hit) groups
        multi_hit = {s: idxs for s, idxs in symbol_groups.items() if len(idxs) >= 2}
        single_hit = {s: idxs for s, idxs in symbol_groups.items() if len(idxs) == 1}

        # Score only multi-hit groups if any exist; otherwise score single-hit
        candidates = multi_hit if multi_hit else single_hit

        # Score = mean(bitscores) = sum(bitscores) / count
        best_symbol_key = None
        best_score = -1.0
        for sym, idxs in candidates.items():
            score = sum(hits_with_names[i][0].bitscore for i in idxs) / len(idxs)
            if score > best_score:
                best_score = score
                best_symbol_key = sym

        winning_sym_indices = candidates[best_symbol_key]

        # Guard: reject a lone symbol in a large cluster where it
        # represents less than 10% of hits — likely a misannotation.
        if not (len(winning_sym_indices) == 1
                and len(winning_indices) >= 10
                and len(winning_sym_indices) / len(winning_indices) < 0.1):
            # Protein name from highest-bitscore hit in winning symbol group
            best_idx = max(winning_sym_indices, key=lambda i: hits_with_names[i][0].bitscore)
            protein_name = hits_with_names[best_idx][1]   # cleaned_name
            gene_symbol = hits_with_names[best_idx][2]    # original-case valid_symbol

            return protein_name, gene_symbol, winning_sym_indices

    # No-symbol path (also reached when single-symbol guard fires above)
    best_idx = max(winning_indices, key=lambda i: hits_with_names[i][0].bitscore)
    protein_name = hits_with_names[best_idx][1]
    gene_symbol = ""

    # NCBI fallback on best hit's accession
    if uniprot_to_geneid is not None and geneid_to_symbol is not None:
        accessions = extract_accession(hits_with_names[best_idx][0].subject_id)
        for acc in accessions:
            ncbi_symbol = get_ncbi_symbol(acc, uniprot_to_geneid, geneid_to_symbol)
            if ncbi_symbol and not bad_name_filter.is_bad_symbol(ncbi_symbol):
                gene_symbol = ncbi_symbol
                break

    # No symbol groups: entire winning cluster is the symbol group
    return protein_name, gene_symbol, list(winning_indices)


###################
# GO TERM MAPPING
###################

def collect_accessions_from_blast(filtered_hits: dict) -> set:
    """
    Collect all possible accession formats from filtered BLAST results.
    Used to filter GO mapping loading to only relevant entries.

    Args:
        filtered_hits: dict of transcript_id -> list of (BlastHit, cleaned_name, valid_symbol) tuples
    """
    accessions = set()
    for hits in filtered_hits.values():
        for hit_tuple in hits:
            accessions.update(extract_accession(hit_tuple[0].subject_id))
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


################################
# NCBI CROSS-REFERENCE LOADING
################################

def load_ncbi_idmapping(idmapping_file: str, accession_filter: set = None) -> dict:
    """
    Load UniProt accession to NCBI GeneID mapping from idmapping_selected.tsv.

    File format (tab-separated):
      Column 1: UniProt accession (e.g., P31946)
      Column 3: NCBI GeneID (e.g., 7529)

    Args:
        idmapping_file: Path to idmapping_selected.tsv
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


###########################
# SQLITE DATABASE LOADING
###########################

def load_go_mapping_sqlite(db_path: str, accession_filter: set = None) -> dict:
    """
    Load accession to GO term mapping from SQLite database.

    Uses a temporary table + JOIN pattern for efficient filtered queries
    when an accession filter is provided.

    Args:
        db_path: Path to the GOAnnotate SQLite database.
        accession_filter: If provided, only load entries for these accessions.

    Returns: dict of accession -> set of GO IDs
    """
    import sqlite3

    logger.info(f"Loading GO mapping from SQLite: {db_path}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter):,} accessions from BLAST results")

    conn = sqlite3.connect(db_path)
    mapping = defaultdict(set)

    if accession_filter:
        conn.execute("CREATE TEMP TABLE _acc_filter (acc TEXT PRIMARY KEY)")
        batch = [(acc,) for acc in accession_filter]
        for i in range(0, len(batch), 100_000):
            conn.executemany(
                "INSERT OR IGNORE INTO _acc_filter VALUES (?)",
                batch[i:i + 100_000]
            )
        cursor = conn.execute(
            "SELECT g.accession, g.go_terms FROM go_mapping g "
            "INNER JOIN _acc_filter f ON g.accession = f.acc"
        )
    else:
        cursor = conn.execute("SELECT accession, go_terms FROM go_mapping")

    for accession, go_terms_raw in cursor:
        for part in go_terms_raw.split('\t'):
            for go_term in re.split(r'[;,]', part):
                go_term = go_term.strip()
                if go_term.startswith('GO:'):
                    mapping[accession].add(sys.intern(go_term))

    conn.close()
    logger.info(f"Loaded GO mappings for {len(mapping):,} accessions")
    return dict(mapping)


def load_ncbi_idmapping_sqlite(db_path: str, accession_filter: set = None) -> dict:
    """
    Load UniProt accession to NCBI GeneID mapping from SQLite database.

    Uses a temporary table + JOIN pattern for efficient filtered queries
    when an accession filter is provided.

    Args:
        db_path: Path to the GOAnnotate SQLite database.
        accession_filter: If provided, only load entries for these accessions.

    Returns: dict of UniProt accession -> NCBI GeneID (as string)
    """
    import sqlite3

    logger.info(f"Loading UniProt->NCBI GeneID mapping from SQLite: {db_path}")
    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter):,} accessions from BLAST results")

    conn = sqlite3.connect(db_path)

    if accession_filter:
        conn.execute("CREATE TEMP TABLE _acc_filter (acc TEXT PRIMARY KEY)")
        batch = [(acc,) for acc in accession_filter]
        for i in range(0, len(batch), 100_000):
            conn.executemany(
                "INSERT OR IGNORE INTO _acc_filter VALUES (?)",
                batch[i:i + 100_000]
            )
        cursor = conn.execute(
            "SELECT m.uniprot_acc, m.gene_id FROM ncbi_idmapping m "
            "INNER JOIN _acc_filter f ON m.uniprot_acc = f.acc"
        )
    else:
        cursor = conn.execute(
            "SELECT uniprot_acc, gene_id FROM ncbi_idmapping"
        )

    mapping = {}
    for uniprot_acc, gene_id in cursor:
        mapping[uniprot_acc] = gene_id

    conn.close()
    logger.info(f"Loaded {len(mapping):,} UniProt->GeneID mappings")
    return mapping


def load_ncbi_geneinfo_sqlite(db_path: str, geneid_filter: set = None) -> dict:
    """
    Load NCBI GeneID to gene symbol mapping from SQLite database.

    Uses a temporary table + JOIN pattern for efficient filtered queries
    when a GeneID filter is provided.

    Args:
        db_path: Path to the GOAnnotate SQLite database.
        geneid_filter: If provided, only load entries for these GeneIDs.

    Returns: dict of GeneID (as string) -> gene symbol
    """
    import sqlite3

    logger.info(f"Loading NCBI GeneID->Symbol mapping from SQLite: {db_path}")
    if geneid_filter:
        logger.info(f"Filtering to {len(geneid_filter):,} GeneIDs")

    conn = sqlite3.connect(db_path)

    if geneid_filter:
        conn.execute("CREATE TEMP TABLE _gid_filter (gid TEXT PRIMARY KEY)")
        batch = [(gid,) for gid in geneid_filter]
        for i in range(0, len(batch), 100_000):
            conn.executemany(
                "INSERT OR IGNORE INTO _gid_filter VALUES (?)",
                batch[i:i + 100_000]
            )
        cursor = conn.execute(
            "SELECT g.gene_id, g.symbol FROM ncbi_geneinfo g "
            "INNER JOIN _gid_filter f ON g.gene_id = f.gid"
        )
    else:
        cursor = conn.execute(
            "SELECT gene_id, symbol FROM ncbi_geneinfo"
        )

    mapping = {}
    for gene_id, symbol in cursor:
        mapping[gene_id] = symbol

    conn.close()
    logger.info(f"Loaded {len(mapping):,} GeneID->Symbol mappings")
    return mapping


###########################
# MAIN ANNOTATION PIPELINE
###########################

# Module-level globals for worker processes (set by _init_worker)
_worker_go_mapping = None
_worker_go_hierarchy = None
_worker_bad_name_filter = None
_worker_uniprot_to_geneid = None
_worker_geneid_to_symbol = None


def _init_worker(go_mapping, go_hierarchy, bad_name_filter,
                 uniprot_to_geneid, geneid_to_symbol):
    """
    Initialize shared read-only data in each worker process.
    Called once per worker when the pool is created.
    """
    global _worker_go_mapping, _worker_go_hierarchy, _worker_bad_name_filter
    global _worker_uniprot_to_geneid, _worker_geneid_to_symbol
    _worker_go_mapping = go_mapping
    _worker_go_hierarchy = go_hierarchy
    _worker_bad_name_filter = bad_name_filter
    _worker_uniprot_to_geneid = uniprot_to_geneid
    _worker_geneid_to_symbol = geneid_to_symbol


def annotate_single_gene(gene_id: str,
                         hits_with_names: list,
                         consensus_threshold: float,
                         namespaces: list) -> tuple:
    """
    Run the full annotation pipeline for a single gene.

    This function is called by worker processes in the multiprocessing pool.
    It accesses shared read-only data via module-level globals set by
    _init_worker().

    Args:
        gene_id: The gene identifier.
        hits_with_names: List of (BlastHit, cleaned_name, valid_symbol,
                         source_transcript) tuples for this gene.
        consensus_threshold: Similarity threshold for clustering.
        namespaces: GO namespace filter list, or None.

    Returns:
        Tuple of (gene_id, result_dict) where result_dict is a plain dict
        containing all annotation fields and statistics flags.
    """
    try:
        go_mapping = _worker_go_mapping
        go_hierarchy = _worker_go_hierarchy
        bad_name_filter = _worker_bad_name_filter
        uniprot_to_geneid = _worker_uniprot_to_geneid
        geneid_to_symbol = _worker_geneid_to_symbol

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
        }

        if not hits_with_names:
            return (gene_id, result)

        result['had_hits'] = True
        result['n_isoforms'] = len(set(t[3] for t in hits_with_names))

        # Store all merged hits sorted by bitscore
        result['hits'] = sorted(
            [t[0] for t in hits_with_names],
            key=lambda h: h.bitscore, reverse=True
        )

        # Phase 4: Cluster
        clusters = cluster_hits_agglomerative(
            hits_with_names,
            similarity_threshold=consensus_threshold
        )

        # Phase 5: Select winning cluster
        winning_indices, cluster_fraction = select_winning_cluster(
            clusters, hits_with_names
        )
        result['cluster_fraction'] = cluster_fraction
        result['n_cluster_hits'] = len(winning_indices)
        result['winning_cluster_names'] = [
            hits_with_names[i][1] for i in winning_indices
        ]
        result['winning_cluster_symbols'] = [
            hits_with_names[i][2] for i in winning_indices if hits_with_names[i][2]
        ]
        if winning_indices:
            result['mean_cluster_bitscore'] = (
                sum(hits_with_names[i][0].bitscore for i in winning_indices)
                / len(winning_indices)
            )

        # Phase 6: Select name and symbol
        ncbi_available = uniprot_to_geneid is not None and geneid_to_symbol is not None
        protein_name, gene_symbol, winning_symbol_indices = select_name_and_symbol(
            winning_indices, hits_with_names,
            bad_name_filter,
            uniprot_to_geneid=uniprot_to_geneid if ncbi_available else None,
            geneid_to_symbol=geneid_to_symbol if ncbi_available else None
        )

        result['consensus_name'] = protein_name
        result['consensus_symbol'] = gene_symbol.lower() if gene_symbol else ""

        if protein_name:
            result['had_protein_name'] = True
        if gene_symbol:
            result['had_symbol'] = True
            # Determine symbol source
            if any(hits_with_names[i][2] for i in winning_indices):
                result['symbol_source'] = 'uniprot'
            else:
                result['symbol_source'] = 'ncbi'

        # Concordance: compare final name/symbol against full merged hit pool
        if protein_name and hits_with_names:
            final_name_norm = normalize_gene_name(protein_name)
            name_match_count = sum(
                1 for _h, cn, _vs, _st in hits_with_names
                if name_similarity(final_name_norm, normalize_gene_name(cn)) >= consensus_threshold
            )
            result['name_concordance'] = name_match_count / len(hits_with_names)

        if gene_symbol and hits_with_names:
            final_sym_lower = gene_symbol.lower()
            symbol_match_count = sum(
                1 for _h, _cn, vs, _st in hits_with_names
                if vs and vs.lower() == final_sym_lower
            )
            result['symbol_concordance'] = symbol_match_count / len(hits_with_names)

        # Set filtered_hits from winning symbol group (sorted by bitscore)
        result['filtered_hits'] = sorted(
            [hits_with_names[i][0] for i in winning_symbol_indices],
            key=lambda h: h.bitscore, reverse=True
        )

        # Phase 7: GO terms from winning symbol group (top 5 by bitscore)
        go_source_indices = sorted(
            winning_symbol_indices,
            key=lambda i: hits_with_names[i][0].bitscore, reverse=True
        )[:5]
        go_hits = [hits_with_names[i][0] for i in go_source_indices]
        go_terms = get_go_terms_for_hits(go_hits, go_mapping)
        go_source = "cluster_top5" if go_terms else ""

        # Phase 7b: SwissProt GO supplement
        winning_name_norm = normalize_gene_name(protein_name) if protein_name else ""
        winning_sym_lower = gene_symbol.lower() if gene_symbol else ""
        go_source_set = set(go_source_indices)

        for idx in winning_indices:
            if idx in go_source_set:
                continue
            hit, cleaned_name, symbol, _src = hits_with_names[idx]
            if not hit.is_swissprot:
                continue
            sym_match = (winning_sym_lower and symbol
                         and symbol.lower() == winning_sym_lower)
            name_match = False
            if not sym_match and winning_name_norm and cleaned_name:
                name_norm = normalize_gene_name(cleaned_name)
                sim = pair_similarity(winning_name_norm, "", name_norm, "")
                name_match = sim >= consensus_threshold
            if sym_match or name_match:
                sp_terms = get_go_terms_for_hits([hit], go_mapping)
                if sp_terms:
                    new_terms = sp_terms - go_terms if go_terms else sp_terms
                    if new_terms:
                        go_source = "cluster_top5+swissprot_supplement"
                    go_terms = go_terms | sp_terms if go_terms else sp_terms

        result['go_source'] = go_source

        # Filter by namespace if specified
        if namespaces and go_terms:
            go_terms = go_hierarchy.filter_by_namespace(go_terms, namespaces)

        result['go_terms'] = go_terms

        if go_terms:
            result['specific_go_terms'] = go_hierarchy.filter_to_specific(go_terms)

        return (gene_id, result)

    except Exception as e:
        raise RuntimeError(f"Error annotating gene {gene_id}: {e}") from e


def annotate_queries(filtered_transcript_hits: dict,
                     transcript_to_gene: dict,
                     gene_to_transcripts: dict,
                     go_mapping: dict,
                     go_hierarchy: GOHierarchy,
                     bad_name_filter: BadNameFilter,
                     all_gene_ids: list,
                     consensus_threshold: float = 0.5,
                     namespaces: list = None,
                     uniprot_to_geneid: dict = None,
                     geneid_to_symbol: dict = None,
                     threads: int = None) -> dict:
    """
    Main annotation pipeline: merge isoform hits, cluster, select name/symbol,
    collect GO terms. Per-gene annotation runs in parallel via multiprocessing.

    For each gene:
      1. Merge hits across isoforms (Phase 3)
      2. Cluster hits by name/symbol similarity (Phase 4)
      3. Score and select winning cluster (Phase 5)
      4. Select protein name and gene symbol (Phase 6)
      5. Set filtered_hits from winning symbol group (Phase 6 output)
      6. Collect GO terms from winning symbol group (Phase 7)

    After all genes: symbol propagation (Phase 8), concordance recalculation,
    and statistics logging (sequential).

    Returns: dict of gene_id -> QueryAnnotation
    """
    # Phase 3: Merge hits by gene
    gene_hits = merge_hits_by_gene(filtered_transcript_hits, transcript_to_gene)

    # Determine worker count
    n_workers = threads if threads is not None else os.cpu_count()
    n_genes_with_hits = sum(1 for gid in all_gene_ids if gid in gene_hits)
    n_workers = min(n_workers, max(1, n_genes_with_hits))

    logger.info(f"Annotating {len(all_gene_ids):,} genes "
                f"({n_genes_with_hits:,} have filtered hits) "
                f"using {n_workers} workers...")

    # Build task list
    tasks = []
    for gene_id in all_gene_ids:
        hits = gene_hits.get(gene_id, [])
        tasks.append((gene_id, hits, consensus_threshold, namespaces))

    # Run per-gene annotation
    start_time = time.time()

    if n_workers == 1:
        # Single-process mode: avoid multiprocessing overhead
        _init_worker(go_mapping, go_hierarchy, bad_name_filter,
                     uniprot_to_geneid, geneid_to_symbol)
        results = [annotate_single_gene(*task) for task in tasks]
    else:
        chunksize = max(1, len(tasks) // (n_workers * 4))
        with multiprocessing.Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=(go_mapping, go_hierarchy, bad_name_filter,
                      uniprot_to_geneid, geneid_to_symbol)
        ) as pool:
            results = pool.starmap(annotate_single_gene, tasks,
                                   chunksize=chunksize)

    elapsed = time.time() - start_time
    genes_per_sec = len(all_gene_ids) / elapsed if elapsed > 0 else 0
    logger.info(f"  Per-gene annotation completed in {elapsed:.1f}s "
                f"({genes_per_sec:.0f} genes/sec)")

    # Collect results into annotations dict and accumulate statistics
    annotations = {}
    cluster_fractions = []
    stats = {
        'total_genes': len(all_gene_ids),
        'genes_with_hits': 0,
        'genes_with_protein_name': 0,
        'genes_with_go_terms': 0,
        'symbols_from_uniprot': 0,
        'symbols_from_ncbi': 0,
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
        annotations[gene_id] = annot

        if result_dict['had_hits']:
            stats['genes_with_hits'] += 1
        if result_dict['had_protein_name']:
            stats['genes_with_protein_name'] += 1
        if result_dict['had_symbol']:
            if result_dict['symbol_source'] == 'uniprot':
                stats['symbols_from_uniprot'] += 1
            elif result_dict['symbol_source'] == 'ncbi':
                stats['symbols_from_ncbi'] += 1
        if result_dict['go_terms']:
            stats['genes_with_go_terms'] += 1
        if result_dict['cluster_fraction'] > 0:
            cluster_fractions.append(result_dict['cluster_fraction'])

    # Phase 8: Symbol propagation
    # For genes sharing an identical normalized protein name, propagate
    # the most common symbol to genes that lack one.
    name_to_genes = defaultdict(list)
    for gene_id, annot in annotations.items():
        if annot.consensus_name:
            norm = normalize_gene_name(annot.consensus_name)
            name_to_genes[norm].append(gene_id)

    symbols_propagated = 0
    propagated_gene_ids = []
    for norm_name, gene_id_list in name_to_genes.items():
        if len(gene_id_list) < 2:
            continue
        existing_symbols = [annotations[gid].consensus_symbol
                            for gid in gene_id_list
                            if annotations[gid].consensus_symbol]
        if not existing_symbols:
            continue
        genes_without = [gid for gid in gene_id_list
                         if not annotations[gid].consensus_symbol]
        if not genes_without:
            continue
        best_symbol = Counter(existing_symbols).most_common(1)[0][0]
        for gid in genes_without:
            annotations[gid].consensus_symbol = best_symbol
            symbols_propagated += 1
            propagated_gene_ids.append(gid)

    stats['symbols_from_propagation'] = symbols_propagated
    if symbols_propagated > 0:
        logger.info(f"  Symbol propagation: {symbols_propagated:,} genes received symbols "
                    f"from genes with identical protein names")

    # Recalculate symbol_concordance for genes that received propagated symbols
    for gid in propagated_gene_ids:
        annot = annotations[gid]
        hits_pool = gene_hits.get(gid, [])
        if annot.consensus_symbol and hits_pool:
            final_sym_lower = annot.consensus_symbol.lower()
            match_count = sum(1 for _h, _cn, vs, _st in hits_pool
                              if vs and vs.lower() == final_sym_lower)
            annot.symbol_concordance = match_count / len(hits_pool)

    # Log statistics
    total = max(1, stats['total_genes'])
    logger.info("Annotation statistics:")
    logger.info(f"  Total genes: {stats['total_genes']:,}")
    logger.info(f"  Genes with filtered hits: {stats['genes_with_hits']:,} "
                f"({100*stats['genes_with_hits']/total:.1f}%)")
    logger.info(f"  Genes with protein name: {stats['genes_with_protein_name']:,} "
                f"({100*stats['genes_with_protein_name']/total:.1f}%)")
    total_symbols = (stats['symbols_from_uniprot'] + stats['symbols_from_ncbi']
                     + stats['symbols_from_propagation'])
    logger.info(f"  Genes with gene symbol: {total_symbols:,} "
                f"({100*total_symbols/total:.1f}%)")
    if stats['symbols_from_uniprot'] > 0:
        logger.info(f"    - From UniProt GN= field: {stats['symbols_from_uniprot']:,}")
    if stats['symbols_from_ncbi'] > 0:
        logger.info(f"    - From NCBI cross-reference: {stats['symbols_from_ncbi']:,}")
    if stats['symbols_from_propagation'] > 0:
        logger.info(f"    - From name-based propagation: {stats['symbols_from_propagation']:,}")
    logger.info(f"  Genes with GO terms: {stats['genes_with_go_terms']:,} "
                f"({100*stats['genes_with_go_terms']/total:.1f}%)")

    # Log cluster fraction statistics
    if cluster_fractions:
        mean_frac = sum(cluster_fractions) / len(cluster_fractions)
        sorted_fracs = sorted(cluster_fractions)
        n = len(sorted_fracs)
        if n % 2 == 0:
            median_frac = (sorted_fracs[n // 2 - 1] + sorted_fracs[n // 2]) / 2
        else:
            median_frac = sorted_fracs[n // 2]
        logger.info(f"  Winning cluster fraction: mean={mean_frac:.2f}, "
                    f"median={median_frac:.2f}, "
                    f"min={sorted_fracs[0]:.2f}, max={sorted_fracs[-1]:.2f}")

    # Log concordance statistics
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


####################
# OUTPUT FUNCTIONS
####################

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


def write_annotated_fasta(input_fasta: str, annotations: dict, output_file: Path,
                          transcript_to_gene: dict = None):
    """
    Write an annotated FASTA file with Name= and product= fields added to headers.
    Replaces any existing Name= or product= fields in the header.

    Args:
        input_fasta: Path to the original FASTA file
        annotations: dict of gene_id -> QueryAnnotation
        output_file: Path to write the annotated FASTA
        transcript_to_gene: Optional mapping from transcript ID to gene ID.
                            If provided, looks up annotation by gene ID.
    """
    logger.info(f"Writing annotated FASTA to {output_file}")

    annotated_count = 0
    total_count = 0

    with open(input_fasta) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                total_count += 1
                header = line[1:].rstrip()
                # Extract first field from header
                first_field = header.split()[0] if header else ""

                # Resolve gene ID via transcript-to-gene mapping if available
                if transcript_to_gene:
                    gene_id = transcript_to_gene.get(first_field, first_field)
                else:
                    gene_id = first_field

                # Look up annotation
                annot = annotations.get(gene_id)

                if annot and (annot.consensus_symbol or annot.consensus_name):
                    annotated_count += 1

                    # Remove any existing Name= or product= fields from header
                    header_cleaned = re.sub(r'\bName=\S+\s*', '', header)
                    header_cleaned = re.sub(r'\bproduct=.*?(?=\s+\w+=|$)', '', header_cleaned)
                    header_cleaned = re.sub(r'  +', ' ', header_cleaned).strip()

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


def _gff3_encode_value(value: str) -> str:
    """
    Percent-encode a GFF3 attribute value per the GFF3 specification.

    Characters requiring encoding: tab, newline, carriage return, %, ;, =, &, comma.
    Spaces are encoded as %20.
    """
    return urllib.parse.quote(value, safe='abcdefghijklmnopqrstuvwxyz'
                              'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
                              '-._~:@!$\'()*+/')


def write_annotated_gff3(input_gff: str, annotations: dict, output_file: Path):
    """
    Write an annotated GFF3 file with cleaned attributes and functional
    annotations following GFF3 conventions.

    Original feature IDs and Parent references are preserved. Custom
    attributes (gene_id, transcript_id) are removed from the output.

    Attributes per feature type:
        gene:  ID, Name, product
        mRNA:  ID, Parent, Name, product
        exon, CDS, five_prime_UTR, three_prime_UTR:  ID, Parent
        start_codon, stop_codon, intron: Parent only

    Args:
        input_gff: Path to the original GFF3 file
        annotations: dict of query_id -> QueryAnnotation
        output_file: Path to write the annotated GFF3
    """
    logger.info(f"Writing annotated GFF3 to {output_file}")

    # Sub-features that keep their ID
    SUB_ID_TYPES = {'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}
    # Sub-features that only receive Parent (no ID)
    PARENT_ONLY_TYPES = {'start_codon', 'stop_codon', 'intron'}

    annotated_genes = 0
    gene_count = 0

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

            feature_type = parts[2]
            attributes = parts[8]

            # Parse attributes into dict for lookup
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr and '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val

            old_id = attr_dict.get('ID')
            old_parent = attr_dict.get('Parent')

            # Determine gene_id for annotation lookup (feature-type-aware)
            if feature_type == 'gene':
                lookup_id = attr_dict.get('gene_id') or old_id
            elif feature_type == 'mRNA':
                # mRNA Parent= points to the gene ID
                lookup_id = old_parent or attr_dict.get('gene_id')
            else:
                lookup_id = attr_dict.get('gene_id')

            annot = annotations.get(lookup_id) if lookup_id else None

            # Build new column-9 attributes based on feature type
            new_attrs = []

            if feature_type == 'gene':
                gene_count += 1
                new_attrs.append(f"ID={old_id}")

                if annot:
                    if annot.consensus_symbol:
                        new_attrs.append(f"Name={annot.consensus_symbol}")
                    if annot.consensus_name:
                        product = strip_uniprot_prefix(annot.consensus_name)
                        new_attrs.append(f"product={product}")
                    if annot.consensus_symbol or annot.consensus_name:
                        annotated_genes += 1

            elif feature_type == 'mRNA':
                new_attrs.append(f"ID={old_id}")
                if old_parent:
                    new_attrs.append(f"Parent={old_parent}")

                if annot:
                    if annot.consensus_symbol:
                        new_attrs.append(f"Name={annot.consensus_symbol}")
                    if annot.consensus_name:
                        product = strip_uniprot_prefix(annot.consensus_name)
                        new_attrs.append(f"product={product}")

            elif feature_type in SUB_ID_TYPES:
                if old_id:
                    new_attrs.append(f"ID={old_id}")
                if old_parent:
                    new_attrs.append(f"Parent={old_parent}")

            elif feature_type in PARENT_ONLY_TYPES:
                if old_parent:
                    new_attrs.append(f"Parent={old_parent}")

            else:
                # Unknown feature type: keep original attributes,
                # remove gene_id/transcript_id
                if old_id:
                    new_attrs.append(f"ID={old_id}")
                if old_parent:
                    new_attrs.append(f"Parent={old_parent}")
                for key, val in attr_dict.items():
                    if key not in ('ID', 'Parent', 'gene_id', 'transcript_id'):
                        new_attrs.append(f"{key}={val}")

            parts[8] = ';'.join(new_attrs) if new_attrs else '.'
            fout.write('\t'.join(parts) + '\n')

    logger.info(f"Annotated {annotated_genes:,} of {gene_count:,} genes in GFF3")


def write_annotation_results(annotations: dict,
                             output_dir: Path,
                             prefix: str,
                             input_fasta: str = None,
                             input_gff: str = None,
                             args: argparse.Namespace = None,
                             blast_stats: BlastStats = None,
                             total_sequences: int = None,
                             transcript_to_gene: dict = None):
    """
    Write annotation results to various output files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Gene-to-symbol mapping
    symbol_file = output_dir / f"{prefix}_gene2symbol.tsv"
    logger.info(f"Writing gene2symbol to {symbol_file}")

    with open(symbol_file, 'w') as f:
        f.write("gene_id\tgene_symbol\n")
        for query_id, annot in sorted(annotations.items()):
            f.write(f"{query_id}\t{annot.consensus_symbol}\n")

    # 2. Gene-to-name mapping
    name_file = output_dir / f"{prefix}_gene2name.tsv"
    logger.info(f"Writing gene2name to {name_file}")

    with open(name_file, 'w') as f:
        f.write("gene_id\tprotein_name\n")
        for query_id, annot in sorted(annotations.items()):
            pname = strip_uniprot_prefix(annot.consensus_name) if annot.consensus_name else ""
            f.write(f"{query_id}\t{pname}\n")

    # 3. Gene-to-GO mapping
    go_file = output_dir / f"{prefix}_gene2go.tsv"
    logger.info(f"Writing gene2go to {go_file}")

    with open(go_file, 'w') as f:
        f.write("gene_id\tgo_terms\n")
        for query_id, annot in sorted(annotations.items()):
            go_str = ';'.join(sorted(annot.specific_go_terms)) if annot.specific_go_terms else ""
            f.write(f"{query_id}\t{go_str}\n")

    # 4. Annotation evidence
    evidence_file = output_dir / f"{prefix}_annotation_evidence.tsv"
    logger.info(f"Writing annotation evidence to {evidence_file}")

    with open(evidence_file, 'w') as f:
        f.write("gene_id\tgene_symbol\tprotein_name\tn_isoforms\tn_total_hits\t"
                "n_cluster_hits\twinning_cluster_names\twinning_cluster_symbols\t"
                "top_hit_sseqid\ttop_hit_bitscore\ttop_hit_source\t"
                "mean_cluster_bitscore\tname_concordance\tsymbol_concordance\t"
                "go_source\tn_go_terms\n")
        for query_id, annot in sorted(annotations.items()):
            pname = strip_uniprot_prefix(annot.consensus_name) if annot.consensus_name else ""
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

            f.write(f"{query_id}\t{annot.consensus_symbol}\t{pname}\t"
                    f"{annot.n_isoforms}\t{len(annot.hits)}\t"
                    f"{annot.n_cluster_hits}\t{cluster_names}\t{cluster_symbols}\t"
                    f"{top_sseqid}\t{top_bitscore}\t{top_source}\t"
                    f"{annot.mean_cluster_bitscore:.2f}\t"
                    f"{annot.name_concordance:.2f}\t{annot.symbol_concordance:.2f}\t"
                    f"{annot.go_source}\t{len(annot.specific_go_terms)}\n")

    # 6. Annotated FASTA file
    if input_fasta:
        fasta_out = output_dir / f"{prefix}_annotated.fasta"
        write_annotated_fasta(input_fasta, annotations, fasta_out,
                              transcript_to_gene=transcript_to_gene)

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
            if args.gff:
                cmd_parts.extend(["--gff", str(args.gff)])
            if args.db:
                cmd_parts.extend(["--db", str(args.db)])
            if args.go_mapping:
                cmd_parts.extend(["--go-mapping", str(args.go_mapping)])
            cmd_parts.extend(["--go-obo", str(args.go_obo)])
            cmd_parts.extend(["--bad-names", str(args.bad_names)])
            cmd_parts.extend(["--output", str(args.output)])
            cmd_parts.extend(["--prefix", str(args.prefix)])
            cmd_parts.extend(["--evalue", str(args.evalue)])
            cmd_parts.extend(["--top-n", str(args.top_n)])
            cmd_parts.extend(["--consensus-threshold", str(args.consensus_threshold)])
            cmd_parts.extend(["--namespace"] + args.namespace)
            f.write(" \\\n    ".join(cmd_parts) + "\n\n")

        # Input files section
        if args:
            f.write("Input Files:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  Transcripts:    {args.transcripts}\n")
            f.write(f"  BLAST results:  {args.blast_results}\n")
            if args.db:
                f.write(f"  SQLite DB:      {args.db}\n")
            if args.go_mapping:
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
            f.write(f"  Mode:                     Annotation only\n")
            f.write(f"  E-value threshold:        {args.evalue}\n")
            f.write(f"  Top N hits:               {args.top_n}\n")
            f.write(f"  Consensus threshold:      {args.consensus_threshold}\n")
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
                f"({100*(1-len(all_specific)/len(all_go)) if all_go else 0:.1f}% reduction)\n\n")

        # Concordance statistics
        name_concs = [a.name_concordance for a in annotations.values() if a.consensus_name]
        sym_concs = [a.symbol_concordance for a in annotations.values() if a.consensus_symbol]

        if name_concs or sym_concs:
            f.write("Concordance Statistics:\n")
            f.write("-" * 70 + "\n")
            if name_concs:
                f.write(f"  Name concordance:   mean={statistics.mean(name_concs):.2f}, "
                        f"median={statistics.median(name_concs):.2f}, "
                        f"min={min(name_concs):.2f}, max={max(name_concs):.2f}\n")
            if sym_concs:
                f.write(f"  Symbol concordance: mean={statistics.mean(sym_concs):.2f}, "
                        f"median={statistics.median(sym_concs):.2f}, "
                        f"min={min(sym_concs):.2f}, max={max(sym_concs):.2f}\n")

    logger.info(f"Results written to {output_dir}")


###################
# CLI
###################

class VerticalHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Help formatter that displays each argument's help on a single line."""

    def __init__(self, prog, indent_increment=2, max_help_position=36, width=None):
        if width is None:
            width = max(shutil.get_terminal_size().columns, 100)
        super().__init__(prog, indent_increment, max_help_position, width)

    def _split_lines(self, text, width):
        return text.splitlines()


def build_parser():
    parser = argparse.ArgumentParser(
        description="GOAnnotate: Gene Ontology annotation pipeline using UniProt BLAST/Diamond results",
        formatter_class=VerticalHelpFormatter,
        epilog="""\
Required: --transcripts, --blast-results, --db (or --go-mapping), --bad-names, --output

  The --db SQLite database (built by build_databases.py) is the recommended way
  to provide GO mappings and NCBI cross-references. Alternatively, use --go-mapping
  with optional --ncbi-idmapping + --ncbi-geneinfo flat files.

  The GO OBO file (--go-obo) is auto-downloaded if not provided or >30 days old.

Examples:
  # Using SQLite database (recommended)
  GOAnnotate.py --transcripts CDS.fasta --blast-results UniProt_results.tsv \\
      --db GOAnnotate_db.sqlite --bad-names bad_names.txt -o output_dir

  # Using flat files
  GOAnnotate.py --transcripts CDS.fasta --blast-results UniProt_results.tsv \\
      --go-mapping GO_mapping.tsv --ncbi-idmapping idmapping_selected.tab \\
      --ncbi-geneinfo gene_info.tsv --bad-names bad_names.txt -o output_dir
"""
    )

    # Input options
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument(
        "--transcripts", metavar="FILE", required=True,
        help="CDS FASTA file with transcript-level headers (gene= field for gene mapping) [required]"
    )
    input_group.add_argument(
        "--blast-results", "--blast", metavar="FILE", required=True,
        help="BLAST/Diamond results file, format 6 with stitle [required]"
    )
    input_group.add_argument(
        "--gff", metavar="FILE",
        help="GFF3 annotation file [optional: for annotated GFF3 output]"
    )

    # GO mapping
    go_group = parser.add_argument_group("GO options")
    go_group.add_argument(
        "--go-mapping", metavar="FILE",
        help="Accession-to-GO mapping file (accession<TAB>GO:xxxx;GO:yyyy) [required unless --db]"
    )
    go_group.add_argument(
        "--db", metavar="FILE",
        help="SQLite database built by build_databases.py "
             "(replaces --go-mapping, --ncbi-idmapping, --ncbi-geneinfo)"
    )
    go_group.add_argument(
        "--go-obo", metavar="FILE",
        help="GO OBO hierarchy file (auto-downloaded if missing or >30 days old)"
    )
    go_group.add_argument(
        "--namespace", nargs='+', default=["BP", "MF", "CC"],
        choices=["BP", "MF", "CC"],
        help="GO namespaces to include (default: BP MF CC)"
    )

    # NCBI cross-reference options
    ncbi_group = parser.add_argument_group("NCBI cross-reference options (for improved gene symbol coverage)")
    ncbi_group.add_argument(
        "--ncbi-idmapping", metavar="FILE",
        help="UniProt ID mapping file (idmapping_selected.tsv) for UniProt->NCBI GeneID lookup"
    )
    ncbi_group.add_argument(
        "--ncbi-geneinfo", metavar="FILE",
        help="NCBI gene_info file for GeneID->Symbol lookup"
    )

    # Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold for filtering hits (default: %(default)s)"
    )
    filter_group.add_argument(
        "--top-n", type=int, default=30,
        help="Number of top BLAST hits to consider per transcript (default: %(default)s)"
    )
    filter_group.add_argument(
        "--bad-names", metavar="FILE",
        help="Bad names pattern file for filtering uninformative protein names [required]"
    )
    filter_group.add_argument(
        "--consensus-threshold", type=float, default=0.5,
        help="Similarity threshold for name clustering, 0-1 (default: %(default)s)"
    )

    # Performance options
    perf_group = parser.add_argument_group("Performance options")
    perf_group.add_argument(
        "--threads", type=int, default=None,
        help="Number of parallel workers for annotation "
             "(default: number of CPUs)"
    )

    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output", "-o", required=True, metavar="DIR",
        help="Output directory (required). Directory name is used as the "
             "file prefix unless overridden with --prefix"
    )
    output_group.add_argument(
        "--prefix", default=None,
        help="Output file prefix (default: derived from output directory name)"
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Initialize components
    logger.info("Initializing GOAnnotate pipeline...")

    ###############################
    # Validate required arguments
    ###############################

    missing = []
    if not args.go_mapping and not args.db:
        missing.append("--go-mapping (or --db)")
    if not args.bad_names:
        missing.append("--bad-names")
    if missing:
        logger.error(f"Missing required arguments: {', '.join(missing)}")
        sys.exit(1)

    # Validate file existence
    if not Path(args.transcripts).exists():
        logger.error(f"Transcripts FASTA file not found: {args.transcripts}")
        sys.exit(1)
    if not Path(args.blast_results).exists():
        logger.error(f"BLAST results file not found: {args.blast_results}")
        sys.exit(1)
    if args.go_mapping and not Path(args.go_mapping).exists():
        logger.error(f"GO mapping file not found: {args.go_mapping}")
        sys.exit(1)
    if args.db and not Path(args.db).exists():
        logger.error(f"SQLite database not found: {args.db}")
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

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Derive prefix from output directory name if not explicitly provided
    if args.prefix is None:
        args.prefix = output_dir.resolve().name

    #######################
    # Resolve GO OBO file
    #######################

    args.go_obo = resolve_go_obo(args.go_obo)

    ############################################################################
    # Phase 0: Parse FASTA headers to build gene universe + transcript-gene map
    ############################################################################

    total_sequences, all_gene_ids, transcript_to_gene, gene_to_transcripts = \
        parse_fasta_headers(args.transcripts)

    #####################################################
    # Phase 1: Parse BLAST results (keyed by transcript)
    #####################################################

    transcript_hits, blast_stats = parse_blast_results(
        args.blast_results,
        evalue_threshold=args.evalue
    )

    ###################################################################
    # Phase 2: Per-transcript cleaning, filtering, and top-N selection
    ###################################################################

    # Initialize bad name filter
    bad_name_filter = BadNameFilter(args.bad_names)

    filtered_transcript_hits = filter_and_select_hits(
        transcript_hits,
        bad_name_filter=bad_name_filter,
        top_n=args.top_n
    )

    ###############################################################
    # Collect accessions from filtered hits for GO mapping loading
    ###############################################################

    blast_accessions = collect_accessions_from_blast(filtered_transcript_hits)

    # Load GO mapping
    if args.db:
        go_mapping = load_go_mapping_sqlite(args.db,
                                            accession_filter=blast_accessions)
    else:
        go_mapping = load_go_mapping(args.go_mapping,
                                     accession_filter=blast_accessions)

    # Load GO hierarchy
    go_hierarchy = GOHierarchy(args.go_obo)

    # Load NCBI cross-reference mappings
    uniprot_to_geneid = None
    geneid_to_symbol = None
    if args.db:
        uniprot_to_geneid = load_ncbi_idmapping_sqlite(
            args.db, accession_filter=blast_accessions
        )
        found_geneids = set(uniprot_to_geneid.values())
        geneid_to_symbol = load_ncbi_geneinfo_sqlite(
            args.db, geneid_filter=found_geneids
        )
    elif args.ncbi_idmapping and args.ncbi_geneinfo:
        uniprot_to_geneid = load_ncbi_idmapping(
            args.ncbi_idmapping,
            accession_filter=blast_accessions
        )
        found_geneids = set(uniprot_to_geneid.values())
        geneid_to_symbol = load_ncbi_geneinfo(
            args.ncbi_geneinfo,
            geneid_filter=found_geneids
        )

    #####################################################
    # Phases 3-7: Annotate (merge, cluster, select, GO)
    #####################################################

    annotations = annotate_queries(
        filtered_transcript_hits=filtered_transcript_hits,
        transcript_to_gene=transcript_to_gene,
        gene_to_transcripts=gene_to_transcripts,
        go_mapping=go_mapping,
        go_hierarchy=go_hierarchy,
        bad_name_filter=bad_name_filter,
        all_gene_ids=all_gene_ids,
        consensus_threshold=args.consensus_threshold,
        namespaces=args.namespace,
        uniprot_to_geneid=uniprot_to_geneid,
        geneid_to_symbol=geneid_to_symbol,
        threads=args.threads
    )

    ################
    # Write outputs
    ################

    write_annotation_results(
        annotations=annotations,
        output_dir=output_dir,
        prefix=args.prefix,
        input_fasta=args.transcripts,
        input_gff=args.gff,
        args=args,
        blast_stats=blast_stats,
        total_sequences=total_sequences,
        transcript_to_gene=transcript_to_gene
    )

    logger.info("GOAnnotate pipeline completed successfully")


if __name__ == "__main__":
    main()
