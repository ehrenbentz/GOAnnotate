#!/usr/bin/env python3
"""
create_go_mapping.py - Generate GO term mappings for GOAnnotate

This script creates accession-to-GO mapping files from various sources:

1. UniProt ID Mapping (idmapping_selected.tab.gz)
   - Best for SwissProt/TrEMBL BLAST databases
   - Download from: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

2. NCBI gene2go
   - Best for RefSeq/nr BLAST databases
   - Download from: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

3. UniProt REST API
   - Query GO terms for accessions extracted from BLAST results
   - Good for smaller datasets or when you don't want to download large files

4. Extract from BLAST hits
   - Parse BLAST results and extract accessions, then query UniProt

Usage examples:

  # From UniProt idmapping file (recommended for UniProt databases)
  python create_go_mapping.py uniprot-idmapping \\
      --idmapping idmapping_selected.tab.gz \\
      --output uniprot2go.tsv

  # From NCBI gene2go (recommended for RefSeq/nr)
  python create_go_mapping.py ncbi-gene2go \\
      --gene2go gene2go.gz \\
      --gene2accession gene2accession.gz \\
      --output ncbi2go.tsv

  # From UniProt API using BLAST results
  python create_go_mapping.py uniprot-api \\
      --blast-results my_blast.tsv \\
      --output blast_hits2go.tsv

  # Download files automatically
  python create_go_mapping.py download --source uniprot --output-dir ./db_files
  python create_go_mapping.py download --source ncbi --output-dir ./db_files

"""

import os
import sys
import re
import gzip
import argparse
import logging
import time
import json
from pathlib import Path
from collections import defaultdict
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
import ssl

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================

DOWNLOAD_URLS = {
    'uniprot_idmapping': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz',
    'uniprot_idmapping_sprot': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/',
    'ncbi_gene2go': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz',
    'ncbi_gene2accession': 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz',
    'go_obo': 'http://purl.obolibrary.org/obo/go.obo',
}


def download_file(url: str, output_path: str, chunk_size: int = 8192):
    """Download a file with progress indication."""
    logger.info(f"Downloading {url}")
    logger.info(f"  -> {output_path}")

    # Create SSL context that doesn't verify (for some FTP servers)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    try:
        req = Request(url, headers={'User-Agent': 'Python/GOAnnotate'})
        with urlopen(req, context=ctx) as response:
            total_size = response.headers.get('Content-Length')
            if total_size:
                total_size = int(total_size)
                logger.info(f"  File size: {total_size / 1e9:.2f} GB")

            downloaded = 0
            last_report = 0

            with open(output_path, 'wb') as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)

                    # Report progress every 100MB
                    if total_size and downloaded - last_report > 100_000_000:
                        pct = 100 * downloaded / total_size
                        logger.info(f"  Progress: {pct:.1f}% ({downloaded / 1e9:.2f} GB)")
                        last_report = downloaded

        logger.info(f"  Download complete: {output_path}")
        return output_path

    except (HTTPError, URLError) as e:
        logger.error(f"Download failed: {e}")
        raise


def download_sources(source: str, output_dir: str):
    """Download required source files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if source == 'uniprot':
        # Download UniProt idmapping (WARNING: ~25GB compressed!)
        logger.warning("UniProt idmapping_selected.tab.gz is ~25GB compressed!")
        logger.warning("Consider using organism-specific files instead.")
        logger.warning("See: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/")

        response = input("Download full idmapping file? (y/n): ")
        if response.lower() == 'y':
            download_file(
                DOWNLOAD_URLS['uniprot_idmapping'],
                str(output_dir / 'idmapping_selected.tab.gz')
            )

    elif source == 'ncbi':
        # Download NCBI gene2go and gene2accession
        download_file(
            DOWNLOAD_URLS['ncbi_gene2go'],
            str(output_dir / 'gene2go.gz')
        )
        logger.info("gene2go downloaded. For RefSeq mapping, also downloading gene2accession...")
        logger.warning("gene2accession.gz is ~2GB compressed")

        response = input("Download gene2accession? (y/n): ")
        if response.lower() == 'y':
            download_file(
                DOWNLOAD_URLS['ncbi_gene2accession'],
                str(output_dir / 'gene2accession.gz')
            )

    elif source == 'go':
        download_file(
            DOWNLOAD_URLS['go_obo'],
            str(output_dir / 'go.obo')
        )

    else:
        logger.error(f"Unknown source: {source}")
        sys.exit(1)


# =============================================================================
# UNIPROT ID MAPPING PARSER
# =============================================================================

def parse_uniprot_idmapping(idmapping_file: str, output_file: str,
                            accession_filter: set = None,
                            taxon_filter: set = None):
    """
    Parse UniProt idmapping_selected.tab.gz file.

    The file has these columns (tab-separated):
    1. UniProtKB-AC
    2. UniProtKB-ID
    3. GeneID (EntrezGene)
    4. RefSeq
    5. GI
    6. PDB
    7. GO
    8. UniRef100
    9. UniRef90
    10. UniRef50
    11. UniParc
    12. PIR
    13. NCBI-taxon
    14. MIM
    15. UniGene
    16. PubMed
    17. EMBL
    18. EMBL-CDS
    19. Ensembl
    20. Ensembl_TRS
    21. Ensembl_PRO
    22. Additional PubMed

    We want column 1 (accession) and column 7 (GO terms, semicolon-separated)
    """
    logger.info(f"Parsing UniProt idmapping file: {idmapping_file}")

    if accession_filter:
        logger.info(f"Filtering to {len(accession_filter)} accessions")
    if taxon_filter:
        logger.info(f"Filtering to taxa: {taxon_filter}")

    # Determine if gzipped
    open_func = gzip.open if idmapping_file.endswith('.gz') else open
    mode = 'rt' if idmapping_file.endswith('.gz') else 'r'

    mapping = {}
    lines_read = 0
    entries_with_go = 0

    with open_func(idmapping_file, mode) as f:
        for line in f:
            lines_read += 1
            if lines_read % 10_000_000 == 0:
                logger.info(f"  Processed {lines_read:,} lines, found {entries_with_go:,} with GO terms")

            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue

            accession = parts[0]

            # Apply filters
            if accession_filter and accession not in accession_filter:
                continue

            if taxon_filter and len(parts) >= 13:
                taxon = parts[12]
                if taxon not in taxon_filter:
                    continue

            # Get GO terms (column 7, 0-indexed = 6)
            go_terms = parts[6].strip()
            if not go_terms:
                continue

            # GO terms are semicolon-separated
            go_list = [g.strip() for g in go_terms.split(';') if g.strip().startswith('GO:')]
            if go_list:
                mapping[accession] = go_list
                entries_with_go += 1

                # Also add UniProtKB-ID mapping (column 2)
                if len(parts) >= 2 and parts[1]:
                    uniprot_id = parts[1]
                    mapping[uniprot_id] = go_list

    logger.info(f"Parsed {lines_read:,} lines, found {len(mapping):,} entries with GO terms")

    # Write output
    write_mapping(mapping, output_file)
    return mapping


def parse_uniprot_dat(dat_file: str, output_file: str):
    """
    Parse UniProt .dat (Swiss-Prot flat file format) to extract GO terms.
    Useful for organism-specific downloads.
    """
    logger.info(f"Parsing UniProt DAT file: {dat_file}")

    open_func = gzip.open if dat_file.endswith('.gz') else open
    mode = 'rt' if dat_file.endswith('.gz') else 'r'

    mapping = {}
    current_acc = None
    current_go = []

    with open_func(dat_file, mode) as f:
        for line in f:
            if line.startswith('AC   '):
                # Accession line - may have multiple, take first
                accs = line[5:].strip().rstrip(';').split(';')
                if accs:
                    current_acc = accs[0].strip()

            elif line.startswith('DR   GO;'):
                # GO cross-reference
                # Format: DR   GO; GO:0005737; C:cytoplasm; IEA:UniProtKB-SubCell.
                parts = line[5:].split(';')
                if len(parts) >= 2:
                    go_term = parts[1].strip()
                    if go_term.startswith('GO:'):
                        current_go.append(go_term)

            elif line.startswith('//'):
                # End of entry
                if current_acc and current_go:
                    mapping[current_acc] = current_go
                current_acc = None
                current_go = []

    logger.info(f"Found {len(mapping)} entries with GO terms")
    write_mapping(mapping, output_file)
    return mapping


# =============================================================================
# GOA GAF FILE PARSER (RECOMMENDED FOR SWISS-PROT)
# =============================================================================

def parse_goa_gaf(gaf_file: str, output_file: str, taxon_filter: set = None,
                  evidence_filter: set = None):
    """
    Parse GO Annotation (GAF) file from EBI/UniProt.

    This is the RECOMMENDED method for Swiss-Prot GO mappings.

    Download from:
      - Swiss-Prot only: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/SWISSPROT/goa_swissprot.gaf.gz
      - All UniProt: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

    GAF 2.2 format columns (tab-separated):
      1. DB              - e.g., UniProtKB
      2. DB_Object_ID    - accession (e.g., P12345)
      3. DB_Object_Symbol - gene symbol
      4. Qualifier       - relationship to GO term
      5. GO_ID           - GO:0000001
      6. DB:Reference    - citation
      7. Evidence_Code   - IEA, IDA, IPI, etc.
      8. With/From       - supporting evidence
      9. Aspect          - F (function), P (process), C (component)
      10. DB_Object_Name  - protein name
      11. DB_Object_Synonym
      12. DB_Object_Type  - protein, gene, etc.
      13. Taxon           - taxon:9606
      14. Date
      15. Assigned_By
      16-17. Optional extension fields

    Args:
        gaf_file: Path to GAF file (.gaf or .gaf.gz)
        output_file: Output mapping file
        taxon_filter: Optional set of taxon IDs to include (e.g., {'9606', '10090'})
        evidence_filter: Optional set of evidence codes to include
                        (e.g., {'EXP', 'IDA', 'IPI'} for experimental only)
    """
    logger.info(f"Parsing GOA GAF file: {gaf_file}")

    if taxon_filter:
        logger.info(f"Filtering to taxa: {taxon_filter}")
    if evidence_filter:
        logger.info(f"Filtering to evidence codes: {evidence_filter}")
    open_func = gzip.open if gaf_file.endswith('.gz') else open
    mode = 'rt' if gaf_file.endswith('.gz') else 'r'

    mapping = defaultdict(set)
    lines_read = 0
    entries_processed = 0

    with open_func(gaf_file, mode) as f:
        for line in f:
            # Skip header lines
            if line.startswith('!'):
                continue

            lines_read += 1
            if lines_read % 1_000_000 == 0:
                logger.info(f"  Processed {lines_read:,} lines, {len(mapping):,} accessions with GO")

            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue

            db = parts[0]
            accession = parts[1]
            go_id = parts[4]
            evidence = parts[6] if len(parts) > 6 else ''
            taxon_field = parts[12] if len(parts) > 12 else ''

            # Only process UniProtKB entries
            if db != 'UniProtKB':
                continue

            # Apply evidence filter
            if evidence_filter and evidence not in evidence_filter:
                continue

            # Apply taxon filter
            if taxon_filter:
                # Taxon field format: "taxon:9606" or "taxon:9606|taxon:10090"
                taxons = re.findall(r'taxon:(\d+)', taxon_field)
                if not any(t in taxon_filter for t in taxons):
                    continue

            # Validate GO ID
            if not go_id.startswith('GO:'):
                continue

            mapping[accession].add(go_id)
            entries_processed += 1

    logger.info(f"Parsed {lines_read:,} lines, found {len(mapping):,} accessions with GO terms")
    logger.info(f"Total GO annotations: {entries_processed:,}")

    # Convert sets to lists for output
    mapping_list = {acc: list(terms) for acc, terms in mapping.items()}
    write_mapping(mapping_list, output_file)

    return mapping_list


# =============================================================================
# NCBI GENE2GO PARSER
# =============================================================================

def parse_ncbi_gene2go(gene2go_file: str, gene2accession_file: str,
                       output_file: str, taxon_filter: set = None):
    """
    Parse NCBI gene2go and gene2accession to create protein accession -> GO mapping.

    gene2go format:
    #tax_id GeneID  GO_ID   Evidence    Qualifier   GO_term PubMed  Category

    gene2accession format:
    #tax_id GeneID  status  RNA_nucleotide_accession.version    RNA_nucleotide_gi   protein_accession.version   protein_gi  genomic_nucleotide_accession.version    genomic_nucleotide_gi   start_position_on_the_genomic_accession end_position_on_the_genomic_accession   orientation assembly    mature_peptide_accession.version    mature_peptide_gi   Symbol
    """
    logger.info("Parsing NCBI gene2go and gene2accession files...")

    # First, parse gene2go to get GeneID -> GO terms
    logger.info(f"Reading gene2go: {gene2go_file}")
    gene_to_go = defaultdict(set)

    open_func = gzip.open if gene2go_file.endswith('.gz') else open
    mode = 'rt' if gene2go_file.endswith('.gz') else 'r'

    with open_func(gene2go_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            tax_id = parts[0]
            gene_id = parts[1]
            go_id = parts[2]

            if taxon_filter and tax_id not in taxon_filter:
                continue

            if go_id.startswith('GO:'):
                gene_to_go[gene_id].add(go_id)

    logger.info(f"Found GO terms for {len(gene_to_go)} genes")

    # Now parse gene2accession to map protein accessions to gene IDs
    logger.info(f"Reading gene2accession: {gene2accession_file}")
    mapping = {}
    lines_read = 0

    open_func = gzip.open if gene2accession_file.endswith('.gz') else open
    mode = 'rt' if gene2accession_file.endswith('.gz') else 'r'

    with open_func(gene2accession_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            lines_read += 1
            if lines_read % 5_000_000 == 0:
                logger.info(f"  Processed {lines_read:,} lines")

            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue

            tax_id = parts[0]
            gene_id = parts[1]
            protein_acc = parts[5]  # protein_accession.version

            if taxon_filter and tax_id not in taxon_filter:
                continue

            if protein_acc == '-' or not protein_acc:
                continue

            if gene_id in gene_to_go:
                go_terms = list(gene_to_go[gene_id])
                mapping[protein_acc] = go_terms
                # Also add without version
                acc_no_version = protein_acc.split('.')[0]
                mapping[acc_no_version] = go_terms

    logger.info(f"Created mapping for {len(mapping)} protein accessions")
    write_mapping(mapping, output_file)
    return mapping


# =============================================================================
# UNIPROT REST API
# =============================================================================

def query_uniprot_api(accessions: list, output_file: str, batch_size: int = 100):
    """
    Query UniProt REST API for GO terms.
    Works well for smaller numbers of accessions (< 10,000).
    """
    logger.info(f"Querying UniProt API for {len(accessions)} accessions")

    mapping = {}
    total = len(accessions)

    # Process in batches
    for i in range(0, total, batch_size):
        batch = accessions[i:i + batch_size]
        logger.info(f"Processing batch {i // batch_size + 1}/{(total + batch_size - 1) // batch_size}")

        # Query UniProt
        try:
            batch_results = _query_uniprot_batch(batch)
            mapping.update(batch_results)
        except Exception as e:
            logger.warning(f"Batch query failed: {e}")

        # Rate limiting
        time.sleep(0.5)

    logger.info(f"Retrieved GO terms for {len(mapping)} accessions")
    write_mapping(mapping, output_file)
    return mapping


def _query_uniprot_batch(accessions: list) -> dict:
    """Query UniProt for a batch of accessions."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # Build query
    acc_query = ' OR '.join(f'accession:{acc}' for acc in accessions)
    params = {
        'query': acc_query,
        'fields': 'accession,go_id',
        'format': 'tsv',
        'size': 500
    }

    url = f"{base_url}?{urlencode(params)}"

    req = Request(url, headers={'User-Agent': 'Python/GOAnnotate'})

    results = {}
    try:
        with urlopen(req, timeout=30) as response:
            content = response.read().decode('utf-8')

            for line in content.strip().split('\n')[1:]:  # Skip header
                parts = line.split('\t')
                if len(parts) >= 2:
                    accession = parts[0]
                    go_terms = [g.strip() for g in parts[1].split(';') if g.strip().startswith('GO:')]
                    if go_terms:
                        results[accession] = go_terms

    except (HTTPError, URLError) as e:
        logger.warning(f"API query failed: {e}")

    return results


def query_uniprot_api_streaming(accessions: list, output_file: str):
    """
    Use UniProt's ID mapping service for larger sets of accessions.
    This is more efficient for >1000 accessions.
    """
    logger.info(f"Using UniProt ID mapping service for {len(accessions)} accessions")

    # Submit job
    submit_url = "https://rest.uniprot.org/idmapping/run"

    data = urlencode({
        'from': 'UniProtKB_AC-ID',
        'to': 'UniProtKB',
        'ids': ','.join(accessions[:100000])  # Max 100k per request
    }).encode()

    req = Request(submit_url, data=data, headers={'User-Agent': 'Python/GOAnnotate'})

    try:
        with urlopen(req, timeout=60) as response:
            result = json.loads(response.read().decode())
            job_id = result.get('jobId')

        if not job_id:
            logger.error("Failed to submit ID mapping job")
            return {}

        logger.info(f"Submitted job: {job_id}")

        # Poll for results
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        while True:
            req = Request(status_url, headers={'User-Agent': 'Python/GOAnnotate'})
            with urlopen(req, timeout=30) as response:
                status = json.loads(response.read().decode())

            if 'jobStatus' in status:
                if status['jobStatus'] == 'FINISHED':
                    break
                elif status['jobStatus'] == 'ERROR':
                    logger.error("Job failed")
                    return {}

            logger.info("Waiting for job to complete...")
            time.sleep(5)

        # Get results
        results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}?fields=accession,go_id&format=tsv"
        req = Request(results_url, headers={'User-Agent': 'Python/GOAnnotate'})

        mapping = {}
        with urlopen(req, timeout=120) as response:
            content = response.read().decode('utf-8')
            for line in content.strip().split('\n')[1:]:
                parts = line.split('\t')
                if len(parts) >= 2:
                    accession = parts[0]
                    go_terms = [g.strip() for g in parts[1].split(';') if g.strip().startswith('GO:')]
                    if go_terms:
                        mapping[accession] = go_terms

        logger.info(f"Retrieved GO terms for {len(mapping)} accessions")
        write_mapping(mapping, output_file)
        return mapping

    except Exception as e:
        logger.error(f"ID mapping failed: {e}")
        return {}


# =============================================================================
# EXTRACT ACCESSIONS FROM BLAST RESULTS
# =============================================================================

def extract_accessions_from_blast(blast_file: str) -> set:
    """
    Extract unique accessions from BLAST results.
    Handles various formats: UniProt, RefSeq, GenBank, etc.
    """
    logger.info(f"Extracting accessions from BLAST results: {blast_file}")

    accessions = set()

    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            subject_id = parts[1]

            # Extract all possible accession formats
            extracted = _extract_accession_variants(subject_id)
            accessions.update(extracted)

    logger.info(f"Extracted {len(accessions)} unique accessions")
    return accessions


def _extract_accession_variants(subject_id: str) -> set:
    """Extract various accession formats from a subject ID."""
    accessions = set()

    # Add the full ID
    accessions.add(subject_id)

    # UniProt format: sp|P12345|GENE_SPECIES or tr|A0A123|GENE_SPECIES
    if '|' in subject_id:
        parts = subject_id.split('|')
        for part in parts:
            if part and part not in ('sp', 'tr', 'ref', 'gi', 'gb', 'emb', 'dbj', 'pir'):
                accessions.add(part)
                # Add without version
                if '.' in part:
                    accessions.add(part.split('.')[0])

    # RefSeq patterns: NP_123456.1, XP_123456.1, WP_123456.1
    for match in re.findall(r'[NXYW]P_\d+\.?\d*', subject_id):
        accessions.add(match)
        accessions.add(match.split('.')[0])

    # GenBank patterns: AAA12345.1
    for match in re.findall(r'[A-Z]{3}\d{5,}\.?\d*', subject_id):
        accessions.add(match)
        accessions.add(match.split('.')[0])

    return accessions


# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

def write_mapping(mapping: dict, output_file: str):
    """Write accession -> GO mapping to file."""
    logger.info(f"Writing mapping to {output_file}")

    with open(output_file, 'w') as f:
        for accession, go_terms in sorted(mapping.items()):
            if isinstance(go_terms, set):
                go_terms = list(go_terms)
            f.write(f"{accession}\t{';'.join(sorted(go_terms))}\n")

    logger.info(f"Wrote {len(mapping)} entries to {output_file}")


# =============================================================================
# CLI
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="Create GO term mappings for GOAnnotate",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Download source files
  python create_go_mapping.py download --source ncbi --output-dir ./db_files

  # Parse UniProt idmapping file
  python create_go_mapping.py uniprot-idmapping \\
      --idmapping idmapping_selected.tab.gz \\
      --output uniprot2go.tsv

  # Parse NCBI gene2go (requires gene2accession for protein mapping)
  python create_go_mapping.py ncbi-gene2go \\
      --gene2go gene2go.gz \\
      --gene2accession gene2accession.gz \\
      --output ncbi2go.tsv

  # Query UniProt API for accessions from BLAST results
  python create_go_mapping.py uniprot-api \\
      --blast-results blast_output.tsv \\
      --output blast2go.tsv

  # Filter to specific taxon (e.g., human = 9606)
  python create_go_mapping.py uniprot-idmapping \\
      --idmapping idmapping_selected.tab.gz \\
      --taxon 9606 \\
      --output human_uniprot2go.tsv
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Download command
    download_parser = subparsers.add_parser('download', help='Download source files')
    download_parser.add_argument(
        '--source', choices=['uniprot', 'ncbi', 'go'], required=True,
        help='Source to download'
    )
    download_parser.add_argument(
        '--output-dir', '-o', default='./',
        help='Output directory for downloaded files'
    )

    # UniProt idmapping command
    uniprot_parser = subparsers.add_parser(
        'uniprot-idmapping',
        help='Parse UniProt idmapping file'
    )
    uniprot_parser.add_argument(
        '--idmapping', required=True,
        help='Path to idmapping_selected.tab.gz'
    )
    uniprot_parser.add_argument(
        '--output', '-o', required=True,
        help='Output mapping file'
    )
    uniprot_parser.add_argument(
        '--taxon', nargs='+',
        help='Filter to specific taxon IDs (e.g., 9606 for human)'
    )
    uniprot_parser.add_argument(
        '--accessions-from-blast',
        help='Only include accessions found in this BLAST file'
    )

    # UniProt DAT file command
    dat_parser = subparsers.add_parser(
        'uniprot-dat',
        help='Parse UniProt .dat flat file'
    )
    dat_parser.add_argument(
        '--dat-file', required=True,
        help='Path to UniProt .dat or .dat.gz file'
    )
    dat_parser.add_argument(
        '--output', '-o', required=True,
        help='Output mapping file'
    )

    # GOA GAF file command (RECOMMENDED for Swiss-Prot)
    gaf_parser = subparsers.add_parser(
        'goa-gaf',
        help='Parse GOA GAF file (RECOMMENDED for Swiss-Prot)'
    )
    gaf_parser.add_argument(
        '--gaf', required=True,
        help='Path to GAF file (e.g., goa_swissprot.gaf.gz)'
    )
    gaf_parser.add_argument(
        '--output', '-o', required=True,
        help='Output mapping file'
    )
    gaf_parser.add_argument(
        '--taxon', nargs='+',
        help='Filter to specific taxon IDs (e.g., 9606 for human)'
    )
    gaf_parser.add_argument(
        '--evidence', nargs='+',
        help='Filter to specific evidence codes (e.g., EXP IDA IPI for experimental)'
    )

    # NCBI gene2go command
    ncbi_parser = subparsers.add_parser(
        'ncbi-gene2go',
        help='Parse NCBI gene2go file'
    )
    ncbi_parser.add_argument(
        '--gene2go', required=True,
        help='Path to gene2go.gz'
    )
    ncbi_parser.add_argument(
        '--gene2accession', required=True,
        help='Path to gene2accession.gz'
    )
    ncbi_parser.add_argument(
        '--output', '-o', required=True,
        help='Output mapping file'
    )
    ncbi_parser.add_argument(
        '--taxon', nargs='+',
        help='Filter to specific taxon IDs'
    )

    # UniProt API command
    api_parser = subparsers.add_parser(
        'uniprot-api',
        help='Query UniProt API for GO terms'
    )
    api_parser.add_argument(
        '--blast-results',
        help='BLAST results file to extract accessions from'
    )
    api_parser.add_argument(
        '--accessions-file',
        help='File with one accession per line'
    )
    api_parser.add_argument(
        '--output', '-o', required=True,
        help='Output mapping file'
    )
    api_parser.add_argument(
        '--batch-size', type=int, default=100,
        help='Batch size for API queries'
    )
    api_parser.add_argument(
        '--use-idmapping', action='store_true',
        help='Use ID mapping service (better for >1000 accessions)'
    )

    # Merge mappings command
    merge_parser = subparsers.add_parser(
        'merge',
        help='Merge multiple mapping files'
    )
    merge_parser.add_argument(
        '--inputs', nargs='+', required=True,
        help='Input mapping files to merge'
    )
    merge_parser.add_argument(
        '--output', '-o', required=True,
        help='Output merged mapping file'
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == 'download':
        download_sources(args.source, args.output_dir)

    elif args.command == 'uniprot-idmapping':
        accession_filter = None
        if args.accessions_from_blast:
            accession_filter = extract_accessions_from_blast(args.accessions_from_blast)

        taxon_filter = set(args.taxon) if args.taxon else None

        parse_uniprot_idmapping(
            args.idmapping,
            args.output,
            accession_filter=accession_filter,
            taxon_filter=taxon_filter
        )

    elif args.command == 'uniprot-dat':
        parse_uniprot_dat(args.dat_file, args.output)

    elif args.command == 'goa-gaf':
        taxon_filter = set(args.taxon) if args.taxon else None
        evidence_filter = set(args.evidence) if args.evidence else None
        parse_goa_gaf(
            args.gaf,
            args.output,
            taxon_filter=taxon_filter,
            evidence_filter=evidence_filter
        )

    elif args.command == 'ncbi-gene2go':
        taxon_filter = set(args.taxon) if args.taxon else None
        parse_ncbi_gene2go(
            args.gene2go,
            args.gene2accession,
            args.output,
            taxon_filter=taxon_filter
        )

    elif args.command == 'uniprot-api':
        # Get accessions
        accessions = []
        if args.blast_results:
            accessions = list(extract_accessions_from_blast(args.blast_results))
        elif args.accessions_file:
            with open(args.accessions_file) as f:
                accessions = [line.strip() for line in f if line.strip()]
        else:
            logger.error("Must provide --blast-results or --accessions-file")
            sys.exit(1)

        if args.use_idmapping:
            query_uniprot_api_streaming(accessions, args.output)
        else:
            query_uniprot_api(accessions, args.output, batch_size=args.batch_size)

    elif args.command == 'merge':
        # Merge multiple mapping files
        merged = defaultdict(set)
        for input_file in args.inputs:
            logger.info(f"Reading {input_file}")
            with open(input_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        accession = parts[0]
                        go_terms = [g.strip() for g in parts[1].split(';') if g.strip()]
                        merged[accession].update(go_terms)

        logger.info(f"Merged {len(merged)} accessions")
        write_mapping(dict(merged), args.output)

    logger.info("Done!")


if __name__ == "__main__":
    main()
