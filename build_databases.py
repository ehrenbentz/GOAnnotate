#!/usr/bin/env python3
"""
build_databases.py - Download and build all databases for the GOAnnotate pipeline.

Produces the following files in the specified output directory:

    Uniprot_All_GO_mapping.tsv      UniProt GO term mapping file
    SwissProt_Diamond_DB.dmnd       Swiss-Prot Diamond BLASTX database
    TrEMBL_Diamond_DB.dmnd          TrEMBL Diamond BLASTX database
    idmapping_selected.tab          UniProt ID mapping (for GOAnnotate --ncbi-idmapping)
    gene_info.tsv                   NCBI gene info (for GOAnnotate --ncbi-geneinfo)
    GOAnnotate_db.sqlite            SQLite index database (for GOAnnotate --db)

WARNING: This script downloads approximately 115 GB of compressed data and may
temporarily require more than 1 TB of disk space during processing. Expect the
full process to take several hours.

Usage:
    ./build_databases.py -o /path/to/databases
    ./build_databases.py -o /path/to/databases --threads 32
    ./build_databases.py -o /path/to/databases --keep-intermediates
"""

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent

# Download URLs
URLS = {
    "goa_gaf": (
        "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
    ),
    "sprot_fasta": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/complete/uniprot_sprot.fasta.gz"
    ),
    "trembl_fasta": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/complete/uniprot_trembl.fasta.gz"
    ),
    "idmapping": (
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        "knowledgebase/idmapping/idmapping_selected.tab.gz"
    ),
    "gene_info": (
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
    ),
}

TOTAL_STEPS = 6


# Utility functions
def format_elapsed(seconds):
    """Format elapsed seconds as a human-readable string."""
    h = int(seconds) // 3600
    m = (int(seconds) % 3600) // 60
    s = int(seconds) % 60
    return f"{h}h {m}m {s}s"


def format_size(nbytes):
    """Format a byte count as a human-readable string."""
    if nbytes >= 1024 ** 3:
        return f"{nbytes / (1024 ** 3):.1f} GB"
    if nbytes >= 1024 ** 2:
        return f"{nbytes / (1024 ** 2):.0f} MB"
    return f"{nbytes / 1024:.0f} KB"


def disk_free_gb(path):
    """Return available disk space in GB for the filesystem containing path."""
    st = os.statvfs(path)
    return (st.f_bavail * st.f_frsize) / (1024 ** 3)


def check_dependencies():
    """Verify that all required external tools and helper scripts are available."""
    missing = []

    for cmd in ["wget", "diamond", "gunzip"]:
        if not shutil.which(cmd):
            missing.append(cmd)

    if missing:
        logger.error("Missing required dependencies:")
        for m in missing:
            logger.error(f"  - {m}")
        sys.exit(1)

    logger.info("All dependencies found.")


def step_banner(step_num, title):
    """Print a prominent step header."""
    logger.info("=" * 72)
    logger.info(f"STEP {step_num}/{TOTAL_STEPS}: {title}")
    logger.info("=" * 72)


def download(url, outdir, label):
    """Download a file with wget. Supports resume via wget -c.

    Files are saved using the filename from the URL. Returns the path
    to the downloaded file.
    """
    filename = url.rsplit("/", 1)[-1]
    dest = outdir / filename

    if dest.exists():
        logger.info(f"[SKIP] {label}: {filename} already exists.")
        return dest

    logger.info(f"Downloading {label} ...")
    logger.info(f"  URL:  {url}")
    logger.info(f"  Dest: {dest}")

    subprocess.run(
        ["wget", "-c", url],
        check=True,
        cwd=str(outdir),
    )

    if not dest.exists():
        logger.error(f"Expected file not found after download: {dest}")
        sys.exit(1)

    return dest


def decompress(gz_path, output_path):
    """Decompress a .gz file to the specified output path using gunzip."""
    if output_path.exists():
        logger.info(f"[SKIP] {output_path.name} already exists.")
        return

    if not gz_path.exists():
        logger.error(f"Cannot decompress: {gz_path} not found.")
        sys.exit(1)

    logger.info(f"Decompressing {gz_path.name} -> {output_path.name} ...")

    with open(output_path, "wb") as fout:
        subprocess.run(["gunzip", "-c", str(gz_path)], stdout=fout, check=True)

    logger.info(f"  Decompressed: {output_path.name} ({format_size(output_path.stat().st_size)})")


def cleanup(path, keep):
    """Remove a file unless keep is True."""
    if keep or not path.exists():
        return
    path.unlink()
    logger.info(f"  Cleaned up intermediate: {path.name}")


# Create GO mapping file from GOA GAF file
def parse_goa_gaf(gaf_file: str, output_file: str, taxon_filter: set = None,
                  evidence_filter: set = None):
    """
    Parse GO Annotation (GAF) file from EBI/UniProt.

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

    Returns:
        dict of accession -> list of GO terms
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


def write_mapping(mapping: dict, output_file: str):
    """Write accession -> GO mapping to file (accession<TAB>GO:xxxx;GO:yyyy)."""
    logger.info(f"Writing mapping to {output_file}")

    with open(output_file, 'w') as f:
        for accession, go_terms in sorted(mapping.items()):
            if isinstance(go_terms, set):
                go_terms = list(go_terms)
            f.write(f"{accession}\t{';'.join(sorted(go_terms))}\n")

    logger.info(f"Wrote {len(mapping)} entries to {output_file}")


# Main pipeline
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download and build all databases required for the GOAnnotate "
            "gene annotation pipeline."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory to store all database files.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 8,
        help="Number of threads for diamond makedb (default: all available CPUs).",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep compressed download files after processing.",
    )
    parser.add_argument(
        "--rebuild-index",
        action="store_true",
        help="Rebuild the SQLite index database even if it already exists.",
    )
    args = parser.parse_args()

    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    keep = args.keep_intermediates
    threads = str(args.threads)

    # Warning and confirmation
    print()
    print("=" * 72)
    print("  GOAnnotate Database Builder")
    print("=" * 72)
    print()
    print("  WARNING: This process will:")
    print("    - Download ~115 GB of compressed data")
    print("    - Temporarily require more than 1 TB of disk space")
    print("    - Take 12-24+ hours to complete")
    print()
    print(f"  Output directory  : {outdir}")
    print(f"  Available disk    : {disk_free_gb(str(outdir)):.0f} GB")
    print(f"  Diamond threads   : {args.threads}")
    print(f"  Keep intermediates: {keep}")
    print()

    response = input("  Proceed? [y/N] ").strip().lower()
    if response != "y":
        print("  Aborted.")
        sys.exit(0)

    print()
    pipeline_start = time.time()

    # Add file handler so logs are saved alongside the databases
    file_handler = logging.FileHandler(str(outdir / "build_databases.log"))
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logging.getLogger().addHandler(file_handler)
    logger.info(f"Logging to {outdir / 'build_databases.log'}")

    check_dependencies()

    # Define output file paths
    go_mapping   = outdir / "Uniprot_All_GO_mapping.tsv"
    sprot_dmnd   = outdir / "SwissProt_Diamond_DB.dmnd"
    trembl_dmnd  = outdir / "TrEMBL_Diamond_DB.dmnd"
    idmap_tab    = outdir / "idmapping_selected.tab"
    geneinfo_tsv = outdir / "gene_info.tsv"
    sqlite_db    = outdir / "GOAnnotate_db.sqlite"

    # UniProt GO Mapping
    step_banner(1, "UniProt GO Mapping (Uniprot_All_GO_mapping.tsv)")
    t0 = time.time()

    if go_mapping.exists():
        logger.info(f"[SKIP] {go_mapping.name} already exists.")
    else:
        gaf_gz = download(URLS["goa_gaf"], outdir, "GOA GAF file (~30 GB)")

        logger.info("Generating GO mapping from GAF file ...")
        gaf_mapping = parse_goa_gaf(str(gaf_gz), str(go_mapping))
        logger.info(f"GO mapping complete: {len(gaf_mapping):,} accessions")

        cleanup(gaf_gz, keep)

    logger.info(f"Step 1 elapsed: {format_elapsed(time.time() - t0)}")

    # Swiss-Prot Diamond Database
    step_banner(2, "Swiss-Prot Diamond Database (SwissProt_Diamond_DB.dmnd)")
    t0 = time.time()

    if sprot_dmnd.exists():
        logger.info(f"[SKIP] {sprot_dmnd.name} already exists.")
    else:
        sprot_gz = download(URLS["sprot_fasta"], outdir, "Swiss-Prot FASTA (~90 MB)")

        logger.info("Building Swiss-Prot Diamond database ...")
        subprocess.run(
            [
                "diamond", "makedb",
                "--in", str(sprot_gz),
                "-d", str(outdir / "SwissProt_Diamond_DB"),
                "--threads", threads,
            ],
            check=True,
        )

        cleanup(sprot_gz, keep)

    logger.info(f"Step 2 elapsed: {format_elapsed(time.time() - t0)}")

    # TrEMBL Diamond Database
    step_banner(3, "TrEMBL Diamond Database (TrEMBL_Diamond_DB.dmnd)")
    t0 = time.time()

    if trembl_dmnd.exists():
        logger.info(f"[SKIP] {trembl_dmnd.name} already exists.")
    else:
        trembl_gz = download(URLS["trembl_fasta"], outdir, "TrEMBL FASTA (~55 GB)")

        logger.info("Building TrEMBL Diamond database ...")
        subprocess.run(
            [
                "diamond", "makedb",
                "--in", str(trembl_gz),
                "-d", str(outdir / "TrEMBL_Diamond_DB"),
                "--threads", threads,
            ],
            check=True,
        )

        cleanup(trembl_gz, keep)

    logger.info(f"Step 3 elapsed: {format_elapsed(time.time() - t0)}")

    # UniProt ID Mapping
    step_banner(4, "UniProt ID Mapping (idmapping_selected.tab)")
    t0 = time.time()

    if idmap_tab.exists():
        logger.info(f"[SKIP] {idmap_tab.name} already exists.")
    else:
        idmap_gz = download(URLS["idmapping"], outdir, "UniProt ID mapping (~25 GB)")

        decompress(idmap_gz, idmap_tab)
        cleanup(idmap_gz, keep)

    logger.info(f"Step 4 elapsed: {format_elapsed(time.time() - t0)}")

    # NCBI Gene Info
    step_banner(5, "NCBI Gene Info (gene_info.tsv)")
    t0 = time.time()

    if geneinfo_tsv.exists():
        logger.info(f"[SKIP] {geneinfo_tsv.name} already exists.")
    else:
        gi_gz = download(URLS["gene_info"], outdir, "NCBI gene_info (~3 GB)")

        decompress(gi_gz, geneinfo_tsv)
        cleanup(gi_gz, keep)

    logger.info(f"Step 5 elapsed: {format_elapsed(time.time() - t0)}")

    # SQLite Index Database
    step_banner(6, "SQLite Index Database (GOAnnotate_db.sqlite)")
    t0 = time.time()

    if sqlite_db.exists() and not args.rebuild_index:
        logger.info(f"[SKIP] {sqlite_db.name} already exists. "
                     "Use --rebuild-index to rebuild.")
    else:
        if sqlite_db.exists():
            sqlite_db.unlink()
            logger.info(f"Removed existing {sqlite_db.name} for rebuild.")

        import sqlite3

        conn = sqlite3.connect(str(sqlite_db))
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=OFF")

        BATCH_SIZE = 100_000

        # go_mapping
        if go_mapping.exists():
            logger.info(f"Building go_mapping table from {go_mapping.name} ...")
            conn.execute(
                "CREATE TABLE go_mapping "
                "(accession TEXT, go_terms TEXT)"
            )
            batch = []
            rows = 0
            with open(go_mapping) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    accession = parts[0].strip()
                    go_terms = '\t'.join(parts[1:])
                    batch.append((accession, go_terms))
                    if len(batch) >= BATCH_SIZE:
                        conn.executemany(
                            "INSERT INTO go_mapping VALUES (?, ?)", batch
                        )
                        rows += len(batch)
                        logger.info(f"  go_mapping: {rows:,} rows inserted...")
                        batch = []
            if batch:
                conn.executemany(
                    "INSERT INTO go_mapping VALUES (?, ?)", batch
                )
                rows += len(batch)
            conn.execute(
                "CREATE INDEX idx_go_mapping_acc ON go_mapping(accession)"
            )
            conn.commit()
            logger.info(f"  go_mapping: {rows:,} total rows, indexed.")
        else:
            logger.warning(
                f"  {go_mapping.name} not found, skipping go_mapping table."
            )

        # ncbi_idmapping
        if idmap_tab.exists():
            logger.info(
                f"Building ncbi_idmapping table from {idmap_tab.name} ..."
            )
            conn.execute(
                "CREATE TABLE ncbi_idmapping "
                "(uniprot_acc TEXT, gene_id TEXT)"
            )
            batch = []
            rows = 0
            with open(idmap_tab) as f:
                for line in f:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    accession = parts[0].strip()
                    gene_id = parts[2].strip()
                    if not gene_id:
                        continue
                    batch.append((accession, gene_id))
                    if len(batch) >= BATCH_SIZE:
                        conn.executemany(
                            "INSERT INTO ncbi_idmapping VALUES (?, ?)", batch
                        )
                        rows += len(batch)
                        logger.info(
                            f"  ncbi_idmapping: {rows:,} rows inserted..."
                        )
                        batch = []
            if batch:
                conn.executemany(
                    "INSERT INTO ncbi_idmapping VALUES (?, ?)", batch
                )
                rows += len(batch)
            conn.execute(
                "CREATE INDEX idx_ncbi_idmapping_acc "
                "ON ncbi_idmapping(uniprot_acc)"
            )
            conn.commit()
            logger.info(f"  ncbi_idmapping: {rows:,} total rows, indexed.")
        else:
            logger.warning(
                f"  {idmap_tab.name} not found, "
                "skipping ncbi_idmapping table."
            )

        # ncbi_geneinfo
        if geneinfo_tsv.exists():
            logger.info(
                f"Building ncbi_geneinfo table from {geneinfo_tsv.name} ..."
            )
            conn.execute(
                "CREATE TABLE ncbi_geneinfo "
                "(gene_id TEXT, symbol TEXT)"
            )
            batch = []
            rows = 0
            with open(geneinfo_tsv) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    gene_id = parts[1].strip()
                    symbol = parts[2].strip()
                    if symbol in ('NEWENTRY', '-'):
                        continue
                    batch.append((gene_id, symbol))
                    if len(batch) >= BATCH_SIZE:
                        conn.executemany(
                            "INSERT INTO ncbi_geneinfo VALUES (?, ?)", batch
                        )
                        rows += len(batch)
                        logger.info(
                            f"  ncbi_geneinfo: {rows:,} rows inserted..."
                        )
                        batch = []
            if batch:
                conn.executemany(
                    "INSERT INTO ncbi_geneinfo VALUES (?, ?)", batch
                )
                rows += len(batch)
            conn.execute(
                "CREATE INDEX idx_ncbi_geneinfo_id "
                "ON ncbi_geneinfo(gene_id)"
            )
            conn.commit()
            logger.info(f"  ncbi_geneinfo: {rows:,} total rows, indexed.")
        else:
            logger.warning(
                f"  {geneinfo_tsv.name} not found, "
                "skipping ncbi_geneinfo table."
            )

        logger.info("Running VACUUM on SQLite database ...")
        conn.execute("VACUUM")
        conn.close()

        logger.info(
            f"  SQLite database: {sqlite_db.name} "
            f"({format_size(sqlite_db.stat().st_size)})"
        )

    logger.info(f"Step 6 elapsed: {format_elapsed(time.time() - t0)}")

    # Summary
    total_elapsed = time.time() - pipeline_start

    print()
    print("=" * 72)
    print("  BUILD COMPLETE")
    print("=" * 72)
    print()
    print("  Output files:")

    outputs = [go_mapping, sprot_dmnd, trembl_dmnd, idmap_tab, geneinfo_tsv,
               sqlite_db]
    all_ok = True
    for f in outputs:
        if f.exists():
            print(f"    [OK]  {f.name}  ({format_size(f.stat().st_size)})")
        else:
            print(f"    [!!]  {f.name}  MISSING")
            all_ok = False

    print()
    print(f"  Total elapsed time  : {format_elapsed(total_elapsed)}")
    print(f"  Remaining disk space: {disk_free_gb(str(outdir)):.0f} GB")
    print()

    if all_ok:
        print("  Use these files with GOAnnotate:")
        print(f"    --db             {sqlite_db}  (recommended)")
        print()
        print("  Or with individual flat files:")
        print(f"    --go-mapping     {go_mapping}")
        print(f"    --ncbi-idmapping {idmap_tab}")
        print(f"    --ncbi-geneinfo  {geneinfo_tsv}")
        print()
        print("  Diamond databases for BLASTX:")
        print(f"    Swiss-Prot: {sprot_dmnd}")
        print(f"    TrEMBL:     {trembl_dmnd}")
        print()

    logger.info("Done.")


if __name__ == "__main__":
    main()
