#!/usr/bin/env python3
"""
build_databases.py - Download, verify and build the reference databases for GOAnnotate.

SQLite database (built unless --skip-sqlite), GOAnnotate_db.sqlite:

    go_mapping          UniProt accession -> GO terms, each with
                        the evidence codes behind it
                        (GO:0005515:IDA,IPI;GO:0005634:IEA)     GOA goa_uniprot_all.gaf.gz
    ncbi_idmapping      UniProt accession -> NCBI GeneID        UniProt idmapping_selected.tab.gz
    ncbi_geneinfo       GeneID -> official symbol               NCBI gene_info.gz
    ncbi_gene_symbols   (taxon, name) -> GeneID, where the name
                        is the official symbol, the locus tag
                        or a synonym (column kind)              NCBI gene_info.gz
    ncbi_orthologs      GeneID -> ortholog GeneID               NCBI gene_orthologs.gz
    taxonomy_nodes      taxon -> parent, rank, name             NCBI taxdump.tar.gz
    taxonomy_merged     old taxon ID -> current taxon ID        NCBI taxdump.tar.gz
    go_taxon_constraints GO term -> in_taxon / never_in_taxon   GO go-computed-taxon-constraints.obo
    go_taxon_unions     GO taxon union class -> member taxa     GO go-taxon-groupings.obo
    sources             one row per downloaded file: URL, size, MD5, how it was
                        verified, release, server date, download date
    metadata            build date, script version, UniProt release, row counts,
                        protein databases built

    --tables go,taxonomy builds only those table groups (go, idmapping,
    gene_info, orthologs, taxonomy); with --update they replace the same
    tables in an existing database, which is otherwise never touched.

Protein databases (with --protein-dbs LINEAGE), for the Diamond searches:

    <lineage>_SwissProt.fasta / .dmnd     UniProtKB/Swiss-Prot entries of the lineage
    <lineage>_TrEMBL.fasta / .dmnd        UniProtKB/TrEMBL entries of the lineage
    <lineage>_*_taxonmap.tsv.gz           accession -> taxon ID maps
    taxdump/                              NCBI taxonomy used for the lineage and the maps
    <lineage>_UniProtKB.dmnd              combined database (only with --combined-db)

    LINEAGE is an NCBI taxon name (metazoa, vertebrata, viridiplantae, ...) or
    a taxon ID (33208); "all" keeps every entry. Several lineages are given
    comma-separated (metazoa,viridiplantae) and built from one download. Each
    lineage FASTA records the UniProt release it was built from
    (<file>.release); when a later run downloads another release, the FASTA,
    its taxon map and its Diamond database are rebuilt. Entries are assigned through
    the OX= taxon of each UniProt header and the NCBI taxonomy tree (merged
    taxon IDs resolved). The Diamond databases carry taxonomy, so searches can
    be restricted further with --taxonlist.

Verification. Every file is downloaded to <name>.part (resumable) and only
renamed when complete. It is then checked against the checksum the provider
publishes: the UniProt RELEASE.metalink (MD5 and size) for the UniProt files,
the .md5 file beside each GOA and NCBI taxonomy file. NCBI gene/DATA publishes
no checksums; gene_info.gz and gene_orthologs.gz are checked against the
server's Content-Length and by decompressing them completely. A verified file
gets a <name>.verified.json record; a file without one is never trusted, so
an interrupted download is never used. The SQLite database is written under
a temporary name and renamed only when every table is complete.

Usage:
    ./build_databases.py -o /path/to/databases --yes
    ./build_databases.py -o /path/to/databases --protein-dbs metazoa --threads 64 --yes
    ./build_databases.py -o /path/to/databases --protein-dbs 7742 --combined-db --keep-downloads
    ./build_databases.py -o /path/to/databases --protein-dbs metazoa,viridiplantae --skip-sqlite --yes
    ./build_databases.py -o /path/to/databases --tables go,taxonomy --update --yes

Disk: about 80 GB of downloads (GAF 29 GB, ID mapping 7 GB, gene_info 1.6 GB;
TrEMBL FASTA 41 GB with --protein-dbs), plus the SQLite database (~30 GB) and,
with --protein-dbs metazoa, ~45 GB of FASTA and Diamond databases. The GAF is
sorted by accession on disk (GNU sort; about 60 GB of temporary files) and
written one accession at a time.
"""

import argparse
import datetime
import email.utils
import gzip
import hashlib
import json
import logging
import os
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import time
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

SCRIPT_VERSION = "1.0.0"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

UNIPROT_KB = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase"
GOA_UNIPROT = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT"
NCBI_GENE = "https://ftp.ncbi.nlm.nih.gov/gene/DATA"
NCBI_TAXONOMY = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
GO_IMPORTS = "http://current.geneontology.org/ontology/imports"

# name -> (URL, checksum method, checksum location)
#   metalink: MD5 and size from a UniProt RELEASE.metalink
#   md5file:  MD5 from a file beside the download
#   none:     no published checksum; Content-Length and full decompression
SOURCES = {
    "goa_gaf": (f"{GOA_UNIPROT}/goa_uniprot_all.gaf.gz", "md5file",
                f"{GOA_UNIPROT}/goa_uniprot_all.gaf.gz.md5"),
    "idmapping": (f"{UNIPROT_KB}/idmapping/idmapping_selected.tab.gz", "metalink",
                  f"{UNIPROT_KB}/idmapping/RELEASE.metalink"),
    "gene_info": (f"{NCBI_GENE}/gene_info.gz", "none", None),
    "gene_orthologs": (f"{NCBI_GENE}/gene_orthologs.gz", "none", None),
    "sprot_fasta": (f"{UNIPROT_KB}/complete/uniprot_sprot.fasta.gz", "metalink",
                    f"{UNIPROT_KB}/complete/RELEASE.metalink"),
    "trembl_fasta": (f"{UNIPROT_KB}/complete/uniprot_trembl.fasta.gz", "metalink",
                     f"{UNIPROT_KB}/complete/RELEASE.metalink"),
    "taxdump": (f"{NCBI_TAXONOMY}/taxdump.tar.gz", "md5file",
                f"{NCBI_TAXONOMY}/taxdump.tar.gz.md5"),
    "go_taxon_constraints": (f"{GO_IMPORTS}/go-computed-taxon-constraints.obo", "none", None),
    "go_taxon_groupings": (f"{GO_IMPORTS}/go-taxon-groupings.obo", "none", None),
}
UNIPROT_RELDATE = f"{UNIPROT_KB}/complete/reldate.txt"

# Tables of the SQLite database by build step (--tables)
TABLE_GROUPS = {
    "go": ["go_mapping"],
    "idmapping": ["ncbi_idmapping"],
    "gene_info": ["ncbi_geneinfo", "ncbi_gene_symbols"],
    "orthologs": ["ncbi_orthologs"],
    "taxonomy": ["taxonomy_nodes", "taxonomy_merged", "go_taxon_constraints", "go_taxon_unions"],
}
SOURCES_FOR_GROUP = {
    "go": ["goa_gaf"], "idmapping": ["idmapping"], "gene_info": ["gene_info"], "orthologs": ["gene_orthologs"],
    "taxonomy": ["taxdump", "go_taxon_constraints", "go_taxon_groupings"],
}

TAXONMAP_HEADER = "accession.version\ttaxid\n"
BATCH_SIZE = 200_000


####################
# UTILITIES
####################

def format_elapsed(seconds):
    h = int(seconds) // 3600
    m = (int(seconds) % 3600) // 60
    s = int(seconds) % 60
    return f"{h}h {m}m {s}s"


def format_size(nbytes):
    if nbytes >= 1024 ** 3:
        return f"{nbytes / (1024 ** 3):.1f} GB"
    if nbytes >= 1024 ** 2:
        return f"{nbytes / (1024 ** 2):.0f} MB"
    return f"{nbytes / 1024:.0f} KB"


def disk_free_gb(path):
    st = os.statvfs(path)
    return (st.f_bavail * st.f_frsize) / (1024 ** 3)


def now_iso():
    return datetime.datetime.now().isoformat(timespec="seconds")


class Steps:
    """Numbered step banners with elapsed time."""

    def __init__(self, titles):
        self.titles = titles
        self.n = 0
        self.t0 = None

    def start(self, title):
        self.finish()
        self.n += 1
        self.t0 = time.time()
        logger.info("=" * 72)
        logger.info(f"STEP {self.n}/{len(self.titles)}: {title}")
        logger.info("=" * 72)

    def finish(self):
        if self.t0 is not None:
            logger.info(f"Step {self.n} elapsed: {format_elapsed(time.time() - self.t0)}")
            self.t0 = None


def check_dependencies(need_diamond):
    missing = [cmd for cmd in ["wget", "md5sum"] if not shutil.which(cmd)]
    if need_diamond and not shutil.which("diamond"):
        missing.append("diamond")
    if missing:
        logger.error("Missing required programs: " + ", ".join(missing))
        sys.exit(1)
    logger.info(f"Decompression: {'pigz' if shutil.which('pigz') else 'gzip'}")
    if need_diamond:
        out = subprocess.run(["diamond", "version"], capture_output=True, text=True)
        logger.info(f"Diamond: {out.stdout.strip() or out.stderr.strip()}")


def decompressor():
    return "pigz" if shutil.which("pigz") else "gzip"


def open_gz_stream(path):
    """Text stream of a .gz file through pigz/gzip (much faster than the gzip module)."""
    proc = subprocess.Popen([decompressor(), "-dc", str(path)], stdout=subprocess.PIPE,
                            bufsize=16 * 1024 * 1024)
    return proc


def close_gz_stream(proc, path):
    proc.stdout.close()
    if proc.wait() != 0:
        logger.error(f"Decompression of {path} failed")
        sys.exit(1)


####################
# REMOTE METADATA
####################

def http_get_text(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": f"GOAnnotate-build_databases/{SCRIPT_VERSION}"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def http_head(url, timeout=120):
    """(size, last_modified ISO date) from the server, or (None, '')."""
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent": f"GOAnnotate-build_databases/{SCRIPT_VERSION}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            size = r.headers.get("Content-Length")
            lm = r.headers.get("Last-Modified", "")
    except Exception as e:  # noqa: BLE001
        logger.warning(f"  HEAD request failed for {url}: {e}")
        return None, ""
    if lm:
        try:
            lm = email.utils.parsedate_to_datetime(lm).date().isoformat()
        except (TypeError, ValueError):
            pass
    return (int(size) if size and size.isdigit() else None), lm


_metalink_cache = {}


def metalink_entry(metalink_url, filename):
    """(md5, size, release version) for a file in a UniProt RELEASE.metalink."""
    if metalink_url not in _metalink_cache:
        _metalink_cache[metalink_url] = ET.fromstring(http_get_text(metalink_url))
    root = _metalink_cache[metalink_url]
    ns = {"m": "http://www.metalinker.org/"}
    version = root.findtext("m:version", default="", namespaces=ns)
    for f in root.iter("{http://www.metalinker.org/}file"):
        if f.get("name") == filename:
            size = f.findtext("m:size", default="", namespaces=ns)
            md5 = ""
            for h in f.iter("{http://www.metalinker.org/}hash"):
                if h.get("type") == "md5":
                    md5 = (h.text or "").strip()
            return md5, (int(size) if size.isdigit() else None), version
    return "", None, version


def expected_for(name):
    """Expected checksum information for a source, fetched from the provider."""
    url, method, where = SOURCES[name]
    filename = url.rsplit("/", 1)[-1]
    size, last_modified = http_head(url)
    info = {"name": name, "url": url, "file": filename, "method": method, "md5": "",
            "size": size, "last_modified": last_modified, "release": ""}
    if method == "metalink":
        md5, msize, version = metalink_entry(where, filename)
        info.update(md5=md5, release=f"UniProt {version}" if version else "",
                    checksum_source=where)
        if msize:
            info["size"] = msize
        if not md5:
            logger.warning(f"  {filename} not found in {where}; falling back to the integrity test")
            info["method"] = "none"
    elif method == "md5file":
        text = http_get_text(where).strip()
        info.update(md5=text.split()[0] if text else "", checksum_source=where)
        if not info["md5"]:
            info["method"] = "none"
    if info["method"] == "none":
        info["checksum_source"] = "none published; Content-Length and full decompression test"
    if not info["release"] and last_modified:
        info["release"] = f"server date {last_modified}"
    return info


####################
# DOWNLOAD AND VERIFY
####################

def md5_of(path):
    out = subprocess.run(["md5sum", str(path)], capture_output=True, text=True, check=True)
    return out.stdout.split()[0]


def gzip_test(path):
    """True when the whole file decompresses without error."""
    if str(path).endswith(".tar.gz") or str(path).endswith(".gz"):
        r = subprocess.run([decompressor(), "-t", str(path)], capture_output=True, text=True)
        return r.returncode == 0
    return True


def verify(path, info):
    """(ok, how) for a local file against the expected information."""
    size = path.stat().st_size
    if info.get("size") and size != info["size"]:
        return False, f"size {size:,} differs from the expected {info['size']:,}"
    if info["method"] in ("metalink", "md5file") and info["md5"]:
        logger.info(f"  Computing MD5 of {path.name} ({format_size(size)}) ...")
        got = md5_of(path)
        if got != info["md5"]:
            return False, f"MD5 {got} differs from the published {info['md5']}"
        return True, f"MD5 {got} matches {info['checksum_source']}"
    if info.get("file", "").endswith(".obo"):
        with open(path, errors="replace") as fh:
            head = fh.read(4096)
        if "format-version:" not in head:
            return False, "not an OBO file (no format-version header)"
        how = "OBO header present"
        if info.get("size"):
            how += "; size matches the server's Content-Length"
        return True, how
    logger.info(f"  Testing decompression of {path.name} ({format_size(size)}) ...")
    if not gzip_test(path):
        return False, "decompression test failed (truncated or corrupt)"
    how = "complete decompression test passed"
    if info.get("size"):
        how += "; size matches the server's Content-Length"
    return True, how


def fetch(name, outdir, allow_existing=True):
    """
    Download one source and verify it. Returns the info dict (with path,
    md5, verification and download date) and writes <file>.verified.json.
    """
    info = expected_for(name)
    dest = outdir / info["file"]
    part = outdir / (info["file"] + ".part")
    record = outdir / (info["file"] + ".verified.json")

    logger.info(f"{info['file']}: {info['release'] or 'release unknown'}; expected "
                f"{format_size(info['size']) if info.get('size') else 'size unknown'}; "
                f"checksum: {info['checksum_source']}")

    if dest.exists() and allow_existing:
        if record.exists():
            saved = json.loads(record.read_text())
            outdated = (info["md5"] and saved.get("md5_expected") and saved["md5_expected"] != info["md5"]) or \
                (not info["md5"] and info["last_modified"] and saved.get("server_date")
                 and saved["server_date"] != info["last_modified"])
            if outdated:
                logger.warning(f"  {dest.name} was verified against an earlier release "
                               f"({saved.get('release')}); the server now has {info['release']}. "
                               f"Downloading the current file.")
                dest.unlink()
                record.unlink()
            else:
                logger.info(f"  [SKIP] {dest.name} already verified ({saved.get('verified_by')})")
                saved["path"] = str(dest)
                return saved
        else:
            logger.info(f"  {dest.name} exists without a verification record; verifying it")
            ok, how = verify(dest, info)
            if ok:
                return _write_record(dest, record, info, how, downloaded=None)
            logger.warning(f"  {dest.name}: {how}; downloading again")
            dest.unlink()

    for attempt in (1, 2):
        logger.info(f"  Downloading {info['url']}")
        if shutil.which("aria2c"):
            # Several connections per file: EBI and UniProt limit the rate of
            # a single connection (about 1 MB/s against 10 MB/s with eight).
            # -c resumes a partial download, including one written by wget.
            cmd = ["aria2c", "-c", "-x", "8", "-s", "8", "-k", "64M", "--file-allocation=none",
                   "--max-tries=20", "--retry-wait=30", "--timeout=120", "--summary-interval=300",
                   "--console-log-level=warn", "--allow-overwrite=false", "--auto-file-renaming=false",
                   "-d", str(part.parent), "-o", part.name, info["url"]]
            tool = "aria2c"
        else:
            cmd = ["wget", "-c", "--progress=dot:giga", "--tries=20", "--waitretry=30",
                   "--read-timeout=120", "-O", str(part), info["url"]]
            tool = "wget"
        r = subprocess.run(cmd)
        if r.returncode != 0:
            logger.error(f"  {tool} failed with exit code {r.returncode}; the partial file is kept as "
                         f"{part.name} and the download resumes on the next run")
            sys.exit(1)
        ok, how = verify(part, info)
        if ok:
            part.rename(dest)
            logger.info(f"  Verified: {how}")
            return _write_record(dest, record, info, how, downloaded=now_iso())
        logger.warning(f"  Attempt {attempt}: {how}")
        part.unlink()
        control = part.with_name(part.name + ".aria2")
        if control.exists():
            control.unlink()
    logger.error(f"  {info['file']} could not be downloaded intact after two attempts")
    sys.exit(1)


def _write_record(dest, record, info, how, downloaded):
    rec = {
        "name": info["name"], "url": info["url"], "file": info["file"], "path": str(dest),
        "size": dest.stat().st_size, "md5_expected": info["md5"],
        "checksum_source": info["checksum_source"], "verified_by": how,
        "release": info["release"], "server_date": info["last_modified"],
        "downloaded": downloaded or datetime.datetime.fromtimestamp(dest.stat().st_mtime).isoformat(
            timespec="seconds"),
        "verified": now_iso(),
    }
    record.write_text(json.dumps(rec, indent=2) + "\n")
    return rec


####################
# SQLITE TABLES
####################

def _format_go_terms(terms):
    """GO:0005515:IDA,IPI;GO:0005634:IEA for a {term: set of codes} dict."""
    return ";".join(f"{t}:{','.join(sorted(c))}" if c else t for t, c in sorted(terms.items()))


def sorted_gaf_stream(gaf_gz, tmpdir, threads):
    """
    (accession, GO ID, evidence code) lines of the UniProtKB annotations in the
    GAF without NOT qualifiers, sorted by accession with GNU sort on disk.
    The GAF is not grouped by accession (annotations of one accession are
    spread over the file), so it is sorted before accessions are written.
    Returns (shell process, counter file with the number of NOT lines).
    """
    not_count = Path(tmpdir) / "gaf_not_count.txt"
    awk = shutil.which("mawk") or "awk"
    cmd = (f"{decompressor()} -dc {shlex.quote(str(gaf_gz))} | "
           f"{awk} -F'\\t' -v OFS='\\t' -v nf={shlex.quote(str(not_count))} "
           f"'$1==\"UniProtKB\" && $5 ~ /^GO:/ {{ if ($4 ~ /NOT/) {{ n++; next }} print $2, $5, $7 }} "
           f"END {{ print n+0 > nf }}' | "
           f"LC_ALL=C sort -t $'\\t' -k1,1 -S 40% --parallel={max(1, int(threads))} "
           f"-T {shlex.quote(str(tmpdir))} --compress-program={decompressor()}")
    logger.info(f"  Sorting the GAF annotations by accession (temporary files in {tmpdir}) ...")
    proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=subprocess.PIPE,
                            bufsize=16 * 1024 * 1024)
    return proc, not_count


def build_go_mapping(conn, gaf_gz, flat_tsv=None, threads=8, tmpdir=None):
    """
    go_mapping (accession, go_terms): GO terms of every UniProtKB accession in
    the GOA GAF with the evidence codes behind each term, as
    GO:0005515:IDA,IPI;GO:0005634:IEA. NOT qualifiers are kept out, so a
    term the source states does not apply is not assigned.

    The annotations are sorted by accession on disk first (sorted_gaf_stream),
    then written one accession at a time, so memory use stays small.
    """
    logger.info(f"Parsing {gaf_gz.name} ...")
    tmpdir = Path(tmpdir or gaf_gz.parent) / "gaf_sort_tmp"
    tmpdir.mkdir(exist_ok=True)
    conn.execute("CREATE TABLE go_mapping (accession TEXT, go_terms TEXT)")
    fh = open(flat_tsv, "w") if flat_tsv else None
    batch = []
    rows = 0
    lines = n_not = n_ann = 0
    current = None
    terms = {}

    def flush():
        nonlocal rows, batch
        if current is None:
            return
        text = _format_go_terms(terms)
        batch.append((current, text))
        if fh:
            fh.write(f"{current}\t{text}\n")
        if len(batch) >= BATCH_SIZE:
            conn.executemany("INSERT INTO go_mapping VALUES (?, ?)", batch)
            rows += len(batch)
            batch = []

    proc, not_count = sorted_gaf_stream(gaf_gz, tmpdir, threads)
    first = True
    for raw in proc.stdout:
        if first:
            logger.info("  Sort finished; writing accessions ...")
            first = False
        lines += 1
        if lines % 200_000_000 == 0:
            logger.info(f"  {lines:,} annotations, {rows + len(batch):,} accessions")
        parts = raw.rstrip(b"\n").split(b"\t")
        if len(parts) < 3:
            continue
        acc = parts[0].decode()
        if acc != current:
            flush()
            current = acc
            terms = {}
        terms.setdefault(sys.intern(parts[1].decode()), set()).add(sys.intern(parts[2].strip().decode()))
        n_ann += 1
    flush()
    if batch:
        conn.executemany("INSERT INTO go_mapping VALUES (?, ?)", batch)
        rows += len(batch)
    proc.wait()
    if proc.returncode != 0:
        logger.error(f"  The GAF decompression/filter/sort pipeline failed (exit {proc.returncode})")
        sys.exit(1)
    try:
        n_not = int(not_count.read_text().strip())
        not_count.unlink()
    except (OSError, ValueError):
        pass
    shutil.rmtree(tmpdir, ignore_errors=True)
    if fh:
        fh.close()
        logger.info(f"  Wrote {flat_tsv}")
    logger.info(f"  {n_ann:,} UniProtKB annotations on {rows:,} accessions; "
                f"{n_not:,} NOT annotations left out")

    logger.info("  Indexing go_mapping ...")
    conn.execute("CREATE INDEX idx_go_mapping_acc ON go_mapping(accession)")
    conn.commit()

    # Accessions split over several groups (the GAF is expected to be grouped)
    dups = [r[0] for r in conn.execute(
        "SELECT accession FROM go_mapping GROUP BY accession HAVING COUNT(*) > 1")]
    if dups:
        logger.info(f"  Merging {len(dups):,} accessions that occurred in more than one group ...")
        item_re = re.compile(r"^(GO:\d{7})(?::([A-Za-z,]*))?$")
        for acc in dups:
            merged = {}
            for (text,) in conn.execute("SELECT go_terms FROM go_mapping WHERE accession = ?", (acc,)):
                for item in text.split(";"):
                    m = item_re.match(item)
                    if not m:
                        continue
                    merged.setdefault(m.group(1), set()).update(c for c in (m.group(2) or "").split(",") if c)
            conn.execute("DELETE FROM go_mapping WHERE accession = ?", (acc,))
            conn.execute("INSERT INTO go_mapping VALUES (?, ?)", (acc, _format_go_terms(merged)))
        conn.commit()
        rows = conn.execute("SELECT COUNT(*) FROM go_mapping").fetchone()[0]
    logger.info(f"  go_mapping: {rows:,} accessions")
    return rows


def build_taxonomy(conn, taxdump_dir, constraints_obo, groupings_obo):
    """
    taxonomy_nodes (tax_id, parent_tax_id, rank, name) and taxonomy_merged
    (old_tax_id, new_tax_id) from the NCBI taxdump, for the lineage of the
    query organism; go_taxon_constraints (go_id, relation, taxon) and
    go_taxon_unions (union_id, member) from the GO taxon constraint files.
    """
    names = {}
    with open(taxdump_dir / "names.dmp") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[int(parts[0])] = parts[1]
    conn.execute("CREATE TABLE taxonomy_nodes (tax_id INTEGER PRIMARY KEY, parent_tax_id INTEGER, "
                 "rank TEXT, name TEXT)")
    batch = []
    n_nodes = 0
    with open(taxdump_dir / "nodes.dmp") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            tid = int(parts[0])
            batch.append((tid, int(parts[1]), parts[2], names.get(tid, "")))
            if len(batch) >= BATCH_SIZE:
                conn.executemany("INSERT INTO taxonomy_nodes VALUES (?, ?, ?, ?)", batch)
                n_nodes += len(batch)
                batch = []
    if batch:
        conn.executemany("INSERT INTO taxonomy_nodes VALUES (?, ?, ?, ?)", batch)
        n_nodes += len(batch)
    conn.execute("CREATE TABLE taxonomy_merged (old_tax_id INTEGER PRIMARY KEY, new_tax_id INTEGER)")
    n_merged = 0
    with open(taxdump_dir / "merged.dmp") as fh:
        rows = []
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            rows.append((int(parts[0]), int(parts[1])))
        conn.executemany("INSERT OR IGNORE INTO taxonomy_merged VALUES (?, ?)", rows)
        n_merged = len(rows)
    logger.info(f"  taxonomy_nodes: {n_nodes:,} taxa; taxonomy_merged: {n_merged:,} merged IDs")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from goannotate.taxon import parse_taxon_constraints_obo, parse_taxon_groupings_obo
    constraints = parse_taxon_constraints_obo(constraints_obo)
    unions = parse_taxon_groupings_obo(groupings_obo)
    conn.execute("CREATE TABLE go_taxon_constraints (go_id TEXT, relation TEXT, taxon TEXT)")
    conn.executemany("INSERT INTO go_taxon_constraints VALUES (?, ?, ?)",
                     [(g, rel, tax) for g, rules in constraints.items() for rel, tax in rules])
    conn.execute("CREATE INDEX idx_go_taxon_constraints ON go_taxon_constraints(go_id)")
    conn.execute("CREATE TABLE go_taxon_unions (union_id TEXT, member TEXT)")
    conn.executemany("INSERT INTO go_taxon_unions VALUES (?, ?)",
                     [(u, m) for u, members in unions.items() for m in members])
    conn.commit()
    n_c = sum(len(v) for v in constraints.values())
    logger.info(f"  go_taxon_constraints: {n_c:,} constraints on {len(constraints):,} terms; "
                f"go_taxon_unions: {len(unions):,} union classes")
    return n_nodes, n_c


def replace_tables(main_db, tmp_db, tables):
    """
    Copy the given tables (with their indexes) from tmp_db into main_db,
    replacing the existing ones. Each table is replaced in one transaction.
    """
    conn = sqlite3.connect(str(main_db))
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("ATTACH DATABASE ? AS newdb", (str(tmp_db),))
    for table in tables:
        row = conn.execute("SELECT sql FROM newdb.sqlite_master WHERE type='table' AND name=?",
                           (table,)).fetchone()
        if row is None:
            continue
        indexes = [r[0] for r in conn.execute(
            "SELECT sql FROM newdb.sqlite_master WHERE type='index' AND tbl_name=? AND sql IS NOT NULL",
            (table,))]
        logger.info(f"  Replacing table {table} in {Path(main_db).name} ...")
        # Names are qualified: an unqualified DROP would find the table in
        # the attached database when the main one does not have it yet.
        conn.execute("BEGIN")
        conn.execute(f"DROP TABLE IF EXISTS main.{table}")
        conn.execute(row[0].replace(f"CREATE TABLE {table}", f"CREATE TABLE main.{table}", 1))
        conn.execute(f"INSERT INTO main.{table} SELECT * FROM newdb.{table}")
        for sql in indexes:
            conn.execute(sql.replace("CREATE INDEX ", "CREATE INDEX main.", 1))
        conn.execute("COMMIT")
    conn.execute("DETACH DATABASE newdb")
    conn.close()


def build_idmapping(conn, idmap_gz):
    """
    ncbi_idmapping (uniprot_acc, gene_id): column 3 (GeneID) of UniProt's
    idmapping_selected.tab; an accession with several GeneIDs gets one row each.
    """
    logger.info(f"Building ncbi_idmapping from {idmap_gz.name} ...")
    conn.execute("CREATE TABLE ncbi_idmapping (uniprot_acc TEXT, gene_id TEXT)")
    proc = open_gz_stream(idmap_gz)
    batch = []
    rows = lines = 0
    for raw in proc.stdout:
        lines += 1
        parts = raw.split(b"\t", 3)
        if len(parts) < 3 or not parts[2].strip():
            continue
        acc = parts[0].decode()
        for gid in parts[2].decode().split(";"):
            gid = gid.strip()
            if gid:
                batch.append((acc, gid))
        if len(batch) >= BATCH_SIZE:
            conn.executemany("INSERT INTO ncbi_idmapping VALUES (?, ?)", batch)
            rows += len(batch)
            batch = []
    if batch:
        conn.executemany("INSERT INTO ncbi_idmapping VALUES (?, ?)", batch)
        rows += len(batch)
    close_gz_stream(proc, idmap_gz)
    logger.info("  Indexing ncbi_idmapping ...")
    conn.execute("CREATE INDEX idx_ncbi_idmapping_acc ON ncbi_idmapping(uniprot_acc)")
    # GeneID -> accessions: the GO annotations of relatives' genes (orthology module)
    conn.execute("CREATE INDEX idx_ncbi_idmapping_gene ON ncbi_idmapping(gene_id)")
    conn.commit()
    logger.info(f"  {lines:,} UniProt accessions read; ncbi_idmapping: {rows:,} rows")
    return lines, rows


def build_gene_info(conn, geneinfo_gz):
    """
    ncbi_geneinfo (gene_id, symbol, type): official symbol and NCBI gene type
    (protein-coding, pseudo, ncRNA, ...) of every GeneID.
    ncbi_gene_symbols (tax_id, symbol, gene_id, kind): lowercase names by
    which a gene of a species can be looked up: kind 'symbol' (Symbol and
    Symbol_from_nomenclature_authority), 'locus_tag' (LocusTag, e.g. the
    HispidOSU_009782 symbols of TrEMBL entries from genome projects) and
    'synonym' (Synonyms, e.g. outdated symbols such as AARS for AARS1).
    """
    logger.info(f"Building ncbi_geneinfo and ncbi_gene_symbols from {geneinfo_gz.name} ...")
    conn.execute("CREATE TABLE ncbi_geneinfo (gene_id TEXT, symbol TEXT, type TEXT)")
    conn.execute("CREATE TABLE ncbi_gene_symbols (tax_id TEXT, symbol TEXT, gene_id TEXT, kind TEXT)")
    proc = open_gz_stream(geneinfo_gz)
    info_batch, sym_batch = [], []
    n_genes = 0
    counts = {"symbol": 0, "locus_tag": 0, "synonym": 0}
    for raw in proc.stdout:
        if raw.startswith(b"#"):
            continue
        p = raw.rstrip(b"\n").decode("utf-8", "replace").split("\t")
        if len(p) < 5:
            continue
        tax_id, gene_id, symbol, locus_tag, synonyms = p[0], p[1], p[2], p[3], p[4]
        authority = p[10] if len(p) > 10 else "-"
        gene_type = p[9] if len(p) > 9 else ""
        n_genes += 1
        if symbol and symbol not in ("-", "NEWENTRY"):
            info_batch.append((gene_id, symbol, gene_type))
        names = {}
        for s in (symbol, authority):
            if s and s not in ("-", "NEWENTRY"):
                names.setdefault(s.lower(), "symbol")
        if locus_tag and locus_tag != "-":
            names.setdefault(locus_tag.lower(), "locus_tag")
        if synonyms and synonyms != "-":
            for s in synonyms.split("|"):
                s = s.strip()
                if s and s != "-":
                    names.setdefault(s.lower(), "synonym")
        for name, kind in names.items():
            sym_batch.append((tax_id, name, gene_id, kind))
            counts[kind] += 1
        if len(sym_batch) >= BATCH_SIZE:
            conn.executemany("INSERT INTO ncbi_gene_symbols VALUES (?, ?, ?, ?)", sym_batch)
            sym_batch = []
        if len(info_batch) >= BATCH_SIZE:
            conn.executemany("INSERT INTO ncbi_geneinfo VALUES (?, ?, ?)", info_batch)
            info_batch = []
        if n_genes % 20_000_000 == 0:
            logger.info(f"  {n_genes:,} genes read")
    if sym_batch:
        conn.executemany("INSERT INTO ncbi_gene_symbols VALUES (?, ?, ?, ?)", sym_batch)
    if info_batch:
        conn.executemany("INSERT INTO ncbi_geneinfo VALUES (?, ?, ?)", info_batch)
    close_gz_stream(proc, geneinfo_gz)
    logger.info("  Indexing ...")
    conn.execute("CREATE INDEX idx_ncbi_geneinfo_id ON ncbi_geneinfo(gene_id)")
    conn.execute("CREATE INDEX idx_ncbi_gene_symbols ON ncbi_gene_symbols(tax_id, symbol)")
    conn.commit()
    logger.info(f"  {n_genes:,} genes; ncbi_gene_symbols: {counts['symbol']:,} symbols, "
                f"{counts['locus_tag']:,} locus tags, {counts['synonym']:,} synonyms")
    return n_genes, counts


def build_orthologs(conn, orthologs_gz):
    """ncbi_orthologs: NCBI gene_orthologs (tax_id, GeneID, relationship, Other_tax_id, Other_GeneID)."""
    logger.info(f"Building ncbi_orthologs from {orthologs_gz.name} ...")
    conn.execute("CREATE TABLE ncbi_orthologs "
                 "(tax_id TEXT, gene_id TEXT, relationship TEXT, other_tax_id TEXT, other_gene_id TEXT)")
    batch = []
    rows = 0
    with gzip.open(orthologs_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            batch.append((p[0], p[1], p[2], p[3], p[4]))
            if len(batch) >= BATCH_SIZE:
                conn.executemany("INSERT INTO ncbi_orthologs VALUES (?, ?, ?, ?, ?)", batch)
                rows += len(batch)
                batch = []
    if batch:
        conn.executemany("INSERT INTO ncbi_orthologs VALUES (?, ?, ?, ?, ?)", batch)
        rows += len(batch)
    conn.execute("CREATE INDEX idx_ncbi_orthologs_gene ON ncbi_orthologs(gene_id)")
    conn.execute("CREATE INDEX idx_ncbi_orthologs_other ON ncbi_orthologs(other_gene_id)")
    conn.commit()
    logger.info(f"  ncbi_orthologs: {rows:,} rows")
    return rows


def write_provenance(conn, records, metadata):
    conn.execute("CREATE TABLE IF NOT EXISTS sources (name TEXT, url TEXT, file TEXT, size_bytes INTEGER, "
                 "md5 TEXT, checksum_source TEXT, verified_by TEXT, release TEXT, server_date TEXT, "
                 "downloaded TEXT, verified TEXT)")
    for r in records:
        conn.execute("DELETE FROM sources WHERE name = ?", (r["name"],))
        conn.execute("INSERT INTO sources VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     (r["name"], r["url"], r["file"], r["size"], r.get("md5_expected", ""),
                      r["checksum_source"], r["verified_by"], r["release"], r.get("server_date", ""),
                      r.get("downloaded", ""), r.get("verified", "")))
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    for k, v in metadata.items():
        conn.execute("INSERT OR REPLACE INTO metadata VALUES (?, ?)", (k, str(v)))
    conn.commit()


def write_build_info(path, records, metadata, sqlite_db=None):
    """
    Human-readable record of the build next to the databases. With an
    existing SQLite database, sources and metadata of earlier runs stored
    there are included, so a later --skip-sqlite run keeps them.
    """
    if sqlite_db is not None and Path(sqlite_db).exists():
        conn = sqlite3.connect(str(sqlite_db))
        try:
            stored = {row[0]: row[1] for row in conn.execute("SELECT key, value FROM metadata")}
            stored.update(metadata)
            metadata = stored
            names = {r["name"] for r in records}
            cols = [c[1] for c in conn.execute("PRAGMA table_info(sources)")]
            for row in conn.execute("SELECT * FROM sources"):
                d = dict(zip(cols, row))
                if d["name"] not in names:
                    records = list(records) + [{"name": d["name"], "file": d["file"], "release": d["release"],
                                                "server_date": d["server_date"], "size": d["size_bytes"],
                                                "md5_expected": d["md5"], "checksum_source": d["checksum_source"],
                                                "verified_by": d["verified_by"], "downloaded": d["downloaded"],
                                                "url": d["url"]}]
        except sqlite3.Error:
            pass
        conn.close()
    with open(path, "w") as fh:
        fh.write("# GOAnnotate database build\n")
        for k, v in metadata.items():
            fh.write(f"# {k}: {v}\n")
        cols = ["name", "file", "release", "server_date", "size", "md5_expected", "checksum_source",
                "verified_by", "downloaded", "url"]
        fh.write("\t".join(cols) + "\n")
        for r in records:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


####################
# PROTEIN DATABASES
####################

def extract_taxdump(tar_path, taxdump_dir):
    needed = ["nodes.dmp", "names.dmp", "merged.dmp"]
    if all((taxdump_dir / n).exists() for n in needed) and \
            (not Path(tar_path).exists() or
             min((taxdump_dir / n).stat().st_mtime for n in needed) >= Path(tar_path).stat().st_mtime):
        return
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        for name in needed:
            tar.extract(tar.getmember(name), path=str(taxdump_dir))
    logger.info(f"  Extracted {', '.join(needed)} to {taxdump_dir}")


def resolve_lineage(lineage, names_dmp):
    """(taxon ID, label) for a taxon name or ID; ('', 'all') for every entry."""
    if lineage.lower() == "all":
        return "", "all"
    if lineage.isdigit():
        taxid = lineage
        label = None
        with open(names_dmp) as f:
            for line in f:
                p = line.split("\t|\t")
                if p[0] == taxid and p[3].startswith("scientific name"):
                    label = p[1]
                    break
        if label is None:
            logger.error(f"Taxon ID {taxid} not found in names.dmp")
            sys.exit(1)
        return taxid, re.sub(r"\W+", "_", label.lower()).strip("_")
    wanted = lineage.lower()
    found = set()
    with open(names_dmp) as f:
        for line in f:
            p = line.split("\t|\t")
            if p[1].lower() == wanted and p[3].startswith("scientific name"):
                found.add(p[0])
    if not found:
        logger.error(f"No NCBI taxon with the scientific name '{lineage}'; give a taxon ID instead")
        sys.exit(1)
    if len(found) > 1:
        logger.error(f"'{lineage}' names several taxa ({', '.join(sorted(found))}); give a taxon ID instead")
        sys.exit(1)
    return found.pop(), re.sub(r"\W+", "_", wanted).strip("_")


def load_taxonomy(nodes_path, merged_path, root_taxid):
    children = defaultdict(list)
    all_ids = set()
    with open(nodes_path) as f:
        for line in f:
            p = line.split("\t|\t")
            all_ids.add(p[0])
            if p[0] != p[1]:
                children[p[1]].append(p[0])
    descendants = None
    if root_taxid:
        descendants = {root_taxid}
        stack = [root_taxid]
        while stack:
            for child in children.get(stack.pop(), []):
                descendants.add(child)
                stack.append(child)
        logger.info(f"  {len(descendants):,} of {len(all_ids):,} taxa at or below taxon {root_taxid}")
    merged = {}
    with open(merged_path) as f:
        for line in f:
            p = line.rstrip("\t|\n").split("\t|\t")
            merged[p[0]] = p[1]
    return descendants, all_ids, merged


def filter_lineage(src_gz, out_fasta, out_map, descendants, all_ids, merged, release="", derived=()):
    """
    Entries of src_gz whose OX= taxon is in descendants (all entries when None), with a taxon map.

    release: release of src_gz, recorded in <out_fasta>.release. Existing
    files recorded with another release are rebuilt, together with the
    derived files (Diamond databases) named in derived; existing files
    without a record are kept, with a warning.
    """
    record = out_fasta.with_name(out_fasta.name + ".release")
    if out_fasta.exists() and out_map.exists():
        built = record.read_text().strip() if record.exists() else None
        if built is None:
            logger.warning(f"  [SKIP] {out_fasta.name} exists without a release record; it is kept as it is. "
                           f"Delete it to rebuild it from {release or 'the current download'}.")
            return None
        if not release or built == release:
            logger.info(f"  [SKIP] {out_fasta.name} and {out_map.name} exist ({built or 'release unknown'})")
            return None
        logger.info(f"  {out_fasta.name} was built from {built}; rebuilding it from {release}")
        for old in (out_fasta, out_map, record, *derived):
            if Path(old).exists():
                Path(old).unlink()
    tmp_fasta = out_fasta.with_name(out_fasta.name + ".tmp")
    tmp_map = out_map.with_name(out_map.name + ".tmp")
    logger.info(f"Filtering {src_gz.name} -> {out_fasta.name} ...")
    ox_re = re.compile(rb" OX=(\d+)")
    total = kept = no_ox = unknown = 0
    keep = False
    proc = open_gz_stream(src_gz)
    with open(tmp_fasta, "wb", buffering=16 * 1024 * 1024) as fout, \
            gzip.open(tmp_map, "wt", compresslevel=1) as fmap:
        fmap.write(TAXONMAP_HEADER)
        for line in proc.stdout:
            if line.startswith(b">"):
                total += 1
                if total % 25_000_000 == 0:
                    logger.info(f"  {total:,} entries read, {kept:,} kept")
                keep = False
                m = ox_re.search(line)
                if m is None:
                    no_ox += 1
                else:
                    taxid = m.group(1).decode()
                    taxid = merged.get(taxid, taxid)
                    if taxid not in all_ids:
                        unknown += 1
                    if descendants is None or taxid in descendants:
                        keep = True
                        kept += 1
                        fmap.write(f"{line[1:].split(b'|', 2)[1].decode()}\t{taxid}\n")
            if keep:
                fout.write(line)
    close_gz_stream(proc, src_gz)
    tmp_fasta.rename(out_fasta)
    tmp_map.rename(out_map)
    record.write_text((release or "") + "\n")
    logger.info(f"  {total:,} entries read; {kept:,} kept; {no_ox:,} without OX=; "
                f"{unknown:,} with a taxon ID absent from the NCBI taxonomy")
    return {"read": total, "kept": kept}


# NCBI taxonomy ranks introduced in 2025-2026 that Diamond 2.1.9 and earlier
# reject ("Invalid taxonomic rank: domain"), mapped to ranks it accepts.
# Diamond uses ranks only for rank-based output, not for --taxonlist
# filtering, so the mapping does not change which entries a search returns.
DIAMOND_RANK_MAP = {"domain": "superkingdom", "realm": "superkingdom",
                    "cellular root": "no rank", "acellular root": "no rank"}


def diamond_nodes(nodes_path):
    """A copy of nodes.dmp with ranks Diamond does not know mapped to ones it does (nodes_diamond.dmp)."""
    out = nodes_path.with_name("nodes_diamond.dmp")
    if out.exists() and out.stat().st_mtime >= nodes_path.stat().st_mtime:
        return out
    changed = defaultdict(int)
    tmp = out.with_name(out.name + ".tmp")
    with open(nodes_path) as fin, open(tmp, "w") as fout:
        for line in fin:
            p = line.split("\t|\t")
            if len(p) > 2 and p[2] in DIAMOND_RANK_MAP:
                changed[p[2]] += 1
                p[2] = DIAMOND_RANK_MAP[p[2]]
                line = "\t|\t".join(p)
            fout.write(line)
    tmp.rename(out)
    if changed:
        logger.info("  Ranks mapped for Diamond in nodes_diamond.dmp: " +
                    ", ".join(f"{k} -> {DIAMOND_RANK_MAP[k]} ({v:,})" for k, v in changed.items()))
    return out


def build_diamond(fasta_paths, db_path, taxonmap, nodes, names, threads):
    if db_path.exists():
        logger.info(f"  [SKIP] {db_path.name} exists")
        return
    tmp_prefix = db_path.with_name(db_path.stem + ".tmp")
    tmp_dmnd = tmp_prefix.with_name(tmp_prefix.name + ".dmnd")
    cmd = ["diamond", "makedb", "-d", str(tmp_prefix), "--taxonmap", str(taxonmap),
           "--taxonnodes", str(diamond_nodes(Path(nodes))), "--taxonnames", str(names), "--threads", threads]
    logger.info(f"Building {db_path.name} from {', '.join(p.name for p in fasta_paths)} ...")
    if len(fasta_paths) == 1:
        subprocess.run(cmd + ["--in", str(fasta_paths[0])], check=True)
    else:
        cat = subprocess.Popen(["cat"] + [str(p) for p in fasta_paths], stdout=subprocess.PIPE)
        subprocess.run(cmd, stdin=cat.stdout, check=True)
        cat.stdout.close()
        if cat.wait() != 0:
            logger.error("Concatenation of the FASTA files failed")
            sys.exit(1)
    tmp_dmnd.rename(db_path)
    logger.info(f"  Built {db_path.name} ({format_size(db_path.stat().st_size)})")


def check_diamond_taxonomy(db_path, query_source, workdir, taxid, threads):
    """Search three sequences of the database against it, restricted with --taxonlist, and require taxon IDs."""
    query = workdir / "taxonomy_check_query.fasta"
    n = 0
    with open(query_source) as fin, open(query, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n += 1
                if n > 3:
                    break
            fout.write(line)
    cmd = ["diamond", "blastp", "-q", str(query), "-d", str(db_path), "--outfmt", "6", "qseqid", "sseqid",
           "staxids", "--max-target-seqs", "5", "--threads", threads, "--quiet"]
    if taxid:
        cmd += ["--taxonlist", taxid]
    r = subprocess.run(cmd, capture_output=True, text=True)
    query.unlink()
    rows = [line.split("\t") for line in r.stdout.splitlines() if line.strip()]
    if r.returncode != 0 or not rows or any(len(x) < 3 or x[2] in ("", "0") for x in rows):
        logger.error(f"Taxonomy check of {db_path.name} failed: {r.stderr.strip()[:300]}")
        sys.exit(1)
    logger.info(f"  Taxonomy check of {db_path.name} passed ({len(rows)} hits, all with taxon IDs)")


####################
# MAIN
####################

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("-o", "--output-dir", required=True, help="Directory for all files")
    ap.add_argument("--protein-dbs", metavar="LINEAGE",
                    help="Also download UniProtKB Swiss-Prot and TrEMBL and build Diamond databases of the\n"
                         "entries in LINEAGE: an NCBI taxon name (metazoa) or taxon ID (33208), or 'all'.\n"
                         "Several lineages, comma-separated, build one pair of databases each\n"
                         "(metazoa,viridiplantae) from a single download")
    ap.add_argument("--combined-db", action="store_true",
                    help="With --protein-dbs, also build <lineage>_UniProtKB.dmnd (Swiss-Prot and TrEMBL together)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 8, help="Threads for Diamond (default: all)")
    ap.add_argument("--keep-downloads", action="store_true",
                    help="Keep the downloaded archives after a successful build (default: delete them)")
    ap.add_argument("--flat-files", action="store_true",
                    help="Also write the flat files used by GOAnnotate without --db: Uniprot_All_GO_mapping.tsv,\n"
                         "idmapping_selected.tab and gene_info.tsv")
    ap.add_argument("--skip-sqlite", action="store_true",
                    help="Do not build the SQLite database (protein databases only)")
    ap.add_argument("--tables", metavar="GROUPS", default="all",
                    help="SQLite table groups to build, comma-separated: go, idmapping, gene_info, orthologs,\n"
                         "taxonomy, or all (default: all)")
    ap.add_argument("--update", action="store_true",
                    help="Replace the --tables groups in an existing GOAnnotate_db.sqlite instead of refusing\n"
                         "to touch it; the new tables are built apart and copied in, one table per transaction")
    ap.add_argument("--yes", "-y", action="store_true", help="Do not ask for confirmation")
    args = ap.parse_args()

    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    threads = str(args.threads)
    sqlite_db = outdir / "GOAnnotate_db.sqlite"

    groups = list(TABLE_GROUPS) if args.tables == "all" else [g.strip() for g in args.tables.split(",") if g.strip()]
    unknown = [g for g in groups if g not in TABLE_GROUPS]
    if unknown:
        logger.error(f"Unknown table group(s): {', '.join(unknown)}; choose from {', '.join(TABLE_GROUPS)}")
        sys.exit(1)
    if args.skip_sqlite:
        groups = []

    if groups and sqlite_db.exists() and not args.update:
        logger.error(f"{sqlite_db} exists. Move or delete it to rebuild (the build never overwrites a database), "
                     f"or pass --update with --tables to replace selected tables in it.")
        sys.exit(1)
    if args.update and not sqlite_db.exists():
        logger.error(f"--update given but {sqlite_db} does not exist")
        sys.exit(1)

    print("\n" + "=" * 72 + "\n  GOAnnotate database builder " + SCRIPT_VERSION + "\n" + "=" * 72)
    print(f"  Output directory : {outdir}")
    print(f"  Available disk   : {disk_free_gb(str(outdir)):.0f} GB")
    print(f"  SQLite tables    : {', '.join(groups) if groups else 'no'}"
          + (" (replacing them in the existing database)" if args.update and groups else ""))
    print(f"  Protein databases: {args.protein_dbs or 'no'}")
    print(f"  Keep downloads   : {args.keep_downloads}")
    print()
    if not args.yes:
        if input("  Proceed? [y/N] ").strip().lower() != "y":
            print("  Aborted.")
            sys.exit(0)

    fh = logging.FileHandler(str(outdir / "build_databases.log"))
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    logging.getLogger().addHandler(fh)
    logger.info(f"build_databases.py {SCRIPT_VERSION}; output {outdir}; "
                f"{disk_free_gb(str(outdir)):.0f} GB free")
    check_dependencies(need_diamond=bool(args.protein_dbs))
    start = time.time()

    step_titles = {"go": "GO mapping", "idmapping": "UniProt to GeneID", "gene_info": "NCBI gene_info",
                   "orthologs": "NCBI gene_orthologs", "taxonomy": "Taxonomy and GO taxon constraints"}
    titles = ["Download and verify"]
    titles += [step_titles[g] for g in groups]
    if groups:
        titles += ["Finish SQLite"]
    lineages = [x.strip() for x in (args.protein_dbs or "").split(",") if x.strip()]
    if lineages:
        titles += ["Taxonomy and lineage"]
        for x in lineages:
            titles += [f"Swiss-Prot database ({x})", f"TrEMBL database ({x})"]
            if args.combined_db:
                titles += [f"Combined database ({x})"]
    steps = Steps(titles)

    # ------------------------------------------------------------------
    steps.start("Download and verify")
    names = []
    for g in groups:
        names += [n for n in SOURCES_FOR_GROUP[g] if n not in names]
    if args.protein_dbs:
        names += [n for n in ("sprot_fasta", "trembl_fasta", "taxdump") if n not in names]
    try:
        reldate = http_get_text(UNIPROT_RELDATE).strip().replace("\n", "; ")
        logger.info(f"UniProt: {reldate}")
    except Exception as e:  # noqa: BLE001
        reldate = ""
        logger.warning(f"Could not read {UNIPROT_RELDATE}: {e}")
    records = {}
    for name in names:
        records[name] = fetch(name, outdir)
    uniprot_releases = {r["release"] for n, r in records.items()
                        if n in ("idmapping", "sprot_fasta", "trembl_fasta") and r["release"]}
    if len(uniprot_releases) > 1:
        logger.warning(f"UniProt files come from different releases: {', '.join(sorted(uniprot_releases))}")

    metadata = {
        "build_script_version": SCRIPT_VERSION,
        "build_started": datetime.datetime.fromtimestamp(start).isoformat(timespec="seconds"),
        "uniprot_release": reldate,
    }

    # ------------------------------------------------------------------
    taxdump_dir = outdir / "taxdump"
    if groups:
        tmp_db = outdir / "GOAnnotate_db.sqlite.tmp"
        if tmp_db.exists():
            tmp_db.unlink()
        conn = sqlite3.connect(str(tmp_db))
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA cache_size=-4000000")
        conn.execute("PRAGMA temp_store=MEMORY")

        if "go" in groups:
            steps.start("GO mapping")
            flat_go = outdir / "Uniprot_All_GO_mapping.tsv" if args.flat_files else None
            metadata["rows_go_mapping"] = build_go_mapping(conn, Path(records["goa_gaf"]["path"]), flat_go,
                                                           threads=args.threads, tmpdir=outdir)
            metadata["go_mapping_evidence_codes"] = "yes"

        if "idmapping" in groups:
            steps.start("UniProt to GeneID")
            n_acc, n_rows = build_idmapping(conn, Path(records["idmapping"]["path"]))
            metadata["idmapping_accessions"] = n_acc
            metadata["rows_ncbi_idmapping"] = n_rows

        if "gene_info" in groups:
            steps.start("NCBI gene_info")
            n_genes, counts = build_gene_info(conn, Path(records["gene_info"]["path"]))
            metadata["gene_info_genes"] = n_genes
            metadata["rows_ncbi_gene_symbols"] = sum(counts.values())
            metadata["ncbi_gene_symbols_kinds"] = ", ".join(f"{k} {v}" for k, v in counts.items())

        if "orthologs" in groups:
            steps.start("NCBI gene_orthologs")
            metadata["rows_ncbi_orthologs"] = build_orthologs(conn, Path(records["gene_orthologs"]["path"]))

        if "taxonomy" in groups:
            steps.start("Taxonomy and GO taxon constraints")
            extract_taxdump(Path(records["taxdump"]["path"]), taxdump_dir)
            n_nodes, n_c = build_taxonomy(conn, taxdump_dir, records["go_taxon_constraints"]["path"],
                                          records["go_taxon_groupings"]["path"])
            metadata["rows_taxonomy_nodes"] = n_nodes
            metadata["rows_go_taxon_constraints"] = n_c

        steps.start("Finish SQLite")
        if args.flat_files:
            for group, name, out in (("idmapping", "idmapping", "idmapping_selected.tab"),
                                     ("gene_info", "gene_info", "gene_info.tsv")):
                if group not in groups:
                    continue
                with open(outdir / out, "wb") as fout:
                    subprocess.run([decompressor(), "-dc", records[name]["path"]], stdout=fout, check=True)
                logger.info(f"  Wrote {out}")
        metadata["build_finished_sqlite"] = now_iso()
        source_names = [n for g in groups for n in SOURCES_FOR_GROUP[g]]
        if args.update:
            conn.commit()
            conn.close()
            tables = [t for g in groups for t in TABLE_GROUPS[g]]
            replace_tables(sqlite_db, tmp_db, tables)
            tmp_db.unlink()
            conn = sqlite3.connect(str(sqlite_db))
            write_provenance(conn, [records[n] for n in source_names], metadata)
            conn.execute("ANALYZE")
            conn.commit()
            conn.close()
            logger.info(f"  {sqlite_db.name}: {format_size(sqlite_db.stat().st_size)} "
                        f"(tables replaced: {', '.join(tables)})")
        else:
            write_provenance(conn, [records[n] for n in source_names], metadata)
            conn.execute("ANALYZE")
            conn.commit()
            conn.close()
            tmp_db.rename(sqlite_db)
            logger.info(f"  {sqlite_db.name}: {format_size(sqlite_db.stat().st_size)}")

    # ------------------------------------------------------------------
    if args.protein_dbs:
        steps.start("Taxonomy and lineage")
        extract_taxdump(Path(records["taxdump"]["path"]), taxdump_dir)
        nodes, names_dmp, merged_dmp = taxdump_dir / "nodes.dmp", taxdump_dir / "names.dmp", taxdump_dir / "merged.dmp"
        resolved = [resolve_lineage(x, names_dmp) for x in lineages]
        for taxid, label in resolved:
            logger.info(f"  Lineage: {label}" + (f" (taxon {taxid})" if taxid else " (every entry)"))
        sp_release, tr_release = records["sprot_fasta"]["release"], records["trembl_fasta"]["release"]
        protein_meta = {"protein_dbs_lineage": "; ".join(f"{label} ({taxid or 'all'})" for taxid, label in resolved),
                        "protein_dbs_release": sp_release}

        for taxid, label in resolved:
            descendants, all_ids, merged = load_taxonomy(nodes, merged_dmp, taxid)
            sp_fasta, tr_fasta = outdir / f"{label}_SwissProt.fasta", outdir / f"{label}_TrEMBL.fasta"
            sp_map, tr_map = outdir / f"{label}_SwissProt_taxonmap.tsv.gz", outdir / f"{label}_TrEMBL_taxonmap.tsv.gz"
            sp_dmnd, tr_dmnd = outdir / f"{label}_SwissProt.dmnd", outdir / f"{label}_TrEMBL.dmnd"
            all_map, all_dmnd = outdir / f"{label}_UniProtKB_taxonmap.tsv.gz", outdir / f"{label}_UniProtKB.dmnd"

            steps.start(f"Swiss-Prot database ({label})")
            c = filter_lineage(Path(records["sprot_fasta"]["path"]), sp_fasta, sp_map, descendants, all_ids, merged,
                               release=sp_release, derived=(sp_dmnd, all_map, all_dmnd))
            if c:
                protein_meta[f"protein_dbs_swissprot_entries_{label}"] = f"{c['kept']} of {c['read']}"
            build_diamond([sp_fasta], sp_dmnd, sp_map, nodes, names_dmp, threads)
            check_diamond_taxonomy(sp_dmnd, sp_fasta, outdir, taxid, threads)

            steps.start(f"TrEMBL database ({label})")
            c = filter_lineage(Path(records["trembl_fasta"]["path"]), tr_fasta, tr_map, descendants, all_ids, merged,
                               release=tr_release, derived=(tr_dmnd, all_map, all_dmnd))
            if c:
                protein_meta[f"protein_dbs_trembl_entries_{label}"] = f"{c['kept']} of {c['read']}"
            build_diamond([tr_fasta], tr_dmnd, tr_map, nodes, names_dmp, threads)
            check_diamond_taxonomy(tr_dmnd, tr_fasta, outdir, taxid, threads)

            if args.combined_db:
                steps.start(f"Combined database ({label})")
                if not all_map.exists():
                    tmp = all_map.with_name(all_map.name + ".tmp")
                    with gzip.open(tmp, "wt", compresslevel=1) as fout:
                        fout.write(TAXONMAP_HEADER)
                        for p in (sp_map, tr_map):
                            with gzip.open(p, "rt") as fin:
                                next(fin)
                                shutil.copyfileobj(fin, fout)
                    tmp.rename(all_map)
                build_diamond([sp_fasta, tr_fasta], all_dmnd, all_map, nodes, names_dmp, threads)
        metadata.update(protein_meta)
        if sqlite_db.exists():
            conn = sqlite3.connect(str(sqlite_db))
            write_provenance(conn, list(records.values()), metadata)
            conn.close()
    steps.finish()

    # ------------------------------------------------------------------
    metadata["build_finished"] = now_iso()
    write_build_info(outdir / "build_info.tsv", list(records.values()), metadata, sqlite_db)
    logger.info(f"Wrote {outdir / 'build_info.tsv'}")

    if not args.keep_downloads:
        for r in records.values():
            p = Path(r["path"])
            if r["name"] in ("taxdump", "go_taxon_constraints", "go_taxon_groupings"):
                continue
            if p.exists():
                p.unlink()
                logger.info(f"  Removed {p.name} (its verification record is kept)")

    logger.info(f"Done in {format_elapsed(time.time() - start)}; {disk_free_gb(str(outdir)):.0f} GB free")


if __name__ == "__main__":
    main()
