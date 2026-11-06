"""Diamond searches against SwissProt and TrEMBL."""

import logging
import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path

from goannotate import GOAnnotateError

logger = logging.getLogger(__name__)

# Format 6 columns. qlen and slen precede stitle so that the title,
# which may contain spaces, stays last; parse_blast_results recognizes files
# written without them (the 1.x layout, 13 columns).
OUTPUT_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "stitle"]


def results_have_lengths(path: Path) -> bool:
    """True when a results file carries the qlen and slen columns (or is empty)."""
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            return len(parts) >= 15 and parts[12].isdigit() and parts[13].isdigit()
    return True

SENSITIVITY_FLAGS = {
    "fast": [],
    "mid-sensitive": ["--mid-sensitive"],
    "sensitive": ["--sensitive"],
    "more-sensitive": ["--more-sensitive"],
    "very-sensitive": ["--very-sensitive"],
    "ultra-sensitive": ["--ultra-sensitive"],
}


def find_diamond(diamond_bin: str = None) -> str:
    path = diamond_bin or shutil.which("diamond")
    if path and not Path(path).exists():
        path = shutil.which(path)
    if not path:
        raise GOAnnotateError("Diamond executable not found; install Diamond or pass --diamond-bin")
    try:
        out = subprocess.run([path, "version"], capture_output=True, text=True, check=True)
        logger.info(f"Using {path}: {out.stdout.strip()}")
    except (subprocess.CalledProcessError, OSError) as e:
        raise GOAnnotateError(f"Cannot run Diamond at {path}: {e}") from e
    return path


def run_diamond(diamond_bin: str, mode: str, query: Path, database: Path, output: Path,
                threads: int, evalue: float, max_target_seqs: int, sensitivity: str,
                extra_args: str = "") -> None:
    """Run one diamond blastp/blastx search, writing format 6 with stitle."""
    if mode not in ("blastp", "blastx"):
        raise GOAnnotateError(f"Unknown Diamond mode: {mode}")
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    cmd = [diamond_bin, mode,
           "--db", str(database),
           "--query", str(query),
           "--out", str(tmp_output),
           "--outfmt", "6"] + OUTPUT_FIELDS + [
           "--max-target-seqs", str(max_target_seqs),
           "--evalue", str(evalue),
           "--threads", str(threads)]
    cmd += SENSITIVITY_FLAGS.get(sensitivity, [])
    if extra_args:
        cmd += shlex.split(extra_args)

    logger.info("Running: " + " ".join(shlex.quote(c) for c in cmd))
    start = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.strip().splitlines()[-15:])
        raise GOAnnotateError(f"Diamond {mode} against {database.name} failed "
                              f"(exit {proc.returncode}):\n{tail}")
    os.replace(tmp_output, output)
    elapsed = time.time() - start
    n_lines = sum(1 for _ in open(output))
    logger.info(f"  {database.name}: {n_lines:,} hits in {elapsed/60:.1f} min")


def run_searches(mode: str, query: Path, output_dir: Path, prefix: str,
                 swissprot_db: str, trembl_db: str, threads: int, evalue: float,
                 max_target_seqs: int = 100, sensitivity: str = "fast",
                 extra_args: str = "", diamond_bin: str = None, force: bool = False,
                 need_lengths: bool = False) -> Path:
    """
    Search the query set against SwissProt and TrEMBL separately and
    concatenate the results. Separate searches keep SwissProt hits from
    being crowded out of the --max-target-seqs slots by TrEMBL entries.

    Existing result files in the output directory are reused unless
    force is set, or they lack the qlen/slen columns while need_lengths
    is set (the hit coverage filter needs them). Returns the path of the
    combined results file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    combined = output_dir / f"{prefix}_diamond_results.tsv"
    # The search mode of the existing results, so that blastx results of
    # the transcripts are not reused as blastp results of their proteins
    mode_file = output_dir / f"{prefix}_diamond_results.mode"
    previous_mode = mode_file.read_text().strip() if mode_file.exists() else ""

    def reusable(path: Path) -> bool:
        if not (path.exists() and path.stat().st_size > 0) or force:
            return False
        if previous_mode and previous_mode != mode:
            logger.info(f"  {path.name} holds {previous_mode} results; the {mode} search is run")
            return False
        if need_lengths and not results_have_lengths(path):
            logger.info(f"  {path.name} was written without the qlen/slen columns that the hit coverage "
                        f"filter needs; the search is rerun (pass --hit-coverage off to reuse it)")
            return False
        return True

    if reusable(combined):
        logger.info(f"Reusing existing Diamond results: {combined} (use --force-search to rerun)")
        return combined

    databases = []
    if swissprot_db:
        databases.append(("SwissProt", Path(swissprot_db), output_dir / f"{prefix}_diamond_SwissProt.tsv"))
    if trembl_db:
        databases.append(("TrEMBL", Path(trembl_db), output_dir / f"{prefix}_diamond_TrEMBL.tsv"))
    if not databases:
        raise GOAnnotateError("No Diamond database given; pass --diamond-swissprot and/or --diamond-trembl, "
                              "or supply --blast-results")
    for label, db, _ in databases:
        if not db.exists():
            raise GOAnnotateError(f"{label} Diamond database not found: {db}")

    diamond_bin = find_diamond(diamond_bin)
    logger.info(f"Diamond {mode}: {query} against {', '.join(d[0] for d in databases)}")

    for label, db, out in databases:
        if reusable(out):
            logger.info(f"  Reusing existing {label} results: {out}")
            continue
        run_diamond(diamond_bin, mode, query, db, out, threads, evalue,
                    max_target_seqs, sensitivity, extra_args)

    tmp = combined.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fout:
        for _, _, out in databases:
            with open(out) as fin:
                shutil.copyfileobj(fin, fout)
    os.replace(tmp, combined)
    mode_file.write_text(mode + "\n")
    logger.info(f"Combined Diamond results written to {combined}")
    return combined
