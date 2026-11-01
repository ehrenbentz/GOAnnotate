"""Base class and shared helpers of the optional annotation modules.

A module wraps an external tool (SignalP, InterProScan, OrthoFinder) that
adds information to the core annotation without changing the symbol and
product vote. Every module follows the same steps:

  detect()     find the tool and its data; report version or the reason
               it cannot run
  run()        run the tool on the query proteins (results in
               intermediate/{prefix}_{module}/ are reused) and parse them
  integrate()  attach per-gene and per-transcript results and return
               counts for the summary

and exposes its results to the writers through gene_columns(),
transcript_columns(), gene_attributes(), transcript_attributes(),
transcript_features() and summary_lines().
"""

import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# Letters SignalP and InterProScan accept; everything else (stop codons,
# selenocysteine U, pyrrolysine O, ambiguity codes) is replaced by X.
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
_NONSTANDARD = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYX]")


@dataclass
class ToolStatus:
    """Result of a module's tool detection."""
    name: str
    available: bool
    path: str = ""
    version: str = ""
    reason: str = ""

    def describe(self) -> str:
        if self.available:
            v = f" (version {self.version})" if self.version else ""
            return f"{self.name}{v}: {self.path}"
        return f"{self.name}: not available ({self.reason})"


@dataclass
class ModuleContext:
    """What a module needs from the pipeline."""
    query_fasta: Path            # protein FASTA that was searched (transcript IDs as headers)
    search_mode: str             # blastp or blastx
    work_dir: Path               # intermediate/
    reports_dir: Path            # reports/
    prefix: str
    threads: int
    force: bool = False
    universe: object = None      # GeneUniverse
    genes: list = None           # parsed GFF3 genes, or None
    annotations: dict = field(default_factory=dict)


class ModuleError(Exception):
    """A module could not complete; the pipeline continues without it."""


class Module:
    """Base of all optional modules. Subclasses set name and description."""

    name = ""
    description = ""
    needs_proteins = True        # skip for blastx (nucleotide) input

    def __init__(self, args):
        self.args = args
        self.status = None       # ToolStatus after detect()
        self.counts = {}         # summary counts after integrate()
        self.runtime = 0.0       # seconds spent in run()
        self.skipped_reason = "" # set when the module was not run

    # ---- tool detection -------------------------------------------------

    def detect(self) -> ToolStatus:
        raise NotImplementedError

    def find_executable(self, option_value: str, default_name: str) -> str:
        """
        Path of the tool: the --<module>-bin value when given (a path to an
        executable, or a name looked up on PATH), else default_name on PATH.
        Returns "" when nothing is found.
        """
        candidate = option_value or default_name
        if os.sep in candidate or candidate.startswith("~"):
            path = Path(os.path.expanduser(candidate))
            if path.is_file() and os.access(path, os.X_OK):
                return str(path)
            return ""
        found = shutil.which(candidate)
        return found or ""

    @staticmethod
    def run_command(cmd: list, timeout: int = 60) -> tuple:
        """(returncode, stdout+stderr) of a short command; (-1, error) when it cannot start."""
        try:
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, timeout=timeout)
            return proc.returncode, proc.stdout
        except (OSError, subprocess.TimeoutExpired) as e:
            return -1, str(e)

    # ---- running ----------------------------------------------------------

    def run(self, ctx: ModuleContext) -> None:
        raise NotImplementedError

    def integrate(self, ctx: ModuleContext) -> dict:
        return {}

    # ---- outputs ----------------------------------------------------------

    def gene_columns(self) -> list:
        """[(column name, fn(gene_key) -> str)] appended to the gene table."""
        return []

    def transcript_columns(self) -> list:
        """[(column name, fn(transcript_key) -> str)] appended to the transcript table."""
        return []

    def gene_attributes(self, gene_key: str) -> dict:
        """GFF3 attributes added to the gene feature."""
        return {}

    def transcript_attributes(self, transcript_key: str) -> dict:
        """GFF3 attributes added to the transcript feature."""
        return {}

    def transcript_features(self, tx, gene) -> list:
        """
        Extra GFF3 feature lines (9 columns, tab-joined, without newline)
        written after the transcript's own line. tx and gene are the parsed
        GFF3 objects.
        """
        return []

    def summary_lines(self) -> list:
        """Lines of the Modules section of the summary (already indented)."""
        return []

    def report_files(self) -> list:
        """Paths of the per-module reports, for the summary."""
        return []


def clean_sequence(seq: str) -> str:
    """Upper case, trailing stop codons removed, internal stops and non-standard letters replaced by X."""
    s = seq.upper().rstrip("*").replace("*", "X")
    return _NONSTANDARD.sub("X", s)


def physical_cores() -> int:
    """
    Number of physical CPU cores (hardware threads counted once), from
    /proc/cpuinfo on Linux; os.cpu_count() where that is not available.
    """
    try:
        cores = set()
        phys = core = None
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("physical id"):
                    phys = line.split(":", 1)[1].strip()
                elif line.startswith("core id"):
                    core = line.split(":", 1)[1].strip()
                elif not line.strip():
                    if core is not None:
                        cores.add((phys, core))
                    phys = core = None
        if core is not None:
            cores.add((phys, core))
        if cores:
            return len(cores)
    except OSError:
        pass
    return os.cpu_count() or 1


def available_memory_gb():
    """MemAvailable from /proc/meminfo in GB, or None where it cannot be read."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024 ** 2
    except (OSError, ValueError, IndexError):
        pass
    return None


def write_clean_proteins(records: list, path: Path, max_length: int = 0) -> int:
    """
    Write protein records as (id, sequence) with every non-standard letter
    replaced by X, trailing stop codons removed, and sequences truncated to
    max_length residues when max_length is set. Records left empty after
    cleaning are skipped. Returns the number written.
    """
    n = 0
    with open(path, "w") as fh:
        for seq_id, seq in records:
            s = clean_sequence(seq)
            if max_length and len(s) > max_length:
                s = s[:max_length]
            if not s:
                continue
            fh.write(f">{seq_id}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
            n += 1
    return n


def protein_to_genomic(tx, aa_start: int, aa_end: int) -> list:
    """
    Genomic intervals [(start, end)] (1-based, inclusive, ascending) of the
    protein residues aa_start..aa_end (1-based, inclusive) of a transcript,
    mapped through its CDS segments in translation order. The first CDS
    segment's phase is skipped, as in translation. Residues beyond the CDS
    are ignored. Returns [] for a non-coding transcript.
    """
    if not tx.cds or aa_end < aa_start:
        return []
    minus = tx.strand == "-"
    ordered = sorted(tx.cds, key=lambda c: c[0], reverse=minus)
    phase = tx.first_phase
    nt_start = (aa_start - 1) * 3 + phase          # 0-based offset in coding sequence
    nt_end = aa_end * 3 - 1 + phase                # inclusive
    intervals = []
    offset = 0
    for s, e, _ in ordered:
        length = e - s + 1
        seg_lo, seg_hi = offset, offset + length - 1
        lo = max(nt_start, seg_lo)
        hi = min(nt_end, seg_hi)
        if lo <= hi:
            if minus:
                g_end = e - (lo - seg_lo)
                g_start = e - (hi - seg_lo)
            else:
                g_start = s + (lo - seg_lo)
                g_end = s + (hi - seg_lo)
            intervals.append((g_start, g_end))
        offset += length
        if offset > nt_end:
            break
    return sorted(intervals)
