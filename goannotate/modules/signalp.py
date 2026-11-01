"""SignalP 6 module: signal peptides and the secretome.

Runs signalp6 (Teufel et al. 2022) on the query proteins and records, per
transcript, whether a Sec/SPI signal peptide was predicted, its
probability and the cleavage site. SignalP 6 is installed by the user
(academic licence, registration at the DTU download page); the module is
skipped with a note when the executable is not found.

SignalP reads only the first 70 residues of a protein, so:
  - the input is truncated to SP_WINDOW residues (results are identical);
  - each distinct 70-residue N-terminus is predicted once and the result
    copied to every transcript that shares it (isoforms often do);
  - the input is split into chunks run as parallel signalp6 processes with
    one PyTorch thread each, one process per physical core (limited by
    --threads and by available memory). On a 32-core machine this was the
    fastest arrangement measured (15 proteins/s against 9.6 for one process
    with 32 threads), with identical predictions.

Optionally (--signalp-alt-starts), transcripts without a signal peptide at
their annotated start are checked for one at an internal methionine within
the first ALT_START_WINDOW residues that is followed by a hydrophobic core
(a possible N-terminal error of the gene model). The result is a QC flag
only; the main call and the gene model are not changed.
"""

import logging
import math
import re
import subprocess
import time
from pathlib import Path

from goannotate.modules.base import (Module, ModuleContext, ModuleError, ToolStatus, available_memory_gb,
                                     clean_sequence, physical_cores, protein_to_genomic)
from goannotate.sequences import read_fasta

logger = logging.getLogger(__name__)

SP_WINDOW = 70             # residues SignalP 6 reads (signalp/utils.py: x[:70])
ALT_START_WINDOW = 60      # internal methionines considered for the alternative-start check
BATCH_SIZE = 100           # SignalP --bsize; 10 (its default) was 1.8 times slower on CPU
PROGRESS_INTERVAL = 60     # seconds between progress lines in the GOAnnotate log
# Memory per signalp6 process: fast measured at about 2.4 GB; slow needs more
# than 14 GB (SignalP documentation); slow-sequential uses the memory of fast.
MEMORY_PER_PROCESS_GB = {"fast": 2.5, "slow-sequential": 2.5, "slow": 16.0}
DOWNLOAD_PAGE = "https://services.healthtech.dtu.dk/services/SignalP-6.0/"

_CS_RE = re.compile(r"CS pos:\s*(\d+)-(\d+)\.?\s*(?:Pr:\s*([0-9.]+))?")
_PROGRESS_RE = re.compile(r"(\d+)/(\d+) \[")
_HYDROPHOBIC = set("AILMFVWC")


class SignalPResult:
    __slots__ = ("prediction", "probability", "sp_probability", "cs_start", "cs_end", "cs_probability")

    def __init__(self, prediction, probability, sp_probability=0.0, cs_start=0, cs_end=0, cs_probability=0.0):
        self.prediction = prediction          # SP or OTHER (eukarya mode)
        self.probability = probability        # probability of the predicted class
        self.sp_probability = sp_probability  # probability of a signal peptide (all non-OTHER classes)
        self.cs_start = cs_start              # last residue of the signal peptide
        self.cs_end = cs_end                  # first residue of the mature protein
        self.cs_probability = cs_probability

    @property
    def has_sp(self) -> bool:
        """SignalP's own call: any class other than OTHER."""
        return self.prediction != "OTHER"

    def shifted(self, offset: int) -> "SignalPResult":
        """The same result with cleavage positions moved by offset residues."""
        return SignalPResult(self.prediction, self.probability, self.sp_probability,
                             self.cs_start + offset if self.cs_start else 0,
                             self.cs_end + offset if self.cs_end else 0, self.cs_probability)


def parse_prediction_results(path: Path) -> dict:
    """
    Parse prediction_results.txt: {sequence ID: SignalPResult}. The header
    line starting with '# ID' names the columns; probability columns lie
    between Prediction and CS Position. Works for the eukarya layout
    (OTHER, SP(Sec/SPI)) and the 'other' layout (more classes).
    """
    results = {}
    columns = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("# ID"):
                    columns = line[2:].split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            seq_id = parts[0].split()[0]
            prediction = parts[1].strip()
            has_cs_col = bool(columns) and columns[-1].strip().lower().startswith("cs")
            prob_cols = parts[2:-1] if has_cs_col else parts[2:]
            names = [c.strip().split("(")[0].strip() for c in columns[2:2 + len(prob_cols)]] if columns else []
            values = []
            for v in prob_cols:
                try:
                    values.append(float(v))
                except ValueError:
                    values.append(0.0)
            probability = 0.0
            sp_probability = 0.0
            if names and len(names) == len(values):
                for name, value in zip(names, values):
                    if name == prediction:
                        probability = value
                    if name != "OTHER":
                        sp_probability = max(sp_probability, value)
            elif values:
                probability = max(values)
                sp_probability = max(values[1:]) if len(values) > 1 else 0.0
            if probability == 0.0 and values:
                probability = max(values)
            cs_start = cs_end = 0
            cs_prob = 0.0
            m = _CS_RE.search(parts[-1]) if len(parts) > 2 else None
            if m:
                cs_start, cs_end = int(m.group(1)), int(m.group(2))
                cs_prob = float(m.group(3)) if m.group(3) else 0.0
            results[seq_id] = SignalPResult(prediction, probability, sp_probability, cs_start, cs_end, cs_prob)
    return results


CORE_LENGTH = 8           # hydrophobic-core window
CORE_MIN_HYDROPHOBIC = 6  # of CORE_LENGTH residues. Calibrated on the prairie vole: 97.8% of the
                          # signal peptides SignalP calls at annotated starts pass from their own Met
                          # (7 of 8: 82.3%, missing cores with S, T or N such as MFQLLLLSNLLW)


def hydrophobic_core(seq: str) -> bool:
    """
    A signal-peptide-like hydrophobic core: at least CORE_MIN_HYDROPHOBIC of
    CORE_LENGTH consecutive residues hydrophobic (A I L M F V W C) within
    residues 3-30 of seq.
    """
    window = seq[2:30]
    return any(sum(c in _HYDROPHOBIC for c in window[i:i + CORE_LENGTH]) >= CORE_MIN_HYDROPHOBIC
               for i in range(max(0, len(window) - CORE_LENGTH + 1)))


def alternative_starts(seq: str, window: int = ALT_START_WINDOW) -> list:
    """
    0-based positions of every internal methionine within the first window
    residues that is followed by a hydrophobic core. Every one is tested:
    taking only the first lets an earlier Met without a signal peptide hide
    a later one with it (16 vole proteins).
    """
    return [i for i in range(1, min(window, len(seq))) if seq[i] == "M" and hydrophobic_core(seq[i:])]


def alternative_start(seq: str, window: int = ALT_START_WINDOW) -> int:
    """0-based position of the first internal methionine passing the gate, or -1."""
    starts = alternative_starts(seq, window)
    return starts[0] if starts else -1


def plan_processes(threads: int, mode: str, requested: int = 0) -> tuple:
    """
    (number of signalp6 processes, reason). requested > 0 is used as given.
    Otherwise one process per physical core, limited by --threads and by
    available memory at MEMORY_PER_PROCESS_GB per process.
    """
    if requested and requested > 0:
        return requested, "--signalp-processes"
    cores = physical_cores()
    n = min(max(1, threads or cores), cores)
    reason = f"{cores} physical cores" if n == cores else f"--threads {threads}"
    mem = available_memory_gb()
    per = MEMORY_PER_PROCESS_GB.get(mode, 2.5)
    if mem is not None:
        cap = max(1, int(mem // per))
        if cap < n:
            n = cap
            reason = f"{mem:.0f} GB available memory at {per:g} GB per process"
    return n, reason


class SignalPModule(Module):
    name = "signalp"
    description = "SignalP 6: Sec/SPI signal peptides and cleavage sites (secretome)"

    def __init__(self, args):
        super().__init__(args)
        self.results = {}            # transcript key -> SignalPResult at the annotated start
        self.alt = {}                # transcript key -> (Met position 1-based, SignalPResult in protein coordinates)
        self.alt_support = {}        # transcript key -> homology note for the alternative start
        self.gene_sp = {}            # gene key -> (n isoforms with SP, n isoforms scored)
        self.mode = getattr(args, "signalp_mode", "fast")
        self.features = getattr(args, "signalp_features", False)
        self.min_probability = getattr(args, "signalp_min_probability", 0.0) or 0.0
        self.check_alt = getattr(args, "signalp_alt_starts", False)
        self.requested_processes = getattr(args, "signalp_processes", 0) or 0
        self.processes = 0
        self.process_reason = ""
        self.n_distinct = 0
        self.n_alt_candidates = 0
        self.n_alt_distinct = 0
        self.report = None
        self.alt_report = None

    # ---- calls ----------------------------------------------------------

    def is_sp(self, r) -> bool:
        """A signal peptide call: SignalP's own call, at --signalp-min-probability or above."""
        return r is not None and r.has_sp and r.sp_probability >= self.min_probability

    # ---- detection ------------------------------------------------------

    def detect(self) -> ToolStatus:
        given = getattr(self.args, "signalp_bin", None)
        path = self.find_executable(given, "signalp6")
        if not path:
            if given:
                reason = f"--signalp-bin {given} is not an executable file or a command on PATH"
            else:
                reason = (f"signalp6 not on PATH (install the academic package from {DOWNLOAD_PAGE} "
                          f"or give --signalp-bin)")
            self.status = ToolStatus("SignalP 6", False, reason=reason)
            return self.status
        rc, out = self.run_command([path, "--version"], timeout=120)
        version = ""
        if rc == 0:
            # "SignalP 6.0 Signal peptide prediction tool 6.0h": the last version token is the release
            found = re.findall(r"\b(\d+\.\d+[a-z0-9.+]*)", out)
            version = found[-1] if found else (out.strip().splitlines()[0] if out.strip() else "")
        else:
            # Versions before 6.0f have no --version; -h must work
            rc2, out2 = self.run_command([path, "-h"], timeout=120)
            if rc2 != 0:
                self.status = ToolStatus("SignalP 6", False, path=path,
                                         reason=f"{path} does not run: {out2.strip()[:200]}")
                return self.status
        missing = self._missing_model(path)
        if missing:
            self.status = ToolStatus("SignalP 6", False, path=path, version=version, reason=missing)
            return self.status
        self.status = ToolStatus("SignalP 6", True, path=path, version=version)
        return self.status

    # Model files each mode needs, inside the installed package's model_weights/
    _MODEL_FILES = {"fast": ["distilled_model_signalp6.pt"],
                    "slow": ["ensemble_model_signalp6.pt", "sequential_models_signalp6"],
                    "slow-sequential": ["sequential_models_signalp6"]}

    def _missing_model(self, path: str) -> str:
        """
        Reason the selected --signalp-mode cannot run because its model
        weights are not installed, or "" when they are (or the location
        cannot be determined, in which case SignalP reports it at run time).
        The installed package is found through the Python interpreter named
        in the signalp6 script's first line, so the environment need not be
        active.
        """
        try:
            with open(path) as fh:
                first = fh.readline().strip()
        except (OSError, UnicodeDecodeError):
            return ""
        if not first.startswith("#!"):
            return ""
        interpreter = first[2:].split()
        if not interpreter:
            return ""
        rc, out = self.run_command(interpreter + ["-c", "import os, signalp; print(os.path.dirname(signalp.__file__))"],
                                   timeout=120)
        if rc != 0 or not out.strip():
            return ""
        weights = Path(out.strip().splitlines()[-1]) / "model_weights"
        wanted = self._MODEL_FILES.get(self.mode, [])
        if wanted and not any((weights / name).exists() for name in wanted):
            return (f"the {self.mode} model is not installed: {weights} has no {' or '.join(wanted)} "
                    f"(copy the model files of the SignalP download into it)")
        return ""

    # ---- running SignalP ------------------------------------------------

    def _predict(self, sequences: dict, out_dir: Path, label: str, force: bool) -> dict:
        """
        Predict {ID: sequence} (already truncated, distinct) with parallel
        signalp6 processes on chunks of the input; return {ID: SignalPResult}.
        Results in out_dir/prediction_results.txt that cover every ID are
        reused unless force is set.
        """
        out_dir.mkdir(parents=True, exist_ok=True)
        merged = out_dir / "prediction_results.txt"
        if merged.exists() and not force:
            previous = parse_prediction_results(merged)
            if all(i in previous for i in sequences):
                logger.info(f"  Reusing SignalP results for the {label} in {merged} ({len(sequences):,} sequences)")
                return previous
        if not sequences:
            return {}

        ids = list(sequences)
        n_proc = max(1, min(self.processes, len(ids)))
        chunks = [ids[i::n_proc] for i in range(n_proc)]
        procs = []
        logs = []
        for k, chunk in enumerate(chunks):
            cdir = out_dir / f"chunk_{k:03d}"
            cdir.mkdir(exist_ok=True)
            fasta = cdir / "input.faa"
            with open(fasta, "w") as fh:
                for i in chunk:
                    fh.write(f">{i}\n{sequences[i]}\n")
            log = cdir / "signalp.log"
            cmd = [self.status.path, "--fastafile", str(fasta), "--output_dir", str(cdir),
                   "--organism", "eukarya", "--format", "none", "--mode", self.mode,
                   "--bsize", str(BATCH_SIZE), "--torch_num_threads", "1", "--write_procs", "1"]
            fh = open(log, "w")
            procs.append((subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT), fh, cdir, len(chunk)))
            logs.append(log)
        logger.info(f"  Running SignalP 6 ({self.mode} mode) on {len(ids):,} {label}: {n_proc} processes "
                    f"x 1 thread ({self.process_reason}); progress in {out_dir}/chunk_*/signalp.log")

        t0 = time.time()
        last = t0
        while any(p.poll() is None for p, _, _, _ in procs):
            time.sleep(5)
            if time.time() - last >= PROGRESS_INTERVAL:
                last = time.time()
                done = sum(_chunk_progress(log) for log in logs)
                logger.info(f"  SignalP {label}: {done:,} of {len(ids):,} ({100 * done / len(ids):.0f}%) "
                            f"after {(last - t0) / 60:.1f} min")
        failed = []
        for p, fh, cdir, _ in procs:
            fh.close()
            if p.returncode != 0 or not (cdir / "prediction_results.txt").exists():
                failed.append(f"{cdir.name} (exit {p.returncode})")
        elapsed = time.time() - t0
        self.runtime += elapsed
        if failed:
            raise ModuleError(f"signalp6 failed for {len(failed)} of {n_proc} chunks: {', '.join(failed[:5])}; "
                              f"see signalp.log in those directories under {out_dir}")

        results = {}
        header = []
        rows = []
        for _, _, cdir, _ in procs:
            path = cdir / "prediction_results.txt"
            results.update(parse_prediction_results(path))
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        if not header or len(header) < 2:
                            header.append(line)
                    elif line.strip():
                        rows.append(line)
        with open(merged, "w") as fh:
            fh.writelines(header[:2])
            fh.writelines(rows)
        missing = [i for i in ids if i not in results]
        if missing:
            raise ModuleError(f"SignalP returned no result for {len(missing):,} sequences (first: {missing[0]})")
        logger.info(f"  SignalP {label} finished in {elapsed / 60:.1f} min ({len(ids) / max(elapsed, 1):.1f} per s)")
        return results

    def run(self, ctx: ModuleContext) -> None:
        if self.status is None:
            self.detect()
        if not self.status.available:
            raise ModuleError(self.status.reason)
        self.processes, self.process_reason = plan_processes(ctx.threads, self.mode, self.requested_processes)
        base_dir = ctx.work_dir / f"{ctx.prefix}_signalp"
        seqs = {}
        for header, seq in read_fasta(str(ctx.query_fasta)):
            s = clean_sequence(seq)
            if s:
                seqs[header.split()[0]] = s

        # One prediction per distinct 70-residue N-terminus
        representative = {}      # N-terminus -> representative ID
        rep_of = {}              # transcript -> representative ID
        for tx, s in seqs.items():
            head = s[:SP_WINDOW]
            rep = representative.setdefault(head, tx)
            rep_of[tx] = rep
        self.n_distinct = len(representative)
        inputs = {rep: head for head, rep in representative.items()}
        logger.info(f"  {len(seqs):,} proteins, {self.n_distinct:,} distinct {SP_WINDOW}-residue N-termini")
        predicted = self._predict(inputs, base_dir / "annotated_starts", "annotated starts", ctx.force)
        self.results = {tx: predicted[rep] for tx, rep in rep_of.items()}

        if self.check_alt:
            candidates = {}          # transcript -> 0-based positions of the gated Mets
            for tx, s in seqs.items():
                if self.is_sp(self.results.get(tx)):
                    continue
                starts = alternative_starts(s)
                if starts:
                    candidates[tx] = starts
            self.n_alt_candidates = len(candidates)
            alt_rep = {}             # 70-residue start -> representative ID (transcript_M<position>)
            for tx, starts in candidates.items():
                for pos in starts:
                    alt_rep.setdefault(seqs[tx][pos:pos + SP_WINDOW], f"{tx}_M{pos + 1}")
            self.n_alt_distinct = len(alt_rep)
            logger.info(f"  Alternative-start check: {len(candidates):,} proteins without a signal peptide have an "
                        f"internal Met within {ALT_START_WINDOW} residues followed by a hydrophobic core "
                        f"({sum(len(v) for v in candidates.values()):,} Mets, {self.n_alt_distinct:,} distinct starts)")
            alt_pred = self._predict({rep: head for head, rep in alt_rep.items()}, base_dir / "alternative_starts",
                                     "alternative starts", ctx.force)
            self.alt = {}
            for tx, starts in candidates.items():
                for pos in starts:      # the first Met (closest to the annotated start) with a signal peptide
                    r = alt_pred.get(alt_rep[seqs[tx][pos:pos + SP_WINDOW]])
                    if self.is_sp(r):
                        self.alt[tx] = (pos + 1, r.shifted(pos))
                        break

    # ---- integration ----------------------------------------------------

    def integrate(self, ctx: ModuleContext) -> dict:
        universe = ctx.universe
        self.gene_sp = {}
        for tx_id, gene_id in universe.transcript_to_gene.items():
            r = self.results.get(tx_id)
            if r is None:
                continue
            n_sp, n_scored = self.gene_sp.get(gene_id, (0, 0))
            self.gene_sp[gene_id] = (n_sp + (1 if self.is_sp(r) else 0), n_scored + 1)
        n_tx_sp = sum(1 for r in self.results.values() if self.is_sp(r))
        n_below = sum(1 for r in self.results.values() if r.has_sp and not self.is_sp(r))
        n_gene_sp = sum(1 for n, _ in self.gene_sp.values() if n)
        self.counts = {"proteins_scored": len(self.results), "distinct_starts": self.n_distinct,
                       "transcripts_with_sp": n_tx_sp, "below_cutoff": n_below,
                       "genes_with_sp": n_gene_sp, "genes_tested": len(self.gene_sp)}
        self.report = ctx.reports_dir / f"{ctx.prefix}_secretome.tsv"
        self._write_secretome(ctx)
        logger.info(f"  Signal peptides: {n_tx_sp:,} proteins on {n_gene_sp:,} genes; secretome list {self.report}"
                    + (f"; {n_below:,} SignalP calls below --signalp-min-probability {self.min_probability:g}"
                       if self.min_probability else ""))
        if self.check_alt:
            self.alt_support = {tx: _homology_note(tx, pos, ctx) for tx, (pos, _) in self.alt.items()}
            n_supported = sum(1 for v in self.alt_support.values() if v)
            self.counts.update({"alt_candidates": self.n_alt_candidates, "alt_distinct": self.n_alt_distinct,
                                "alt_sp": len(self.alt), "alt_supported": n_supported})
            self.alt_report = ctx.reports_dir / f"{ctx.prefix}_signalp_alternative_starts.tsv"
            self._write_alt_report(ctx)
            logger.info(f"  Alternative starts with a signal peptide: {len(self.alt):,} proteins "
                        f"({n_supported:,} with homologs aligned from that start); {self.alt_report}")
        return self.counts

    def _write_secretome(self, ctx: ModuleContext) -> None:
        from goannotate.output import product_for
        rows = []
        for tx_id, r in self.results.items():
            if not self.is_sp(r):
                continue
            gene_id = ctx.universe.transcript_to_gene.get(tx_id, tx_id)
            a = ctx.annotations.get(gene_id)
            rows.append((gene_id, tx_id, a.consensus_symbol if a else "", product_for(a) if a else "",
                         f"{r.sp_probability:.3f}", f"{r.cs_start}-{r.cs_end}" if r.cs_start else "",
                         f"{r.cs_probability:.3f}" if r.cs_start else ""))
        rows.sort()
        with open(self.report, "w") as fh:
            fh.write("gene_id\ttranscript_id\tsymbol\tproduct\tsp_probability\tcleavage_site\tcs_probability\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")

    def _write_alt_report(self, ctx: ModuleContext) -> None:
        from goannotate.output import product_for
        rows = []
        for tx_id, (pos, r) in self.alt.items():
            gene_id = ctx.universe.transcript_to_gene.get(tx_id, tx_id)
            a = ctx.annotations.get(gene_id)
            annotated = self.results.get(tx_id)
            rows.append((gene_id, tx_id, a.consensus_symbol if a else "", product_for(a) if a else "",
                         f"{annotated.sp_probability:.3f}" if annotated else "", f"M{pos}",
                         f"{r.sp_probability:.3f}", f"{r.cs_start}-{r.cs_end}" if r.cs_start else "",
                         self.alt_support.get(tx_id, "")))
        rows.sort()
        with open(self.alt_report, "w") as fh:
            fh.write("gene_id\ttranscript_id\tsymbol\tproduct\tsp_probability_annotated_start\t"
                     "alternative_start\tsp_probability_alternative_start\tcleavage_site\thomology_support\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")

    # ---- outputs --------------------------------------------------------

    def _gene_value(self, gene_key: str) -> str:
        v = self.gene_sp.get(gene_key)
        if v is None:
            return ""
        n_sp, n_scored = v
        if n_sp == 0:
            return "no"
        return "yes" if n_sp == n_scored else f"yes ({n_sp}/{n_scored} isoforms)"

    def _tx_value(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        if r is None:
            return ""
        return "yes" if self.is_sp(r) else "no"

    def _tx_probability(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        return f"{r.sp_probability:.3f}" if self.is_sp(r) else ""

    def _tx_cleavage(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        return f"{r.cs_start}-{r.cs_end}" if self.is_sp(r) and r.cs_start else ""

    def _tx_alt(self, tx_key: str) -> str:
        v = self.alt.get(tx_key)
        if v is None:
            return ""
        pos, r = v
        text = f"M{pos}: SP {r.sp_probability:.2f}"
        if r.cs_start:
            text += f", cleavage {r.cs_start}-{r.cs_end}"
        note = self.alt_support.get(tx_key)
        if note:
            text += f"; {note}"
        return text

    def gene_columns(self) -> list:
        return [("signal_peptide", self._gene_value)]

    def transcript_columns(self) -> list:
        cols = [("signal_peptide", self._tx_value), ("sp_probability", self._tx_probability),
                ("sp_cleavage", self._tx_cleavage)]
        if self.check_alt:
            cols.append(("sp_alternative_start", self._tx_alt))
        return cols

    def gene_attributes(self, gene_key: str) -> dict:
        v = self.gene_sp.get(gene_key)
        if v and v[0]:
            return {"goannotate_signal_peptide": "yes"}
        return {}

    def transcript_attributes(self, transcript_key: str) -> dict:
        attrs = {}
        r = self.results.get(transcript_key)
        if self.is_sp(r):
            attrs["goannotate_signal_peptide"] = "yes"
            if r.cs_start:
                attrs["goannotate_sp_cleavage"] = f"{r.cs_start}-{r.cs_end}"
        alt = self.alt.get(transcript_key)
        if alt:
            attrs["goannotate_sp_alternative_start"] = f"M{alt[0]}"
        return attrs

    def transcript_features(self, tx, gene) -> list:
        """sig_peptide feature lines (INSDC type) under the mRNA, with --signalp-features."""
        if not self.features:
            return []
        r = self.results.get(tx.key)
        if not self.is_sp(r) or not r.cs_start or not tx.coding:
            return []
        intervals = protein_to_genomic(tx, 1, r.cs_start)
        if not intervals:
            return []
        from goannotate.gff3 import encode_attribute_value
        fid = encode_attribute_value(f"sig_peptide-{tx.id}")
        parent = encode_attribute_value(tx.id)
        lines = []
        for start, end in intervals:
            attrs = f"ID={fid};Parent={parent};Note=signal peptide predicted by SignalP 6 " \
                    f"(probability {r.sp_probability:.2f}, cleavage after residue {r.cs_start})"
            lines.append("\t".join([tx.seqid, "GOAnnotate", "sig_peptide", str(start), str(end),
                                    ".", tx.strand, ".", attrs]))
        return lines

    def summary_lines(self) -> list:
        c = self.counts
        lines = [f"  SignalP 6 ({self.status.version or 'version unknown'}, {self.mode} mode): "
                 f"{c.get('proteins_scored', 0):,} proteins, {c.get('distinct_starts', 0):,} distinct "
                 f"{SP_WINDOW}-residue N-termini"
                 + (f", {self.runtime / 60:.1f} min on {self.processes} processes" if self.runtime
                    else " (results reused)")]
        cutoff = (f"SignalP's call at probability >= {self.min_probability:g}" if self.min_probability
                  else "SignalP's own call")
        lines.append(f"    Signal peptides:            {c.get('transcripts_with_sp', 0):,} proteins on "
                     f"{c.get('genes_with_sp', 0):,} of {c.get('genes_tested', 0):,} genes ({cutoff})")
        if self.min_probability:
            lines.append(f"    Below the cutoff:           {c.get('below_cutoff', 0):,} SignalP calls")
        if self.report:
            lines.append(f"    Secretome list:             reports/{self.report.name}")
        lines.append(f"    GFF3 sig_peptide features:  {'written' if self.features else 'off (--signalp-features)'}")
        if self.check_alt:
            lines.append(f"    Alternative starts:         {c.get('alt_candidates', 0):,} candidates, "
                         f"{c.get('alt_sp', 0):,} with a signal peptide ({c.get('alt_supported', 0):,} with "
                         f"homologs aligned from that start); reports/{self.alt_report.name if self.alt_report else ''}")
        return lines

    def report_files(self) -> list:
        return [p for p in (self.report, self.alt_report) if p]


def _chunk_progress(log: Path) -> int:
    """Sequences done in one chunk, from the last tqdm counter in its log."""
    try:
        with open(log, "rb") as fh:
            fh.seek(0, 2)
            size = fh.tell()
            fh.seek(max(0, size - 4000))
            tail = fh.read().decode("utf-8", "replace")
    except OSError:
        return 0
    found = _PROGRESS_RE.findall(tail)
    return int(found[-1][0]) if found else 0


def _homology_note(tx: str, met_pos: int, ctx: ModuleContext) -> str:
    """
    Homology support for an alternative start at residue met_pos (1-based):
    the best hit of the transcript that aligns from the database protein's
    first residues (sstart <= 3) begins its alignment at or after the
    internal methionine (qstart >= met_pos - 2), so the residues before it
    have no counterpart in the homolog. Returns "" without support.
    """
    gene = ctx.universe.transcript_to_gene.get(tx, tx)
    annot = ctx.annotations.get(gene)
    if annot is None or not annot.hits:
        return ""
    best = None
    for h in annot.hits:
        if h.query_id != tx or min(h.sstart, h.send) > 3:
            continue
        if best is None or h.bitscore > best.bitscore:
            best = h
    if best is None:
        return ""
    qstart = min(best.qstart, best.qend)
    if qstart >= met_pos - 2:
        return f"homolog {best.accession} aligns from residue {qstart}"
    return ""
