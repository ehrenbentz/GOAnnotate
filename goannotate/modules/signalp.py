"""SignalP 6 module: signal peptides and the secretome.

Runs signalp6 (Teufel et al. 2022) on the query proteins and records, per
transcript, whether a Sec/SPI signal peptide was predicted, its
probability and the cleavage site. SignalP 6 is installed by the user
(academic licence, registration at the DTU download page); the module is
skipped with a note when the executable is not found.

Only the first 70 residues of a protein are used by SignalP, so the
proteins are truncated to SP_INPUT_LENGTH residues before the run.
"""

import logging
import re
import time
from pathlib import Path

from goannotate.modules.base import (Module, ModuleContext, ModuleError, ToolStatus,
                                     protein_to_genomic, write_clean_proteins)
from goannotate.sequences import read_fasta

logger = logging.getLogger(__name__)

SP_INPUT_LENGTH = 200
DOWNLOAD_PAGE = "https://services.healthtech.dtu.dk/services/SignalP-6.0/"

_CS_RE = re.compile(r"CS pos:\s*(\d+)-(\d+)\.?\s*(?:Pr:\s*([0-9.]+))?")


class SignalPResult:
    __slots__ = ("prediction", "probability", "cs_start", "cs_end", "cs_probability")

    def __init__(self, prediction, probability, cs_start=0, cs_end=0, cs_probability=0.0):
        self.prediction = prediction          # SP or OTHER (eukarya mode)
        self.probability = probability        # probability of the predicted class
        self.cs_start = cs_start              # last residue of the signal peptide
        self.cs_end = cs_end                  # first residue of the mature protein
        self.cs_probability = cs_probability

    @property
    def has_sp(self) -> bool:
        return self.prediction != "OTHER"


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
            prob_cols = parts[2:-1] if columns and columns[-1].strip().lower().startswith("cs") else parts[2:]
            # Probability of the predicted class: the column whose name
            # starts with the prediction label, else the maximum of the
            # non-OTHER columns (SP), or the OTHER column.
            probability = 0.0
            if columns:
                names = [c.strip() for c in columns[2:2 + len(prob_cols)]]
                for name, value in zip(names, prob_cols):
                    base = name.split("(")[0].strip()
                    if base == prediction:
                        try:
                            probability = float(value)
                        except ValueError:
                            probability = 0.0
                        break
            if probability == 0.0:
                try:
                    values = [float(v) for v in prob_cols]
                    probability = max(values) if values else 0.0
                except ValueError:
                    probability = 0.0
            cs_start = cs_end = 0
            cs_prob = 0.0
            m = _CS_RE.search(parts[-1]) if len(parts) > 2 else None
            if m:
                cs_start, cs_end = int(m.group(1)), int(m.group(2))
                cs_prob = float(m.group(3)) if m.group(3) else 0.0
            results[seq_id] = SignalPResult(prediction, probability, cs_start, cs_end, cs_prob)
    return results


class SignalPModule(Module):
    name = "signalp"
    description = "SignalP 6: Sec/SPI signal peptides and cleavage sites (secretome)"

    def __init__(self, args):
        super().__init__(args)
        self.results = {}            # transcript key -> SignalPResult
        self.gene_sp = {}            # gene key -> (n isoforms with SP, n coding isoforms tested)
        self.mode = getattr(args, "signalp_mode", "fast")
        self.features = getattr(args, "signalp_features", False)
        self.report = None

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

    # ---- run ------------------------------------------------------------

    def _out_dir(self, ctx: ModuleContext) -> Path:
        return ctx.work_dir / f"{ctx.prefix}_signalp"

    def run(self, ctx: ModuleContext) -> None:
        if self.status is None:
            self.detect()
        if not self.status.available:
            raise ModuleError(self.status.reason)
        out_dir = self._out_dir(ctx)
        out_dir.mkdir(parents=True, exist_ok=True)
        results_file = out_dir / "prediction_results.txt"
        clean_fasta = out_dir / "proteins_sp_input.faa"

        records = [(h.split()[0], s) for h, s in read_fasta(str(ctx.query_fasta))]
        ids = {r[0] for r in records}
        reuse = False
        if results_file.exists() and not ctx.force:
            previous = parse_prediction_results(results_file)
            missing = ids - set(previous)
            if not missing:
                reuse = True
                logger.info(f"  Reusing SignalP results in {results_file} ({len(previous):,} proteins)")
                self.results = {k: v for k, v in previous.items() if k in ids}
            else:
                logger.info(f"  SignalP results in {results_file} lack {len(missing):,} proteins; rerunning")
        if not reuse:
            n = write_clean_proteins(records, clean_fasta, max_length=SP_INPUT_LENGTH)
            threads = max(1, ctx.threads or 1)
            cmd = [self.status.path, "--fastafile", str(clean_fasta), "--output_dir", str(out_dir),
                   "--organism", "eukarya", "--format", "none", "--mode", self.mode,
                   "--bsize", "100", "--torch_num_threads", str(threads),
                   "--write_procs", str(min(threads, 8))]
            logger.info(f"  Running SignalP 6 ({self.mode} mode) on {n:,} proteins "
                        f"(first {SP_INPUT_LENGTH} residues): {' '.join(cmd)}")
            t0 = time.time()
            rc, out = self.run_command(cmd, timeout=7 * 24 * 3600)
            self.runtime = time.time() - t0
            (out_dir / "signalp.log").write_text(out)
            if rc != 0 or not results_file.exists():
                raise ModuleError(f"signalp6 failed (exit {rc}); see {out_dir / 'signalp.log'}")
            self.results = parse_prediction_results(results_file)
            logger.info(f"  SignalP finished in {self.runtime / 60:.1f} min; "
                        f"{len(self.results):,} proteins scored")

    # ---- integration ----------------------------------------------------

    def integrate(self, ctx: ModuleContext) -> dict:
        universe = ctx.universe
        self.gene_sp = {}
        for tx_id, gene_id in universe.transcript_to_gene.items():
            r = self.results.get(tx_id)
            if r is None:
                continue
            n_sp, n_tested = self.gene_sp.get(gene_id, (0, 0))
            self.gene_sp[gene_id] = (n_sp + (1 if r.has_sp else 0), n_tested + 1)
        n_tx_sp = sum(1 for r in self.results.values() if r.has_sp)
        n_gene_sp = sum(1 for n, _ in self.gene_sp.values() if n)
        self.counts = {"proteins_scored": len(self.results), "transcripts_with_sp": n_tx_sp,
                       "genes_with_sp": n_gene_sp, "genes_tested": len(self.gene_sp)}
        self.report = ctx.reports_dir / f"{ctx.prefix}_secretome.tsv"
        self._write_secretome(ctx)
        logger.info(f"  Signal peptides: {n_tx_sp:,} proteins on {n_gene_sp:,} genes; "
                    f"secretome list {self.report}")
        return self.counts

    def _write_secretome(self, ctx: ModuleContext) -> None:
        from goannotate.output import product_for
        rows = []
        for tx_id, r in self.results.items():
            if not r.has_sp:
                continue
            gene_id = ctx.universe.transcript_to_gene.get(tx_id, tx_id)
            a = ctx.annotations.get(gene_id)
            rows.append((gene_id, tx_id, a.consensus_symbol if a else "", product_for(a) if a else "",
                         f"{r.probability:.3f}", f"{r.cs_start}-{r.cs_end}" if r.cs_start else "",
                         f"{r.cs_probability:.3f}" if r.cs_start else ""))
        rows.sort()
        with open(self.report, "w") as fh:
            fh.write("gene_id\ttranscript_id\tsymbol\tproduct\tsp_probability\tcleavage_site\tcs_probability\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")

    # ---- outputs --------------------------------------------------------

    def _gene_value(self, gene_key: str) -> str:
        v = self.gene_sp.get(gene_key)
        if v is None:
            return ""
        n_sp, n_tested = v
        if n_sp == 0:
            return "no"
        return "yes" if n_sp == n_tested else f"yes ({n_sp}/{n_tested} isoforms)"

    def _tx_value(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        if r is None:
            return ""
        return "yes" if r.has_sp else "no"

    def _tx_probability(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        return f"{r.probability:.3f}" if r is not None and r.has_sp else ""

    def _tx_cleavage(self, tx_key: str) -> str:
        r = self.results.get(tx_key)
        return f"{r.cs_start}-{r.cs_end}" if r is not None and r.has_sp and r.cs_start else ""

    def gene_columns(self) -> list:
        return [("signal_peptide", self._gene_value)]

    def transcript_columns(self) -> list:
        return [("signal_peptide", self._tx_value), ("sp_probability", self._tx_probability),
                ("sp_cleavage", self._tx_cleavage)]

    def gene_attributes(self, gene_key: str) -> dict:
        v = self.gene_sp.get(gene_key)
        if v and v[0]:
            return {"goannotate_signal_peptide": "yes"}
        return {}

    def transcript_attributes(self, transcript_key: str) -> dict:
        r = self.results.get(transcript_key)
        if r is None or not r.has_sp:
            return {}
        attrs = {"goannotate_signal_peptide": "yes"}
        if r.cs_start:
            attrs["goannotate_sp_cleavage"] = f"{r.cs_start}-{r.cs_end}"
        return attrs

    def transcript_features(self, tx, gene) -> list:
        """sig_peptide feature lines (INSDC type) under the mRNA, with --signalp-features."""
        if not self.features:
            return []
        r = self.results.get(tx.key)
        if r is None or not r.has_sp or not r.cs_start or not tx.coding:
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
                    f"(probability {r.probability:.2f}, cleavage after residue {r.cs_start})"
            lines.append("\t".join([tx.seqid, "GOAnnotate", "sig_peptide", str(start), str(end),
                                    ".", tx.strand, ".", attrs]))
        return lines

    def summary_lines(self) -> list:
        c = self.counts
        lines = [f"  SignalP 6 ({self.status.version or 'version unknown'}, {self.mode} mode): "
                 f"{c.get('proteins_scored', 0):,} proteins scored"
                 + (f" in {self.runtime / 60:.1f} min" if self.runtime else " (results reused)")]
        lines.append(f"    Signal peptides:            {c.get('transcripts_with_sp', 0):,} proteins on "
                     f"{c.get('genes_with_sp', 0):,} of {c.get('genes_tested', 0):,} genes")
        if self.report:
            lines.append(f"    Secretome list:             reports/{self.report.name}")
        lines.append(f"    GFF3 sig_peptide features:  {'written' if self.features else 'off (--signalp-features)'}")
        return lines

    def report_files(self) -> list:
        return [self.report] if self.report else []
