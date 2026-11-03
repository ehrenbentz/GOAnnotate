"""InterProScan 6 module: protein domains, InterPro GO terms, domain-based product names.

Runs the InterProScan 6 Nextflow workflow (EMBL-EBI; Blum et al. 2025,
Nucleic Acids Research) on the query proteins and adds, per gene and
transcript, the InterPro entries and Pfam families found, the GO terms of
the InterPro entries (InterPro2GO), and, for coding genes that the core
annotation left without a product, a domain-based name such as
"uncharacterized EGF-like domain-containing protein".

InterProScan 6 is a Nextflow workflow run in containers: it needs Nextflow,
a container runtime (Docker, Singularity, Apptainer or Podman) and the
InterPro data for one release, downloaded beforehand. The module never lets
InterProScan download data during a run: InterProScan 6.0.3.0 downloads
member databases that share a directory (PROSITE patterns and profiles,
CATH-Gene3D and CATH-FunFam) with a race that can start a search on an
empty directory.
"""

import json
import logging
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

from goannotate.modules.base import Module, ModuleContext, ModuleError, ToolStatus, clean_sequence
from goannotate.sequences import read_fasta

logger = logging.getLogger(__name__)

DEFAULT_WORKFLOW = "ebi-pf-team/interproscan6"
PROGRESS_INTERVAL = 60
OUTPREFIX = "interproscan"
RESOURCES_CONFIG = "goannotate_resources.config"

# CPUs requested by the InterProScan steps that run several threads although the workflow
# requests one CPU per task (MobiDB-lite runs several disorder predictors at once).
# Nextflow counts these CPUs against --max-workers, so the parallel tasks fit the threads given.
STEP_CPUS = {"RUN_MOBIDBLITE": 4}
# Container runtimes that take --cpus: each task's container is held to the CPUs it requested
CPU_LIMIT_RUNTIMES = ("docker", "podman")


def resources_config(profile: str) -> str:
    """Nextflow configuration with the CPUs of the multi-threaded steps and, for Docker and
    Podman, a CPU limit on each task's container equal to the CPUs it requested."""
    lines = ["// Written by GOAnnotate: CPUs of the InterProScan steps", "process {"]
    for step, cpus in STEP_CPUS.items():
        lines.append(f"    withName: '{step}' {{ cpus = {cpus} }}")
    if profile in CPU_LIMIT_RUNTIMES:
        lines.append('    containerOptions = { "--cpus ${task.cpus}" }')
    lines.append("}")
    return "\n".join(lines) + "\n"
DOMAIN_TYPES = ("Domain", "Repeat")       # InterPro entry types used for product names
MAX_NAMED_DOMAINS = 2


def _normalise(name: str) -> str:
    """Application name as InterProScan compares them: no spaces or hyphens, lower case."""
    return re.sub(r"[\s\-]+", "", name).lower()


def parse_applications_config(path: Path) -> dict:
    """
    {application: directory} for the applications with data
    (has_data = true) in InterProScan 6's conf/applications.config.
    """
    apps = {}
    current = None
    block = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^\s{8}(\w+)\s*\{\s*$", line)
            if m:
                current = m.group(1)
                block = {}
                continue
            if current is None:
                continue
            m = re.match(r'^\s+(\w+)\s*=\s*"?([^"\n]*)"?\s*$', line)
            if m:
                block[m.group(1)] = m.group(2).strip()
            if re.match(r"^\s{8}\}\s*$", line):
                if block.get("has_data") == "true" and block.get("dir"):
                    apps[current] = block["dir"]
                current = None
    return apps


# Entries named after an uncharacterized protein ("UPF0606 protein KIAA1549, alpha/beta domain 1",
# "C2orf74 domain"): without the identifier only a generic fold remains, so they give no name.
_IDENTIFIER_ENTRY = re.compile(r"\b(?:UPF\d{4}|KIAA\d{3,4}L?|C(?:\d{1,2}|X|Y)orf\d+[A-Z]?)\b", re.I)


def domain_label(description: str, entry_type: str) -> str:
    """
    'EGF-like domain' -> 'EGF-like domain'; 'Protein kinase domain' ->
    'protein kinase domain'; 'WD40 repeat' -> 'WD40 repeat'. A capitalized
    ordinary first word is put in lower case (acronyms and mixed-case names
    keep their case); commas are removed; the entry type's word is appended
    when the description does not end with it. Entries named after an
    uncharacterized protein (UPF, KIAA and chromosome-ORF identifiers) give ''.
    """
    if _IDENTIFIER_ENTRY.search(description):
        return ""
    suffix = "domain" if entry_type == "Domain" else "repeat"
    m = re.match(r"(?i)^domain of unknown function,?\s*(DUF\d+)", description.strip())
    if m:                                   # "Domain of unknown function DUF4515" -> "DUF4515 domain"
        return f"{m.group(1).upper()} domain"
    text = description.strip().replace(",", "")
    # "ABC transporter type 1, transmembrane domain MetI-like" -> "... transmembrane domain"
    m = re.search(rf"\b{suffix}\b", text, re.I)
    if m:
        text = text[:m.end()]
    words = text.split()
    if words and re.fullmatch(r"[A-Z][a-z]+", words[0]):
        words[0] = words[0].lower()
    text = " ".join(words)
    if not text.lower().endswith(suffix):
        text = f"{text} {suffix}"
    return text


def domain_product(entries: list) -> str:
    """
    Product name from (start, description, type) of the InterPro domain and
    repeat entries of a protein, in order of position:
    'uncharacterized EGF-like domain-containing protein', with two entries
    'uncharacterized X and Y domain-containing protein'. Domains come before
    repeats; at most MAX_NAMED_DOMAINS distinct names.
    """
    labels = []
    for kind in DOMAIN_TYPES:
        for _, desc, etype in sorted(e for e in entries if e[2] == kind):
            label = domain_label(desc, etype)
            if label and label not in labels:
                labels.append(label)
    labels = labels[:MAX_NAMED_DOMAINS]
    if not labels:
        return ""
    if len(labels) == 1:
        return f"uncharacterized {labels[0]}-containing protein"
    first, second = labels
    # "X domain and Y domain" -> "X and Y domain" when both are domains (or both repeats)
    for word in ("domain", "repeat"):
        if first.endswith(" " + word) and second.endswith(" " + word):
            return f"uncharacterized {first[:-len(word) - 1]} and {second}-containing protein"
    return f"uncharacterized {first} and {second}-containing protein"


def parse_tsv(path: Path) -> dict:
    """
    InterProScan TSV -> {protein ID: [match dict]}. Columns: protein,
    md5, length, analysis, signature accession, signature description,
    start, end, score, status, date, InterPro accession, InterPro
    description, GO terms, pathways (the last four may be '-' or absent).
    """
    matches = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ipr = parts[11] if len(parts) > 11 and parts[11] not in ("-", "") else ""
            go = []
            if len(parts) > 13 and parts[13] not in ("-", ""):
                go = [g.split("(")[0] for g in parts[13].split("|") if g.startswith("GO:")]
            try:
                start, end = int(parts[6]), int(parts[7])
            except ValueError:
                start = end = 0
            matches.setdefault(parts[0], []).append({
                "analysis": parts[3], "signature": parts[4],
                "signature_description": parts[5] if parts[5] != "-" else "",
                "start": start, "end": end, "score": parts[8],
                "interpro": ipr, "interpro_description": parts[12] if len(parts) > 12 and parts[12] != "-" else "",
                "go": go})
    return matches


class InterProScanModule(Module):
    name = "interproscan"
    description = "InterProScan 6: InterPro entries, Pfam families, InterPro2GO terms, domain-based product names"

    def __init__(self, args):
        super().__init__(args)
        self.nextflow = getattr(args, "interproscan_nextflow", None)
        self.workflow = getattr(args, "interproscan_workflow", None)
        self.revision = getattr(args, "interproscan_revision", None)
        self.datadir = getattr(args, "interproscan_datadir", None)
        self.interpro = getattr(args, "interproscan_interpro", None)
        self.profile = getattr(args, "interproscan_profile", "auto") or "auto"
        self.apps = getattr(args, "interproscan_apps", None)
        self.local = getattr(args, "interproscan_local", False)
        self.nxf_home = getattr(args, "nxf_home", None)
        self.go_mode = getattr(args, "interpro_go", "add")
        self.product_mode = getattr(args, "interpro_products", "fill")
        self.keep_work = getattr(args, "interproscan_keep_work", False)
        self.workers = getattr(args, "interproscan_workers", 0) or 0
        self.ips_version = ""
        self.matches = {}            # transcript -> [match]
        self.entries = {}            # InterPro accession -> (description, type)
        self.gene_entries = {}       # gene -> [(first start, accession)]
        self.tx_entries = {}         # transcript -> [(first start, accession)]
        self.gene_pfam = {}
        self.tx_pfam = {}
        self.gene_go_added = {}      # gene -> sorted specific terms that came from InterPro only
        self.report = None

    # ---- environment ------------------------------------------------------

    def _env(self) -> dict:
        env = dict(os.environ)
        if self.nxf_home:
            env["NXF_HOME"] = str(Path(self.nxf_home).expanduser())
        return env

    def _local_workflow(self):
        if self.workflow and Path(self.workflow).expanduser().is_dir():
            return Path(self.workflow).expanduser()
        return None

    def _profile(self) -> str:
        if self.profile != "auto":
            return self.profile
        if shutil.which("docker"):
            rc, _ = self.run_command(["docker", "info"], timeout=60)
            if rc == 0:
                return "docker"
        for runtime in ("singularity", "apptainer", "podman"):
            if shutil.which(runtime):
                return runtime
        return ""

    def _interpro_version(self, datadir: Path) -> str:
        if self.interpro:
            return self.interpro
        versions = sorted(p.name for p in (datadir / "interpro").glob("*") if (p / "databases.json").exists())
        return versions[-1] if len(versions) == 1 else ""

    # ---- detection ----------------------------------------------------------

    def detect(self) -> ToolStatus:
        name = "InterProScan 6"
        nxf = self.find_executable(self.nextflow, "nextflow")
        if not nxf:
            self.status = ToolStatus(name, False, reason="nextflow not found (give --interproscan-nextflow)")
            return self.status
        if not self.workflow:
            self.status = ToolStatus(name, False, path=nxf,
                                     reason="give --interproscan-workflow (a local clone of ebi-pf-team/interproscan6 "
                                            "at a release tag, or the GitHub name with --interproscan-revision)")
            return self.status
        local = self._local_workflow()
        if local is None and not self.revision:
            self.status = ToolStatus(name, False, path=nxf, reason=f"{self.workflow} is not a local directory; "
                                     f"a remote workflow needs --interproscan-revision (for example 6.0.3.0)")
            return self.status
        if not self.datadir or not Path(self.datadir).is_dir():
            self.status = ToolStatus(name, False, path=nxf,
                                     reason="give --interproscan-datadir, the directory with the InterPro data")
            return self.status
        datadir = Path(self.datadir)
        version = self._interpro_version(datadir)
        if not version:
            found = sorted(p.name for p in (datadir / "interpro").glob("*"))
            self.status = ToolStatus(name, False, path=nxf, reason=(
                f"no InterPro release in {datadir}/interpro" if not found else
                f"several InterPro releases in {datadir}/interpro ({', '.join(found)}); give --interproscan-interpro"))
            return self.status
        self.interpro = version
        ipr_dir = datadir / "interpro" / version
        for f in ("databases.json", "entries.json", "goterms.ipr.json"):
            if not (ipr_dir / f).exists():
                self.status = ToolStatus(name, False, path=nxf, reason=f"{ipr_dir / f} is missing")
                return self.status
        if local is not None:
            missing = self._missing_data(local, datadir, ipr_dir)
            if missing:
                self.status = ToolStatus(name, False, path=nxf, reason=(
                    f"InterPro {version} data missing in {datadir} for: {', '.join(missing)} "
                    f"(download it beforehand; see the README)"))
                return self.status
            self.ips_version = _workflow_version(local)
        else:
            self.ips_version = self.revision
        if not self._profile():
            self.status = ToolStatus(name, False, path=nxf,
                                     reason="no container runtime found (Docker, Singularity, Apptainer or Podman)")
            return self.status
        rc, out = self.run_command([nxf, "-version"], timeout=300)
        if rc != 0:
            self.status = ToolStatus(name, False, path=nxf, reason=f"nextflow does not run: {out.strip()[:200]}")
            return self.status
        m = re.search(r"version\s+(\d+\.\d+\.\d+)", out)
        nf_version = m.group(1) if m else ""
        self.status = ToolStatus(name, True, path=nxf, version=(
            f"{self.ips_version or 'unknown'}, InterPro {version}, Nextflow {nf_version or 'unknown'}, "
            f"profile {self._profile()}"))
        return self.status

    def _selected_apps(self, available: dict) -> list:
        if not self.apps:
            return list(available)
        wanted = [_normalise(a) for a in self.apps.replace(" ", ",").split(",") if a.strip()]
        return [a for a in available if _normalise(a) in wanted]

    def _missing_data(self, workflow: Path, datadir: Path, ipr_dir: Path) -> list:
        """Applications whose data directory for this InterPro release is missing or empty."""
        config = workflow / "conf" / "applications.config"
        if not config.exists():
            return []
        apps = parse_applications_config(config)
        versions = {_normalise(k): v for k, v in json.load(open(ipr_dir / "databases.json")).items()}
        missing = []
        for app in self._selected_apps(apps):
            version = versions.get(_normalise(app))
            if version is None:
                continue
            first, _, sub = apps[app].partition("/")
            path = datadir / first / version / sub
            if not path.is_dir() or not any(path.iterdir()):
                missing.append(f"{app} {version}")
        return missing

    # ---- running --------------------------------------------------------------

    def run(self, ctx: ModuleContext) -> None:
        if self.status is None:
            self.detect()
        if not self.status.available:
            raise ModuleError(self.status.reason)
        # Absolute paths: Nextflow runs with out_dir as its working directory
        out_dir = (ctx.work_dir / f"{ctx.prefix}_interproscan").resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        tsv = out_dir / f"{OUTPREFIX}.tsv"
        marker = out_dir / f"{OUTPREFIX}.done"

        fasta = out_dir / "proteins.faa"
        n = 0
        with open(fasta, "w") as fh:
            for header, seq in read_fasta(str(ctx.query_fasta)):
                s = clean_sequence(seq)
                if s:
                    fh.write(f">{header.split()[0]}\n{s}\n")
                    n += 1
        settings = {"input_size": fasta.stat().st_size, "proteins": n, "interpro": self.interpro,
                    "workflow": str(self.workflow), "revision": self.revision or "", "apps": self.apps or "all",
                    "local": bool(self.local)}
        if tsv.exists() and marker.exists() and not ctx.force:
            try:
                if json.loads(marker.read_text()) == settings:
                    logger.info(f"  Reusing InterProScan results in {tsv}")
                    self._load(tsv)
                    return
            except (json.JSONDecodeError, OSError):
                pass

        # CPUs: --max-workers is the CPU budget of Nextflow's local executor (--threads by default).
        # Most tasks request one CPU; the multi-threaded steps request their real use (STEP_CPUS)
        # and, with Docker or Podman, each container is limited to the CPUs it requested, so the
        # load stays within the budget. The configuration is added with Nextflow's -c option.
        profile = self._profile()
        workers = self.workers or max(1, ctx.threads or 2)
        resources = out_dir / RESOURCES_CONFIG
        resources.write_text(resources_config(profile))
        local = self._local_workflow()
        cmd = [self.status.path, "-c", str(resources), "run", str(local) if local else self.workflow]
        if not local:
            cmd += ["-r", self.revision]
        cmd += ["-profile", profile, "-ansi-log", "false", "-w", str(out_dir / "work"),
                "-resume",            # an interrupted run reuses its finished tasks (same work directory)
                "--input", str(fasta), "--datadir", str(Path(self.datadir).resolve()), "--interpro", self.interpro,
                "--goterms", "--formats", "tsv", "--outdir", str(out_dir), "--outprefix", OUTPREFIX,
                "--max-workers", str(workers)]
        if self.apps:
            cmd += ["--applications", self.apps.replace(" ", "")]
        if self.local:
            cmd += ["--no-matches-api", "--skip-interpro-version-check"]
        log = out_dir / "nextflow.out"
        logger.info(f"  Running InterProScan on {n:,} proteins"
                    f"{' (all analyses run locally)' if self.local else ' (precalculated matches looked up at EBI by checksum)'}: "
                    f"{' '.join(cmd)}")
        t0 = time.time()
        last = t0
        with open(log, "w") as fh:
            proc = subprocess.Popen(cmd, cwd=str(out_dir), stdout=fh, stderr=subprocess.STDOUT, env=self._env())
            while proc.poll() is None:
                time.sleep(5)
                if time.time() - last >= PROGRESS_INTERVAL:
                    last = time.time()
                    submitted, done = _nextflow_progress(out_dir / ".nextflow.log")
                    logger.info(f"  InterProScan: {done} of {submitted} submitted tasks done after "
                                f"{(last - t0) / 60:.1f} min")
        self.runtime = time.time() - t0
        if proc.returncode != 0 or not tsv.exists():
            raise ModuleError(f"InterProScan failed (exit {proc.returncode}); see {log} and {out_dir}/.nextflow.log")
        marker.write_text(json.dumps(settings))
        logger.info(f"  InterProScan finished in {self.runtime / 60:.1f} min")
        if not self.keep_work:
            shutil.rmtree(out_dir / "work", ignore_errors=True)
        self._load(tsv)

    def _load(self, tsv: Path) -> None:
        self.matches = parse_tsv(tsv)
        ipr_dir = Path(self.datadir) / "interpro" / self.interpro
        needed = {m["interpro"] for ms in self.matches.values() for m in ms if m["interpro"]}
        with open(ipr_dir / "entries.json") as fh:
            all_entries = json.load(fh)
        self.entries = {acc: (all_entries[acc].get("description") or all_entries[acc].get("name") or "",
                              all_entries[acc].get("type") or "")
                        for acc in needed if acc in all_entries}
        del all_entries
        with open(ipr_dir / "goterms.ipr.json") as fh:
            ipr_go = json.load(fh)
        self.ipr_go = {acc: ipr_go.get(acc, []) for acc in needed}

    # ---- integration ------------------------------------------------------------

    def integrate(self, ctx: ModuleContext) -> dict:
        universe = ctx.universe
        seq_len = {}
        for header, seq in read_fasta(str(ctx.query_fasta)):
            seq_len[header.split()[0]] = len(seq.rstrip("*"))

        # Entries and Pfam families per transcript and gene
        for tx, ms in self.matches.items():
            first = {}
            for m in ms:
                if m["interpro"]:
                    first[m["interpro"]] = min(first.get(m["interpro"], m["start"]), m["start"])
            self.tx_entries[tx] = sorted((s, acc) for acc, s in first.items())
            self.tx_pfam[tx] = sorted({m["signature"] for m in ms if m["analysis"] == "Pfam"})
        gene_tx = {}
        for tx in self.matches:
            gene_tx.setdefault(universe.transcript_to_gene.get(tx, tx), []).append(tx)
        for gene, txs in gene_tx.items():
            # positions from the longest protein with matches
            longest = max(txs, key=lambda t: seq_len.get(t, 0))
            first = {acc: s for s, acc in self.tx_entries.get(longest, [])}
            for t in txs:
                for s, acc in self.tx_entries.get(t, []):
                    first.setdefault(acc, s)
            self.gene_entries[gene] = sorted((s, acc) for acc, s in first.items())
            self.gene_pfam[gene] = sorted({p for t in txs for p in self.tx_pfam.get(t, [])})

        # GO terms of the InterPro entries
        n_genes_go = n_terms_added = n_genes_only_ipr = 0
        if self.go_mode != "off" and ctx.go_hierarchy is not None:
            gh = ctx.go_hierarchy
            from goannotate.go import GO_ROOT_TERMS
            for gene, entries in self.gene_entries.items():
                annot = ctx.annotations.get(gene)
                if annot is None:
                    continue
                if self.go_mode == "fill" and annot.go_terms:
                    continue
                ipr_terms = {}
                for _, acc in entries:
                    for t in self.ipr_go.get(acc, []):
                        norm = gh.normalize_id(t)
                        if not norm or norm in GO_ROOT_TERMS or norm in ctx.forbidden_go:
                            continue
                        ipr_terms.setdefault(norm, set()).add(acc)
                if ctx.namespaces and ipr_terms:
                    keep = gh.filter_by_namespace(set(ipr_terms), ctx.namespaces)
                    ipr_terms = {t: v for t, v in ipr_terms.items() if t in keep}
                if not ipr_terms:
                    continue
                had_go = bool(annot.go_terms)
                before = set(annot.go_terms)
                annot.go_terms = set(annot.go_terms) | set(ipr_terms)
                annot.specific_go_terms = gh.filter_to_specific(annot.go_terms)
                for t, accs in ipr_terms.items():
                    annot.go_interpro[t] = sorted(accs)
                    if t not in before:
                        annot.go_evidence[t] = "IEA"
                # Terms in the specific set that came from InterPro alone; terms already
                # carried by the hits, or ancestors of the gene's terms, add nothing.
                added = sorted(t for t in annot.specific_go_terms if t in ipr_terms and t not in before)
                if added:
                    self.gene_go_added[gene] = added
                    n_genes_go += 1
                    n_terms_added += len(added)
                    if not had_go:
                        n_genes_only_ipr += 1

        # Domain-based product names for coding genes without a product
        n_products = 0
        if self.product_mode != "off":
            for gene, entries in self.gene_entries.items():
                annot = ctx.annotations.get(gene)
                if annot is None or annot.consensus_name:
                    continue
                if universe.coding and not universe.coding.get(gene, True):
                    continue
                named = [(s, self.entries[acc][0], self.entries[acc][1]) for s, acc in entries
                         if acc in self.entries and self.entries[acc][1] in DOMAIN_TYPES]
                product = domain_product(named)
                if product:
                    annot.consensus_name = product
                    annot.product_note = "interpro_domain"
                    n_products += 1

        genes_with_entries = sum(1 for v in self.gene_entries.values() if v)
        self.counts = {"proteins_with_matches": len(self.matches), "genes_with_matches": len(gene_tx),
                       "genes_with_interpro": genes_with_entries,
                       "genes_with_pfam": sum(1 for v in self.gene_pfam.values() if v),
                       "genes_go_added": n_genes_go, "go_terms_added": n_terms_added,
                       "genes_go_only_interpro": n_genes_only_ipr, "products_filled": n_products}
        self.report = ctx.reports_dir / f"{ctx.prefix}_interproscan.tsv"
        self._write_report(ctx)
        logger.info(f"  InterPro entries on {genes_with_entries:,} genes; GO terms added to {n_genes_go:,} genes "
                    f"({n_terms_added:,} specific terms; {n_genes_only_ipr:,} genes had no GO from hits); "
                    f"domain-based products for {n_products:,} genes")
        return self.counts

    def _write_report(self, ctx: ModuleContext) -> None:
        with open(self.report, "w") as fh:
            fh.write("gene_id\ttranscript_id\tanalysis\tsignature\tsignature_description\tstart\tend\tscore\t"
                     "interpro\tinterpro_description\tinterpro_type\tgo_terms\n")
            for tx in sorted(self.matches):
                gene = ctx.universe.transcript_to_gene.get(tx, tx)
                for m in sorted(self.matches[tx], key=lambda m: (m["start"], m["analysis"])):
                    desc, etype = self.entries.get(m["interpro"], ("", ""))
                    go = self.ipr_go.get(m["interpro"], []) if m["interpro"] else []
                    fh.write("\t".join([gene, tx, m["analysis"], m["signature"], m["signature_description"],
                                        str(m["start"]), str(m["end"]), m["score"], m["interpro"],
                                        m["interpro_description"] or desc, etype, "|".join(go)]) + "\n")

    # ---- outputs ------------------------------------------------------------------

    def _entries_text(self, entries: list) -> str:
        return ";".join(f"{acc}:{self.entries.get(acc, ('', ''))[0]}" for _, acc in entries)

    def gene_columns(self) -> list:
        return [("interpro", lambda g: self._entries_text(self.gene_entries.get(g, []))),
                ("pfam", lambda g: ",".join(self.gene_pfam.get(g, []))),
                ("go_from_interpro", lambda g: ";".join(self.gene_go_added.get(g, [])))]

    def transcript_columns(self) -> list:
        return [("interpro", lambda t: self._entries_text(self.tx_entries.get(t, []))),
                ("pfam", lambda t: ",".join(self.tx_pfam.get(t, [])))]

    def gene_attributes(self, gene_key: str) -> dict:
        xrefs = [f"InterPro:{acc}" for _, acc in self.gene_entries.get(gene_key, [])]
        xrefs += [f"Pfam:{p}" for p in self.gene_pfam.get(gene_key, [])]
        return {"Dbxref": ",".join(xrefs)} if xrefs else {}

    def summary_lines(self) -> list:
        c = self.counts
        lines = [f"  InterProScan ({self.status.version}): {c.get('proteins_with_matches', 0):,} proteins with matches"
                 + (f", {self.runtime / 60:.1f} min" if self.runtime else " (results reused)")
                 + ("; all analyses run locally" if self.local else "; precalculated matches looked up at EBI")]
        lines.append(f"    Genes with InterPro entries: {c.get('genes_with_interpro', 0):,}; with Pfam families: "
                     f"{c.get('genes_with_pfam', 0):,}")
        lines.append(f"    InterPro GO terms ({self.go_mode}): {c.get('go_terms_added', 0):,} specific terms added to "
                     f"{c.get('genes_go_added', 0):,} genes ({c.get('genes_go_only_interpro', 0):,} without GO from hits)")
        lines.append(f"    Domain-based products ({self.product_mode}): {c.get('products_filled', 0):,} genes")
        if self.report:
            lines.append(f"    Matches:                     reports/{self.report.name}")
        return lines

    def report_files(self) -> list:
        return [self.report] if self.report else []


def _workflow_version(path: Path) -> str:
    """Release tag of a local InterProScan 6 clone, or its manifest version."""
    try:
        out = subprocess.run(["git", "-C", str(path), "describe", "--tags", "--exact-match"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=30)
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        m = re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]", (path / "nextflow.config").read_text())
        if m:
            return m.group(1)
    except OSError:
        pass
    return ""


def _nextflow_progress(log: Path) -> tuple:
    """(submitted, completed) task counts from .nextflow.log in the run directory."""
    submitted = done = 0
    try:
        with open(log, errors="replace") as fh:
            for line in fh:
                if "Submitted process >" in line:
                    submitted += 1
                elif "Task completed >" in line:
                    done += 1
    except OSError:
        pass
    return submitted, done
