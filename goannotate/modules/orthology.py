"""Orthology module: symbols from close relatives through OrthoFinder orthologs.

Close relatives with NCBI RefSeq annotations are chosen automatically from
the query organism's lineage (--taxid), or named with --relative. Their
proteomes and annotations are downloaded with the NCBI datasets
command-line tool; one protein per gene (the longest) of the query and of
each relative is given to OrthoFinder (Emms & Kelly 2019, Genome Biology).

For every query gene the module records its orthogroup and its orthologs in
each relative with their symbols. With --ortholog-symbol fill (default), a
gene without a symbol, or with a lineage symbol kept only because no
reference ortholog was found (foreign_kept), takes the symbol of its
one-to-one orthologs when the relatives that name them agree; when it has
no product either, it takes their product. Genes that already carry a
supported symbol are compared with their one-to-one orthologs' symbols, as
a check of both.

Relatives' symbols are put in the reference nomenclature first, as the core
does for hit symbols: a relative's gene (its NCBI GeneID) is mapped through
NCBI gene_orthologs to its ortholog in the reference species, whose symbol
is used (mouse Zfp558 -> ZNF558 with a human reference). A gene without a
reference ortholog keeps its own symbol.

With --same-species, the RefSeq annotation of the query species itself (a
previous annotation of the same genome or of another assembly) is added to
the relatives, beyond --relatives, so that its names carry over: its
one-to-one ortholog's symbol is used when the other relatives agree with it
or name no one-to-one ortholog; when they agree on a different symbol,
nothing is filled and the gene is counted as a disagreement. Supported core
symbols are also compared with it, as a continuity check.

Choosing relatives automatically: the lineage is walked outward from the
genus; at each level the NCBI reference genomes with a RefSeq annotation are
listed (newest annotation first), excluding the query species, and a
candidate is accepted when at least --relative-min-named of its
protein-coding genes have a symbol rather than an identifier (LOC numbers;
the symbol patterns), measured from its GFF3. The first --relatives
accepted are used. Every candidate examined is listed in the report.
"""

import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import time
import zipfile
from collections import Counter
from pathlib import Path

from goannotate.modules.base import Module, ModuleContext, ModuleError, ToolStatus, clean_sequence
from goannotate.sequences import read_fasta

logger = logging.getLogger(__name__)

PROGRESS_INTERVAL = 60
MAX_CANDIDATES_PER_LEVEL = 10
FILL_STATUSES = {"", "none", "foreign_kept"}          # symbol statuses the relatives may fill
SUPPORTED_STATUSES = {"rbh", "rbh_override", "rbh_added", "vote_supported", "reference_best_hit",
                      "foreign_ortholog"}              # symbols checked against the relatives
_GENOME_SUFFIX = re.compile(r"(?:\.\d+|L?[AB]|L\d*)$", re.I)


def species_label(name: str) -> str:
    """'Microtus oregoni' -> 'Microtus_oregoni' (file and column names)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def symbol_base(symbol: str) -> str:
    """Upper case without a genome suffix (.2, a, b, l, la, lb, l2): foxo3b -> FOXO3."""
    s = symbol.upper()
    stripped = _GENOME_SUFFIX.sub("", s)
    return stripped or s


def agreed_symbol(symbols: list) -> tuple:
    """
    (symbol, how) from the informative symbols of a gene's one-to-one
    orthologs in the relatives that name them: 'agree' when they are
    identical (case-insensitive), 'agree_base' when they differ only by a
    genome suffix (the most frequent form is used, the suffix-free form on a
    tie), ('', 'disagree') otherwise, ('', 'none') without symbols.
    """
    if not symbols:
        return "", "none"
    upper = Counter(s.upper() for s in symbols)
    if len(upper) == 1:
        return symbols[0], "agree"
    if len({symbol_base(s) for s in symbols}) == 1:
        (top, n), *rest = upper.most_common()
        if rest and rest[0][1] == n:
            return symbol_base(symbols[0]), "agree_base"
        return next(s for s in symbols if s.upper() == top), "agree_base"
    return "", "disagree"


def parse_refseq_gff(path: Path) -> dict:
    """
    RefSeq GFF3 -> {protein_id: (gene_id, symbol, product, biotype, NCBI GeneID)}.
    The CDS's protein_id, gene (symbol), product and Dbxref GeneID are used;
    the gene is found through the CDS parent (an mRNA whose parent is the
    gene, or the gene itself), and its biotype from the gene feature.
    """
    biotype = {}
    parent_of = {}
    cds = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            ftype, attrs = parts[2], parts[8]
            fid = _attr(attrs, "ID")
            if ftype in ("gene", "pseudogene"):
                biotype[fid] = _attr(attrs, "gene_biotype") or ftype
            elif ftype == "CDS":
                pid = _attr(attrs, "protein_id")
                if pid and pid not in cds:
                    m = re.search(r"GeneID:(\d+)", _attr(attrs, "Dbxref"))
                    cds[pid] = (_attr(attrs, "Parent").split(",")[0], _attr(attrs, "gene"),
                                _unquote(_attr(attrs, "product")), m.group(1) if m else "")
            elif fid:
                parent = _attr(attrs, "Parent")
                if parent:
                    parent_of[fid] = parent.split(",")[0]
    proteins = {}
    for pid, (parent, symbol, product, geneid) in cds.items():
        gene = parent if parent in biotype else parent_of.get(parent, parent)
        proteins[pid] = (gene, symbol, product, biotype.get(gene, ""), geneid)
    return proteins


def _attr(attrs: str, key: str) -> str:
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]*)", attrs)
    return m.group(1) if m else ""


def _unquote(value: str) -> str:
    from urllib.parse import unquote
    return unquote(value)


def named_fraction(gff: Path, is_bad_symbol) -> tuple:
    """(protein-coding genes, fraction with an informative symbol) of a RefSeq GFF3."""
    n = named = 0
    with open(gff) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "gene" or "gene_biotype=protein_coding" not in parts[8]:
                continue
            n += 1
            sym = _attr(parts[8].rstrip("\n"), "gene")
            if sym and not is_bad_symbol(sym):
                named += 1
    return n, (named / n if n else 0.0)


def parse_orthologues(path: Path) -> list:
    """OrthoFinder query__v__species.tsv -> [(orthogroup, [query genes], [relative genes])]."""
    rows = []
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q = [g.strip() for g in parts[1].split(",") if g.strip()]
            r = [g.strip() for g in parts[2].split(",") if g.strip()]
            rows.append((parts[0], q, r))
    return rows


def parse_orthologs_by_species(path: Path) -> dict:
    """OrthoFinder 3 Orthologs/query.tsv(.gz) -> {species: [(orthogroup, [query genes], [relative genes])]}
    with one row per query gene, that gene first in the query list.

    OrthoFinder 3 writes one long table per species (columns Orthogroup,
    Species, query, Orthologs) in which the orthologs of one query gene may
    be split over several rows. Rows are therefore read as gene pairs and
    regrouped per query gene: the relative genes are all of its orthologs,
    and the query genes are every query gene that shares one of them, so
    that a 1:1 ortholog has one gene on each side, as in OrthoFinder 2."""
    opener = gzip.open if str(path).endswith(".gz") else open
    pairs = {}                                     # species -> (query -> {rel}, rel -> {query}, query -> og)
    with opener(path, "rt") as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            og, species = parts[0], parts[1]
            q_to_r, r_to_q, og_of = pairs.setdefault(species, ({}, {}, {}))
            qs = [g.strip() for g in parts[2].split(",") if g.strip()]
            rs = [g.strip() for g in parts[3].split(",") if g.strip()]
            for q in qs:
                og_of[q] = og
                q_to_r.setdefault(q, set()).update(rs)
            for r in rs:
                r_to_q.setdefault((og, r), set()).update(qs)
    out = {}
    for species, (q_to_r, r_to_q, og_of) in pairs.items():
        rows = []
        for q in sorted(q_to_r):
            rel = sorted(q_to_r[q])
            co = sorted(set().union(*(r_to_q[(og_of[q], r)] for r in rel)) - {q})
            rows.append((og_of[q], [q] + co, rel))
        out[species] = rows
    return out


def parse_orthogroups(path: Path, species: str) -> dict:
    """Orthogroups.tsv -> {gene of the species: orthogroup}."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if species not in header:
            return out
        col = header.index(species)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if col < len(parts):
                for g in parts[col].split(","):
                    g = g.strip()
                    if g:
                        out[g] = parts[0]
    return out


def cpu_temperature() -> str:
    """Highest CPU package temperature reported by `sensors` (Tctl, Tdie or Package), or ''."""
    if not shutil.which("sensors"):
        return ""
    try:
        out = subprocess.run(["sensors"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                             timeout=10).stdout
    except (OSError, subprocess.TimeoutExpired):
        return ""
    temps = [float(m) for m in re.findall(r"(?:Tctl|Tdie|Package id \d+):\s*\+?([0-9.]+)", out)]
    return f"{max(temps):.0f} C" if temps else ""


class RelGene:
    """A relative's representative protein: its gene, symbol, product, NCBI GeneID and the symbol in
    the reference nomenclature (the reference-species ortholog's symbol, else its own)."""
    __slots__ = ("gene", "symbol", "product", "geneid", "ref_symbol")

    def __init__(self, gene="", symbol="", product="", geneid="", ref_symbol=None):
        self.gene, self.symbol, self.product, self.geneid = gene, symbol, product, geneid
        self.ref_symbol = symbol if ref_symbol is None else ref_symbol


_NO_GENE = RelGene()


class Relative:
    def __init__(self, accession, organism, taxid, level="", annotation="", release="", protein_coding=0):
        self.accession = accession
        self.organism = organism
        self.taxid = taxid
        self.level = level                  # lineage level at which it was found
        self.annotation = annotation
        self.release = release
        self.protein_coding = protein_coding
        self.named = None                   # fraction of protein-coding genes with a symbol
        self.status = ""                    # chosen, rejected_named, not_needed, failed
        self.label = species_label(organism)
        self.proteins = {}                  # representative protein -> RelGene
        self.same_species = False           # the query species' own annotation (--same-species)
        self.n_genes = 0


class OrthologyModule(Module):
    name = "orthology"
    description = "OrthoFinder orthologs in close relatives (NCBI RefSeq); symbols for unnamed genes"

    def __init__(self, args):
        super().__init__(args)
        self.orthofinder = getattr(args, "orthofinder_bin", None)
        self.datasets = getattr(args, "datasets_bin", None)
        self.n_relatives = getattr(args, "relatives", 3) or 3
        self.manual = getattr(args, "relative", None) or []
        self.min_named = getattr(args, "relative_min_named", 0.7)
        self.cache = getattr(args, "relatives_dir", None)
        self.of_threads = getattr(args, "orthofinder_threads", 0) or 0
        self.mode = getattr(args, "ortholog_symbol", "fill")
        self.taxid = getattr(args, "taxid", None)
        self.db = getattr(args, "db", None)
        self.with_same_species = getattr(args, "same_species", False)
        self.same = None                    # Relative: the query species' own annotation
        self.candidates = []                # every Relative examined
        self.relatives = []                 # chosen Relatives
        self.orthogroup = {}                # query gene -> orthogroup
        self.orthologs = {}                 # query gene -> {relative label: (type, [(gene, symbol)])}
        self.gene_symbol = {}               # query gene -> (agreed 1:1 symbol, how)
        self.filled = {}                    # query gene -> (old symbol, old status)
        self.check = Counter()
        self.same_check = Counter()         # supported symbols against the same-species annotation
        self.normalized = Counter()         # relatives' symbols put in the reference nomenclature
        self.report = None

    # ---- detection ------------------------------------------------------------

    def _env(self, of_path: str) -> dict:
        env = dict(os.environ)
        # OrthoFinder's script runs with `env python3` and needs its environment's tools (diamond, mcl)
        env["PATH"] = str(Path(of_path).parent) + os.pathsep + env.get("PATH", "")
        return env

    def detect(self) -> ToolStatus:
        name = "OrthoFinder + NCBI datasets"
        of = self.find_executable(self.orthofinder, "orthofinder")
        if not of:
            self.status = ToolStatus(name, False, reason="orthofinder not found (give --orthofinder-bin)")
            return self.status
        ds = self.find_executable(self.datasets, "datasets")
        if not ds:
            self.status = ToolStatus(name, False, path=of, reason="NCBI datasets not found (give --datasets-bin)")
            return self.status
        if not self.manual and not self.taxid:
            self.status = ToolStatus(name, False, path=of, reason="give --taxid (relatives are chosen from its "
                                                                  "lineage) or --relative")
            return self.status
        if self.with_same_species and not self.taxid:
            self.status = ToolStatus(name, False, path=of, reason="--same-species needs --taxid")
            return self.status
        try:
            out = subprocess.run([of, "-h"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                 timeout=120, env=self._env(of)).stdout
        except (OSError, subprocess.TimeoutExpired) as e:
            out = str(e)
        if "orthofinder" not in out.lower() or "Traceback" in out:
            self.status = ToolStatus(name, False, path=of, reason=f"{of} does not run: {out.strip()[-200:]}")
            return self.status
        try:
            vout = subprocess.run([of, "-v"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                  timeout=120, env=self._env(of)).stdout
        except (OSError, subprocess.TimeoutExpired):
            vout = ""
        vout = re.sub(r"\x1b\[[0-9;]*m", "", vout)          # "OrthoFinder:v3.2" in colour
        m = re.search(r"OrthoFinder:?\s*v?(\d+\.\d+(?:\.\d+)?)", vout)
        rc, dsv = self.run_command([ds, "--version"], timeout=60)
        dsv = dsv.strip().split()[-1] if rc == 0 and dsv.strip() else ""
        self.datasets = ds
        self.status = ToolStatus(name, True, path=of, version=f"OrthoFinder {m.group(1) if m else '?'}, "
                                                              f"datasets {dsv or '?'}")
        return self.status

    # ---- relatives --------------------------------------------------------------

    def _summary(self, taxon: str) -> list:
        """Annotated RefSeq reference genomes of a taxon (datasets summary), newest annotation first."""
        rc, out = self.run_command([self.datasets, "summary", "genome", "taxon", str(taxon), "--reference",
                                    "--annotated", "--assembly-source", "RefSeq", "--as-json-lines"], timeout=600)
        if rc != 0:
            logger.warning(f"  datasets summary for taxon {taxon} failed: {out.strip()[:200]}")
            return []
        rows = []
        for line in out.splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            d = json.loads(line)
            ann = d.get("annotation_info", {})
            rows.append(Relative(d.get("accession", ""), d.get("organism", {}).get("organism_name", ""),
                                 int(d.get("organism", {}).get("tax_id", 0) or 0), annotation=ann.get("name", ""),
                                 release=ann.get("release_date", ""),
                                 protein_coding=ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", 0)))
        rows.sort(key=lambda r: (r.release, r.protein_coding), reverse=True)
        return rows

    def _download(self, rel: Relative, cache: Path, parts: str) -> Path:
        """Download (or reuse) a RefSeq annotation into cache/<accession>/; returns that directory."""
        d = cache / rel.accession
        need = {"gff3": d / "genomic.gff", "protein": d / "protein.faa"}
        wanted = [p for p in parts.split(",")]
        if all(need[p].exists() and need[p].stat().st_size > 0 for p in wanted):
            return d
        d.mkdir(parents=True, exist_ok=True)
        zpath = d / f"{rel.accession}_{'_'.join(wanted)}.zip"
        rc, out = self.run_command([self.datasets, "download", "genome", "accession", rel.accession,
                                    "--include", ",".join(wanted), "--filename", str(zpath), "--no-progressbar"],
                                   timeout=3600)
        if rc != 0 or not zpath.exists():
            raise ModuleError(f"datasets download of {rel.accession} failed: {out.strip()[:200]}")
        names = {"gff3": "genomic.gff", "protein": "protein.faa"}
        with zipfile.ZipFile(zpath) as z:
            for p in wanted:
                member = f"ncbi_dataset/data/{rel.accession}/{names[p]}"
                with z.open(member) as src, open(need[p], "wb") as dst:
                    shutil.copyfileobj(src, dst)
        zpath.unlink()
        return d

    def _choose_relatives(self, ctx: ModuleContext, cache: Path, is_bad_symbol) -> None:
        self._choose_other_relatives(ctx, cache, is_bad_symbol)
        if self.with_same_species:
            self._choose_same_species(cache, is_bad_symbol)

    def _choose_same_species(self, cache: Path, is_bad_symbol) -> None:
        """The newest RefSeq annotation of the query species itself, added beyond --relatives."""
        from goannotate import taxon
        query_name = taxon.taxon_name_sqlite(self.db, int(self.taxid)) if self.db else ""
        rows = [r for r in self._summary(self.taxid)
                if r.taxid == int(self.taxid) or (query_name and r.organism.startswith(query_name))]
        rows = [r for r in rows if r.accession not in {c.accession for c in self.relatives}]
        if not rows:
            logger.warning(f"  --same-species: no RefSeq annotation of taxon {self.taxid} was found")
            return
        rel = rows[0]
        rel.level = "same species"
        rel.same_species = True
        self.candidates.append(rel)
        try:
            self._evaluate(rel, cache, is_bad_symbol)
        except ModuleError as e:
            rel.status = "failed"
            logger.warning(f"  {e}")
            return
        rel.status = "chosen"
        if any(r.label == rel.label for r in self.relatives):      # a subspecies named like a relative
            rel.label += "_same_species"
        self.same = rel
        self.relatives.append(rel)

    def _choose_other_relatives(self, ctx: ModuleContext, cache: Path, is_bad_symbol) -> None:
        from goannotate import taxon
        if self.manual:
            for value in self.manual:
                value = value.strip()
                if re.match(r"^GC[AF]_\S+$", value):
                    rows = [r for r in self._summary_accession(value)]
                else:
                    rows = self._summary(value)[:1]
                if not rows:
                    raise ModuleError(f"--relative {value}: no annotated RefSeq genome found")
                rel = rows[0]
                rel.level = "--relative"
                self._evaluate(rel, cache, is_bad_symbol)
                rel.status = "chosen"
                self.candidates.append(rel)
                self.relatives.append(rel)
            return
        lineage = taxon.lineage_sqlite(self.db, int(self.taxid)) if self.db else []
        if not lineage:
            raise ModuleError(f"no lineage for taxon {self.taxid} (needs the taxonomy tables in --db)")
        query_name = taxon.taxon_name_sqlite(self.db, lineage[0]) or ""
        seen = set()
        for level in lineage[1:]:
            level_name = taxon.taxon_name_sqlite(self.db, level) or str(level)
            rows = [r for r in self._summary(level) if r.accession not in seen]
            examined = 0
            for rel in rows:
                seen.add(rel.accession)
                if rel.taxid == lineage[0] or (query_name and rel.organism.startswith(query_name)):
                    continue
                if any(c.taxid == rel.taxid for c in self.candidates):
                    continue
                if examined >= MAX_CANDIDATES_PER_LEVEL:
                    break
                examined += 1
                rel.level = level_name
                self.candidates.append(rel)
                try:
                    self._evaluate(rel, cache, is_bad_symbol)
                except ModuleError as e:
                    rel.status = "failed"
                    logger.warning(f"  {e}")
                    continue
                if rel.named >= self.min_named:
                    rel.status = "chosen"
                    self.relatives.append(rel)
                    if len(self.relatives) >= self.n_relatives:
                        return
                else:
                    rel.status = "rejected_named"
            if len(self.relatives) >= self.n_relatives:
                return

    def _summary_accession(self, accession: str) -> list:
        rc, out = self.run_command([self.datasets, "summary", "genome", "accession", accession, "--as-json-lines"],
                                   timeout=300)
        rows = []
        if rc == 0:
            for line in out.splitlines():
                if line.strip().startswith("{"):
                    d = json.loads(line)
                    ann = d.get("annotation_info", {})
                    rows.append(Relative(d.get("accession", accession), d.get("organism", {}).get("organism_name", ""),
                                         int(d.get("organism", {}).get("tax_id", 0) or 0),
                                         annotation=ann.get("name", ""), release=ann.get("release_date", ""),
                                         protein_coding=ann.get("stats", {}).get("gene_counts", {})
                                         .get("protein_coding", 0)))
        return rows

    def _evaluate(self, rel: Relative, cache: Path, is_bad_symbol) -> None:
        d = self._download(rel, cache, "gff3")
        n, frac = named_fraction(d / "genomic.gff", is_bad_symbol)
        rel.named = frac
        if not rel.protein_coding:
            rel.protein_coding = n
        logger.info(f"  Candidate relative {rel.organism} ({rel.accession}, {rel.level}): {n:,} protein-coding "
                    f"genes, {100 * frac:.1f}% with a symbol")

    def _load_relative(self, rel: Relative, cache: Path) -> Path:
        """Representative proteins (longest per gene) of a chosen relative; returns its FASTA for OrthoFinder."""
        d = self._download(rel, cache, "gff3,protein")
        info = parse_refseq_gff(d / "genomic.gff")
        best = {}
        seqs = {}
        for header, seq in read_fasta(str(d / "protein.faa")):
            pid = header.split()[0]
            rec = info.get(pid)
            if rec is None or (rec[3] and rec[3] != "protein_coding"):
                continue
            gene = rec[0]
            s = clean_sequence(seq)
            if gene not in best or len(s) > len(seqs[best[gene]]):
                best[gene] = pid
                seqs[pid] = s
        rel.proteins = {pid: RelGene(info[pid][0], info[pid][1], info[pid][2], info[pid][4])
                        for pid in best.values()}
        rel.n_genes = len(best)
        # only the longest protein of each gene (seqs also holds proteins that were longest earlier)
        return {pid: seqs[pid] for pid in best.values()}

    def _normalize_symbols(self, is_bad=None) -> None:
        """
        Put the relatives' symbols in the reference nomenclature: each gene's
        NCBI GeneID is mapped through gene_orthologs (table ncbi_orthologs) to
        its ortholog in the reference species, whose symbol is used (mouse
        Zfp558 -> ZNF558). The relative's own symbol is kept when it has no
        reference ortholog, when that ortholog is a pseudogene (CDHR18P for a
        coding Cdhr18) or a readthrough locus (KLRC4-KLRK1, whose parts are
        both reference symbols). Former names are updated (vole Sept2 ->
        SEPTIN2); comparisons with the core's symbols treat a former name of
        the same reference gene as agreement (_equivalent), since the core's
        UniProt names and NCBI gene_info are not always renamed at the same
        time (TNFAIP8L3 / TIPE3). A LOC identifier is replaced by its reference ortholog's
        symbol under the same rules (counted apart); otherwise identifiers
        are not counted.
        """
        from goannotate import mappings
        args = self.args
        ref = str(getattr(args, "reference_taxid", "") or "")
        anchor = str(getattr(args, "ortholog_anchor_taxid", "") or "")
        if getattr(args, "no_resolve_orthologs", False):
            logger.info("  Relatives' symbols kept as given (--no-resolve-orthologs)")
            return
        if not (self.db and ref.isdigit() and anchor.isdigit()):
            logger.info("  Relatives' symbols kept as given (no reference taxon or --db)")
            return
        if not (mappings.sqlite_has_table(self.db, "ncbi_orthologs")
                and mappings.sqlite_has_table(self.db, "ncbi_geneinfo")):
            logger.warning("  The database has no ncbi_orthologs or ncbi_geneinfo table; relatives' symbols "
                           "are kept as given (run build_databases.py --add-orthologs)")
            return
        is_bad = is_bad or (lambda s: s.upper().startswith("LOC"))
        genes = [g for r in self.relatives if str(r.taxid) != ref for g in r.proteins.values() if g.geneid]
        geneids = {g.geneid for g in genes}
        to_ref = mappings.load_ncbi_orthologs_sqlite(self.db, geneids, anchor, ref) if geneids else {}
        ref_symbol = mappings.load_ncbi_geneinfo_sqlite(self.db, set(to_ref.values())) if to_ref else {}
        pseudo = _pseudogenes(self.db, set(ref_symbol))
        # official reference symbols, to recognise readthroughs (both parts are symbols)
        lookup = set()
        for sym in ref_symbol.values():
            if "-" in sym:
                lookup.update(p.lower() for p in sym.split("-"))
        names = _reference_names(self.db, ref, lookup)
        for g in genes:
            if not g.symbol:
                continue
            ref_gid = to_ref.get(g.geneid, "")
            new = ref_symbol.get(ref_gid, "")
            if not new:
                outcome = "no_reference_ortholog"
            elif ref_gid in pseudo:
                outcome = "reference_pseudogene"
            elif new.upper() == g.symbol.upper():
                outcome = "same"
            elif "-" in new and all(names.get(p.lower(), {}).get("symbol") for p in new.split("-")):
                outcome = "reference_readthrough"
            else:
                outcome = "changed"
                g.ref_symbol = new
            if is_bad(g.symbol):
                if outcome == "changed":
                    self.normalized["identifier_replaced"] += 1
                continue
            self.normalized[outcome] += 1
        n = self.normalized
        logger.info(f"  Relatives' symbols in the reference nomenclature (taxon {ref}): {n['same']:,} unchanged, "
                    f"{n['changed']:,} replaced by the reference ortholog's symbol "
                    f"(and {n['identifier_replaced']:,} LOC identifiers); kept: "
                    f"{n['no_reference_ortholog']:,} without a reference ortholog, "
                    f"{n['reference_pseudogene']:,} with a pseudogene and {n['reference_readthrough']:,} with a "
                    f"readthrough as reference ortholog")

    # ---- running -------------------------------------------------------------------

    def run(self, ctx: ModuleContext) -> None:
        if self.status is None:
            self.detect()
        if not self.status.available:
            raise ModuleError(self.status.reason)
        base = (ctx.work_dir / f"{ctx.prefix}_orthology").resolve()
        base.mkdir(parents=True, exist_ok=True)
        cache = Path(self.cache).expanduser().resolve() if self.cache else base / "relatives"
        cache.mkdir(parents=True, exist_ok=True)
        is_bad = ctx.bad_name_filter.is_bad_symbol if ctx.bad_name_filter else (lambda s: s.startswith("LOC"))

        self._choose_relatives(ctx, cache, is_bad)
        if not self.relatives:
            raise ModuleError("no relative with an annotation naming at least "
                              f"{100 * self.min_named:.0f}% of its protein-coding genes was found")
        logger.info("  Relatives: " + "; ".join(f"{r.organism} ({r.accession}, {r.level}, "
                                                 f"{100 * r.named:.1f}% named)" for r in self.relatives))

        # OrthoFinder input: one protein per gene
        indir = base / "input"
        if indir.exists():
            shutil.rmtree(indir)
        indir.mkdir()
        query = {}
        for header, seq in read_fasta(str(ctx.query_fasta)):
            tx = header.split()[0]
            gene = ctx.universe.transcript_to_gene.get(tx, tx)
            s = clean_sequence(seq)
            if s and (gene not in query or len(s) > len(query[gene])):
                query[gene] = s
        _write_fasta(indir / "query.faa", query)
        for rel in self.relatives:
            _write_fasta(indir / f"{rel.label}.faa", self._load_relative(rel, cache))
        logger.info(f"  OrthoFinder input: query {len(query):,} genes; "
                    + ", ".join(f"{r.label} {r.n_genes:,}" for r in self.relatives))
        self._normalize_symbols(is_bad)

        signature = json.dumps({"query": len(query), "relatives": [r.accession for r in self.relatives]})
        results = base / "orthofinder" / "Results_goannotate"
        marker = base / "orthofinder.done"
        if not (results.exists() and marker.exists() and marker.read_text() == signature and not ctx.force):
            if (base / "orthofinder").exists():
                shutil.rmtree(base / "orthofinder")
            threads = self.of_threads or max(1, ctx.threads or 1)
            cmd = [self.status.path, "-f", str(indir), "-o", str(base / "orthofinder"), "-n", "goannotate",
                   "-t", str(threads), "-a", str(max(1, threads // 4)), "-S", "diamond"]
            log = base / "orthofinder.log"
            logger.info(f"  Running OrthoFinder with {threads} threads: {' '.join(cmd)}")
            t0 = time.time()
            last = t0
            with open(log, "w") as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=self._env(self.status.path))
                while proc.poll() is None:
                    time.sleep(5)
                    if time.time() - last >= PROGRESS_INTERVAL:
                        last = time.time()
                        temp = cpu_temperature()
                        logger.info(f"  OrthoFinder running, {(last - t0) / 60:.1f} min"
                                    + (f"; CPU {temp}" if temp else "") + f"; {_last_line(log)}")
            self.runtime = time.time() - t0
            if proc.returncode != 0:
                raise ModuleError(f"OrthoFinder failed (exit {proc.returncode}); see {log}")
            found = list((base / "orthofinder").glob("**/Orthogroups/Orthogroups.tsv"))
            if not found:
                raise ModuleError(f"OrthoFinder wrote no Orthogroups.tsv; see {log}")
            results = found[0].parent.parent
            marker.write_text(signature)
            logger.info(f"  OrthoFinder finished in {self.runtime / 60:.1f} min")
        else:
            logger.info(f"  Reusing OrthoFinder results in {results}")
        self._read_results(results)

    def _read_results(self, results: Path) -> None:
        og = list(results.glob("**/Orthogroups/Orthogroups.tsv"))
        self.orthogroup = parse_orthogroups(og[0], "query") if og else {}
        # OrthoFinder 3: Orthologs/query.tsv.gz holds every relative; OrthoFinder 2:
        # Orthologues/Orthologues_query/query__v__<species>.tsv, one file per relative.
        v3 = (list(results.glob("**/Orthologs/query.tsv.gz")) + list(results.glob("**/Orthologs/query.tsv")))
        by_species = parse_orthologs_by_species(v3[0]) if v3 else {}
        for rel in self.relatives:
            if v3:
                rows = by_species.get(rel.label, [])
            else:
                files = list(results.glob(f"**/Orthologues_query/query__v__{rel.label}.tsv"))
                rows = parse_orthologues(files[0]) if files else []
            if not rows:
                logger.warning(f"  No OrthoFinder orthologs for {rel.label}")
                continue
            for _, qgenes, rgenes in rows:
                kind = f"{'1' if len(qgenes) == 1 else 'many'}:{'1' if len(rgenes) == 1 else 'many'}"
                named = [(pid, rel.proteins.get(pid, _NO_GENE).ref_symbol) for pid in rgenes]
                for q in (qgenes[:1] if v3 else qgenes):     # OrthoFinder 3 rows are per query gene
                    self.orthologs.setdefault(q, {})[rel.label] = (kind, named)

    # ---- integration ------------------------------------------------------------------

    def integrate(self, ctx: ModuleContext) -> dict:
        is_bad = ctx.bad_name_filter.is_bad_symbol if ctx.bad_name_filter else (lambda s: s.startswith("LOC"))
        fmt = ctx.symbol_case_fn or (lambda s: s)
        rel_by_label = {r.label: r for r in self.relatives}
        same_label = self.same.label if self.same else None
        symbols_seen = {a.consensus_symbol.lower() for a in ctx.annotations.values() if a.consensus_symbol}
        symbols_seen |= {sym.lower() for by in self.orthologs.values() for _, named in by.values()
                         for _, sym in named if sym}
        ref = str(getattr(self.args, "reference_taxid", "") or "")
        self._names = _reference_names(self.db, ref, symbols_seen) if (self.db and ref.isdigit()) else {}
        n_filled = n_products = 0
        for gene, by_rel in self.orthologs.items():
            symbols = []
            products = []
            same_sym, same_prod = "", ""
            for label, (kind, named) in by_rel.items():
                if kind != "1:1":
                    continue
                pid, sym = named[0]
                if sym and not is_bad(sym):
                    prod = rel_by_label[label].proteins.get(pid, _NO_GENE).product
                    if label == same_label:
                        same_sym, same_prod = sym, prod
                        continue
                    symbols.append(sym)
                    if prod:
                        products.append(prod)
            sym, how = agreed_symbol(symbols)
            other_sym = sym                 # the check of supported symbols uses the other relatives only
            if same_sym:
                sym, how, products = self._with_same_species(same_sym, same_prod, sym, how, products)
            self.gene_symbol[gene] = (sym, how)
            annot = ctx.annotations.get(gene)
            if annot is None:
                continue
            if not sym and how in ("disagree", "same_species_disagree"):
                self.check["relatives_disagree"] += 1
            status = annot.symbol_status or ""
            if same_sym and annot.consensus_symbol and status in SUPPORTED_STATUSES:
                self.same_check["compared"] += 1
                self.same_check[self._equivalent(annot.consensus_symbol, same_sym) or "different"] += 1
            if annot.consensus_symbol and status in SUPPORTED_STATUSES:
                if other_sym:
                    self.check["compared"] += 1
                    self.check[self._equivalent(annot.consensus_symbol, other_sym) or "different"] += 1
                continue
            if not sym:
                continue
            if self.mode != "fill":
                continue
            if annot.consensus_symbol and status not in FILL_STATUSES:
                continue
            if annot.consensus_symbol and self._equivalent(annot.consensus_symbol, sym) in ("same", "synonym"):
                self.check["kept_same"] += 1
                continue
            old = annot.consensus_symbol
            self.filled[gene] = (old, status)
            if old and old not in annot.symbol_aliases:
                annot.symbol_aliases = list(annot.symbol_aliases) + [old]
            annot.consensus_symbol = fmt(sym)
            annot.symbol_status = "ortholog_fill"
            n_filled += 1
            if not annot.consensus_name and products:
                cleaned = [p for p in (ctx.bad_name_filter.process_name(x) if ctx.bad_name_filter else x
                                       for x in products) if p]
                if cleaned:
                    annot.consensus_name = Counter(cleaned).most_common(1)[0][0]
                    annot.product_note = "ortholog"
                    n_products += 1

        coding = sum(1 for v in ctx.universe.coding.values() if v) if ctx.universe.coding else len(ctx.annotations)
        with_og = sum(1 for g in self.orthogroup)
        per_rel = {}
        for rel in self.relatives:
            one = sum(1 for by in self.orthologs.values() if by.get(rel.label, ("",))[0] == "1:1")
            anyo = sum(1 for by in self.orthologs.values() if rel.label in by)
            per_rel[rel.label] = (one, anyo)
        self.per_relative = per_rel
        self.counts = {"relatives": len(self.relatives), "candidates_examined": len(self.candidates),
                       "query_genes": coding, "genes_in_orthogroups": with_og,
                       "genes_with_orthologs": len(self.orthologs),
                       "genes_with_agreed_symbol": sum(1 for s, _ in self.gene_symbol.values() if s),
                       "symbols_filled": n_filled,
                       "filled_unnamed": sum(1 for o, _ in self.filled.values() if not o),
                       "filled_foreign_kept": sum(1 for o, _ in self.filled.values() if o),
                       "products_filled": n_products, **{f"check_{k}": v for k, v in self.check.items()},
                       **{f"same_species_{k}": v for k, v in self.same_check.items()},
                       **{f"normalized_{k}": v for k, v in self.normalized.items()}}
        self.report = ctx.reports_dir / f"{ctx.prefix}_orthology.tsv"
        self._write_report(ctx)
        self._write_relatives(ctx)
        c = self.check
        agree = (c["same"] + c["same_base"] + c["synonym"]) / c["compared"] if c["compared"] else 0
        logger.info(f"  Orthologs: {len(self.orthologs):,} genes with orthologs in the relatives; symbols filled "
                    f"for {n_filled:,} genes ({self.counts['filled_unnamed']:,} unnamed, "
                    f"{self.counts['filled_foreign_kept']:,} lineage symbols), products for {n_products:,}; "
                    f"check against supported symbols: {100 * agree:.1f}% of {c['compared']:,} agree")
        if self.same:
            sc = self.same_check
            sagree = (sc["same"] + sc["same_base"] + sc["synonym"]) / sc["compared"] if sc["compared"] else 0
            logger.info(f"  Same species ({self.same.organism}, {self.same.accession}): supported symbols agree "
                        f"with it for {100 * sagree:.1f}% of {sc['compared']:,}; "
                        f"{sum(1 for _, h in self.gene_symbol.values() if h == 'same_species'):,} symbols from it "
                        f"alone, {sum(1 for _, h in self.gene_symbol.values() if h == 'same_species_disagree'):,} "
                        f"where the other relatives agree on a different symbol")
        return self.counts

    def _equivalent(self, a: str, b: str) -> str:
        """
        How two symbols name the same gene: 'same' (case-insensitive), 'same_base' (they differ by a
        genome suffix), 'synonym' (one is the reference species' official symbol of a gene that the
        other names as symbol or synonym: TIPE3 and its former name TNFAIP8L3), or '' when they differ.
        """
        if a.upper() == b.upper():
            return "same"
        if symbol_base(a) == symbol_base(b):
            return "same_base"
        na = getattr(self, "_names", {}).get(a.lower(), {})
        nb = getattr(self, "_names", {}).get(b.lower(), {})
        if (na.get("symbol", set()) & nb.get("any", set())) or (nb.get("symbol", set()) & na.get("any", set())):
            return "synonym"
        return ""

    def _with_same_species(self, same_sym: str, same_prod: str, sym: str, how: str, products: list) -> tuple:
        """
        Combine the same-species ortholog's symbol with the other relatives'
        result: it is used when they agree with it (same suffix-free form), name
        no one-to-one ortholog, or disagree among themselves; when they agree on
        a different symbol, nothing is used ('same_species_disagree').
        """
        prods = [same_prod] if same_prod else []
        if not sym and how in ("none", "disagree"):
            return same_sym, "same_species", prods
        if self._equivalent(sym, same_sym):
            return sym, "agree_same_species", prods + products
        return "", "same_species_disagree", []

    def _write_report(self, ctx: ModuleContext) -> None:
        labels = [r.label for r in self.relatives]
        with open(self.report, "w") as fh:
            fh.write("gene_id\tsymbol\tsymbol_status\torthogroup\tagreed_symbol\tagreement\tfilled_from\t"
                     + "\t".join(f"{l}_type\t{l}_orthologs" for l in labels) + "\n")
            for gene in sorted(set(self.orthogroup) | set(self.orthologs)):
                annot = ctx.annotations.get(gene)
                sym, how = self.gene_symbol.get(gene, ("", "none"))
                old = self.filled.get(gene)
                row = [gene, annot.consensus_symbol if annot else "", annot.symbol_status if annot else "",
                       self.orthogroup.get(gene, ""), sym, how,
                       (old[0] or "(none)") if old else ""]
                for l in labels:
                    kind, named = self.orthologs.get(gene, {}).get(l, ("", []))
                    rel = next(r for r in self.relatives if r.label == l)
                    row += [kind, ",".join(_ortholog_text(pid, rel.proteins.get(pid, _NO_GENE)) for pid, _ in named)]
                fh.write("\t".join(row) + "\n")

    def _write_relatives(self, ctx: ModuleContext) -> None:
        path = ctx.reports_dir / f"{ctx.prefix}_orthology_relatives.tsv"
        with open(path, "w") as fh:
            fh.write("accession\torganism\ttaxid\tlevel\tannotation\trelease\tprotein_coding_genes\tnamed_fraction\t"
                     "status\tgenes_used\tquery_genes_1to1\tquery_genes_with_orthologs\n")
            for r in self.candidates:
                one, anyo = self.per_relative.get(r.label, ("", "")) if r.status == "chosen" else ("", "")
                fh.write("\t".join(str(x) for x in [r.accession, r.organism, r.taxid, r.level, r.annotation, r.release,
                                                    r.protein_coding, f"{r.named:.3f}" if r.named is not None else "",
                                                    r.status, r.n_genes or "", one, anyo]) + "\n")

    # ---- outputs ----------------------------------------------------------------------

    def _orthologs_text(self, gene: str) -> str:
        parts = []
        for label, (kind, named) in self.orthologs.get(gene, {}).items():
            syms = ",".join(s for _, s in named if s) or ",".join(p for p, _ in named)
            parts.append(f"{label}:{syms}({kind})")
        return ";".join(parts)

    def gene_columns(self) -> list:
        return [("orthogroup", lambda g: self.orthogroup.get(g, "")),
                ("relative_orthologs", self._orthologs_text),
                ("ortholog_symbol", lambda g: self.gene_symbol.get(g, ("", ""))[0])]

    def gene_attributes(self, gene_key: str) -> dict:
        attrs = {}
        if gene_key in self.orthogroup:
            attrs["goannotate_orthogroup"] = self.orthogroup[gene_key]
        sym = self.gene_symbol.get(gene_key, ("", ""))[0]
        if sym:
            attrs["goannotate_ortholog_symbol"] = sym
        return attrs

    def summary_lines(self) -> list:
        c = self.counts
        lines = [f"  Orthology ({self.status.version}): {len(self.relatives)} relatives"
                 + (f", OrthoFinder {self.runtime / 60:.1f} min" if self.runtime else " (OrthoFinder results reused)")]
        for r in self.relatives:
            one, anyo = self.per_relative.get(r.label, (0, 0))
            lines.append(f"    {r.organism} ({r.accession}, {r.level}): {100 * r.named:.1f}% of {r.protein_coding:,} "
                         f"genes named; orthologs for {anyo:,} query genes, 1:1 for {one:,}"
                         + (" (same species, beyond --relatives)" if r.same_species else ""))
        lines.append(f"    Candidates examined:         {c.get('candidates_examined', 0)} "
                     f"(reports/{Path(self.report).name.replace('.tsv', '_relatives.tsv')})")
        lines.append(f"    Genes in orthogroups:        {c.get('genes_in_orthogroups', 0):,}")
        lines.append(f"    Symbols filled ({self.mode}):        {c.get('symbols_filled', 0):,} "
                     f"({c.get('filled_unnamed', 0):,} unnamed genes, {c.get('filled_foreign_kept', 0):,} lineage symbols); "
                     f"products {c.get('products_filled', 0):,}")
        compared = c.get("check_compared", 0)
        if compared:
            agree = c.get("check_same", 0) + c.get("check_same_base", 0) + c.get("check_synonym", 0)
            lines.append(f"    Check, supported symbols:    {agree:,} of {compared:,} agree with the relatives' "
                         f"1:1 symbol ({100 * agree / compared:.1f}%; {c.get('check_same_base', 0):,} by "
                         f"suffix-free form, {c.get('check_synonym', 0):,} as former names of the same gene); "
                         f"{c.get('check_different', 0):,} differ")
        sc = self.same_check
        if self.same and sc["compared"]:
            agree = sc["same"] + sc["same_base"] + sc["synonym"]
            lines.append(f"    Check, previous annotation:  {agree:,} of {sc['compared']:,} supported symbols agree with "
                         f"the {self.same.organism} annotation ({self.same.accession}) "
                         f"({100 * agree / sc['compared']:.1f}%); {sc['different']:,} differ")
        n = self.normalized
        if sum(n.values()):
            kept = n["no_reference_ortholog"] + n["reference_pseudogene"] + n["reference_readthrough"]
            lines.append(f"    Relatives' symbols:          {n['changed']:,} put in the reference nomenclature "
                         f"(and {n['identifier_replaced']:,} LOC identifiers), {n['same']:,} already in it, "
                         f"{kept:,} kept (no coding, non-readthrough reference ortholog)")
        lines.append(f"    Relatives disagree:          {c.get('check_relatives_disagree', 0):,} genes")
        if self.report:
            lines.append(f"    Per gene:                    reports/{Path(self.report).name}")
        return lines

    def report_files(self) -> list:
        return [self.report] if self.report else []


def _pseudogenes(db: str, geneids: set) -> set:
    """GeneIDs among geneids that NCBI gene_info types as pseudogenes."""
    import sqlite3
    if not geneids:
        return set()
    conn = sqlite3.connect(db)
    try:
        conn.execute("CREATE TEMP TABLE _pseudo_filter (gid TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO _pseudo_filter VALUES (?)", [(g,) for g in geneids])
        return {r[0] for r in conn.execute("SELECT g.gene_id FROM ncbi_geneinfo g INNER JOIN _pseudo_filter f "
                                           "ON g.gene_id = f.gid WHERE g.type = 'pseudo'")}
    finally:
        conn.close()


def _reference_names(db: str, taxid: str, symbols: set) -> dict:
    """
    {lower-case name: {"any": GeneIDs it names (symbol or synonym), "symbol": GeneIDs whose official
    symbol it is}} for the reference species, from ncbi_gene_symbols; {} without that table.
    """
    import sqlite3
    from goannotate import mappings
    if not symbols or not mappings.sqlite_has_table(db, "ncbi_gene_symbols"):
        return {}
    conn = sqlite3.connect(db)
    out = {}
    try:
        conn.execute("CREATE TEMP TABLE _name_filter (sym TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO _name_filter VALUES (?)", [(s,) for s in symbols])
        rows = conn.execute("SELECT s.symbol, s.gene_id, s.kind FROM _name_filter f INNER JOIN ncbi_gene_symbols s "
                            "ON s.tax_id = ? AND s.symbol = f.sym", (taxid,))
        for sym, gid, kind in rows:
            d = out.setdefault(sym, {"any": set(), "symbol": set()})
            if kind in ("symbol", "synonym"):
                d["any"].add(gid)
            if kind == "symbol":
                d["symbol"].add(gid)
    finally:
        conn.close()
    return out


def _ortholog_text(pid: str, g: RelGene) -> str:
    """'XP_1:Zfp558>ZNF558' when the reference nomenclature changed the symbol, 'XP_1:Avp', or 'XP_1'."""
    if not g.symbol:
        return pid
    if g.ref_symbol and g.ref_symbol != g.symbol:
        return f"{pid}:{g.symbol}>{g.ref_symbol}"
    return f"{pid}:{g.symbol}"


def _write_fasta(path: Path, seqs: dict) -> None:
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def _last_line(log: Path) -> str:
    try:
        with open(log, "rb") as fh:
            fh.seek(0, 2)
            fh.seek(max(0, fh.tell() - 2000))
            text = fh.read().decode("utf-8", "replace")
    except OSError:
        return ""
    text = re.sub(r"\x1b\[[0-9;?]*[A-Za-z]", "", text)       # colour and cursor codes
    lines = []
    for line in re.split(r"[\r\n]+", text):
        line = re.sub(r"[\u2500-\u257f]+", " ", line)         # progress-bar box characters
        line = re.sub(r"\s+", " ", line).strip()
        if re.search(r"[A-Za-z]", line):                        # skip "-----" and "=====" rules
            lines.append(line)
    return lines[-1][:120] if lines else ""
