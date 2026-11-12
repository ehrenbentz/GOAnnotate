"""
Checks of the orthology module on small inputs, without OrthoFinder or
NCBI datasets: stand-in scripts answer `datasets summary`/`download` from
fixtures and write OrthoFinder result tables; a miniature taxonomy table
gives the lineage.

Run: python3 tests/test_orthology.py
"""

import json
import os
import sqlite3
import stat
import subprocess
import sys
import tempfile
import zipfile
from argparse import Namespace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.filters import BadNameFilter, default_bad_names_path  # noqa: E402
from goannotate.models import GeneUniverse, QueryAnnotation  # noqa: E402
from goannotate.modules import ModuleContext, run_modules, select_modules  # noqa: E402
from goannotate.modules.orthology import (OrthologyModule, _last_line, agreed_symbol,  # noqa: E402
                                          parse_orthologs_by_species, parse_refseq_gff, species_label,
                                          symbol_base)

failures = 0


def check(cond, label):
    global failures
    print(("ok  " if cond else "FAIL") + " " + label)
    if not cond:
        failures += 1


def gff(genes, acc=""):
    """genes: [(gene id, symbol, [(protein id, product, length)])] -> RefSeq-style GFF3."""
    lines = ["##gff-version 3"]
    for gid, sym, prots in genes:
        geneid = GENEIDS.get((acc, gid), "")
        dbx = f"Dbxref=GeneID:{geneid},Genbank:x;" if geneid else ""
        lines.append(f"chr1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=gene-{gid};{dbx}gene={sym};gene_biotype=protein_coding")
        for i, (pid, prod, _) in enumerate(prots):
            lines.append(f"chr1\tGnomon\tmRNA\t1\t100\t.\t+\t.\tID=rna-{pid};Parent=gene-{gid};gene={sym}")
            lines.append(f"chr1\tGnomon\tCDS\t1\t100\t.\t+\t0\tID=cds-{pid};Parent=rna-{pid};{dbx}gene={sym};"
                         f"product={prod.replace(' ', '%20') if False else prod};protein_id={pid}")
    return "\n".join(lines) + "\n"


def faa(genes):
    out = []
    for _, _, prots in genes:
        for pid, prod, n in prots:
            out.append(f">{pid} {prod}\n" + "M" + "A" * (n - 1))
    return "\n".join(out) + "\n"


GENOMES = {
    # rejected: one of three genes named
    "GCF_A": ("Microtus a", 111, "2025-08-01", [("Avp", "Avp", [("XP_A1", "vasopressin", 120)]),
                                                ("LOC1", "LOC1", [("XP_A2", "uncharacterized protein LOC1", 90)]),
                                                ("LOC2", "LOC2", [("XP_A3", "uncharacterized protein LOC2", 90)])]),
    # chosen: two of three named; Kcnq5 has two isoforms
    "GCF_B": ("Microtus b", 112, "2024-01-01", [("Avp", "Avp", [("XP_B1", "vasopressin-neurophysin 2-copeptin", 150)]),
                                                ("Kcnq5", "Kcnq5", [("XP_B2", "potassium channel KCNQ5 isoform X2", 200),
                                                                    ("XP_B3", "potassium channel KCNQ5 isoform X1", 900)]),
                                                ("LOC3", "LOC3", [("XP_B4", "uncharacterized protein LOC3", 80)]),
                                                ("Zfp558", "Zfp558", [("XP_B5", "zinc finger protein 558", 300)]),
                                                ("Cdhr18", "Cdhr18", [("XP_B6", "cadherin-related family member 18",
                                                                       400)]),
                                                ("Oxt", "Oxt", [("XP_B7", "oxytocin-neurophysin 1", 125)])]),
    "GCF_C": ("Myodes c", 113, "2026-01-01", [("Avp", "Avp", [("XP_C1", "vasopressin-neurophysin 2-copeptin", 150)]),
                                              ("Kcnq5", "Kcnq5", [("XP_C2", "potassium channel KCNQ5", 900)]),
                                              ("Grm1", "Grm1", [("XP_C3", "metabotropic glutamate receptor 1", 1100)]),
                                              ("Znf558", "Znf558", [("XP_C4", "zinc finger protein 558", 300)]),
                                              ("Tnfaip8l3", "Tnfaip8l3", [("XP_C5", "TNF alpha-induced protein 8-like 3",
                                                                           200)]),
                                              ("LOC77", "LOC77", [("XP_C6", "uncharacterized protein LOC77", 200)]),
                                              ("LOC78", "LOC78", [("XP_C7", "uncharacterized protein LOC78", 150)]),
                                              ("Oxt", "Oxt", [("XP_C8", "oxytocin-neurophysin 1", 125)]),
                                              ("Cunh19orf81", "Cunh19orf81",
                                               [("XP_C9", "uncharacterized protein C19orf81 homolog", 250)])]),
    # a distant relative (another family), named with --relative
    "GCF_E": ("Mus e", 115, "2025-01-01", [("Oxt", "Oxt", [("XP_E1", "oxytocin-neurophysin 1", 125)]),
                                           ("Sstr1", "Sstr1", [("XP_E2", "somatostatin receptor type 1", 390)]),
                                           ("Avp", "Avp", [("XP_E3", "vasopressin", 150)])]),
    # half named: passes the floor, fails the ratio against Myodes c at the same level
    "GCF_D": ("Myodes d", 114, "2020-01-01", [("Avp", "Avp", [("XP_D1", "vasopressin", 150)]),
                                              ("LOC9", "LOC9", [("XP_D2", "uncharacterized protein LOC9", 100)])]),
    # the query species' own annotation (--same-species)
    "GCF_Q": ("Microtus ochrogaster", 79684, "2026-06-01", [("Avp", "Avp", [("XP_Q1", "vasopressin", 150)]),
                                                           ("Grm1", "Grm1", [("XP_Q2", "glutamate receptor", 1100)]),
                                                           ("Oldname", "Oldname", [("XP_Q3", "old protein", 200)]),
                                                           ("Sstr2", "Sstr2", [("XP_Q4", "somatostatin receptor 2",
                                                                                 300)]),
                                                           ("Zfp99", "Zfp99", [("XP_Q5", "zinc finger 99", 300)])]),
}
SUMMARY = {"10053": ["GCF_Q", "GCF_A", "GCF_B"], "39087": ["GCF_B", "GCF_C", "GCF_D"], "79684": ["GCF_Q"]}
# NCBI GeneIDs of the relatives' genes; gene_orthologs links them to human (9606) genes
GENEIDS = {("GCF_B", "Avp"): "2001", ("GCF_B", "Kcnq5"): "2002", ("GCF_B", "Zfp558"): "2005",
           ("GCF_B", "Cdhr18"): "2006", ("GCF_C", "Avp"): "3001", ("GCF_C", "Kcnq5"): "3002", ("GCF_C", "Grm1"): "3003", ("GCF_C", "Znf558"): "3004",
           ("GCF_C", "Tnfaip8l3"): "3005", ("GCF_C", "LOC77"): "3006", ("GCF_C", "LOC78"): "3007",
           ("GCF_Q", "Avp"): "4001", ("GCF_Q", "Grm1"): "4002", ("GCF_Q", "Zfp99"): "4005",
           ("GCF_B", "Oxt"): "2007", ("GCF_C", "Oxt"): "3008", ("GCF_C", "Cunh19orf81"): "3009",
           ("GCF_E", "Oxt"): "5001", ("GCF_E", "Sstr1"): "5002", ("GCF_E", "Avp"): "5003"}
HUMAN = [("9001", "AVP", [("112", "2001"), ("113", "3001"), ("79684", "4001"), ("115", "5003")]),
         ("9012", "OXT", [("112", "2007"), ("113", "3008"), ("115", "5001")]),
         ("9013", "SSTR1", [("115", "5002")]),
         ("9014", "C19orf81", [("113", "3009")]),
         ("9002", "KCNQ5", [("112", "2002"), ("113", "3002")]),
         ("9005", "ZNF558", [("112", "2005"), ("113", "3004")]),
         ("9006", "CDHR18P", [("112", "2006")]),
         ("9007", "TIPE3", [("113", "3005")]),                  # renamed; TNFAIP8L3 is a synonym
         ("9008", "KLRC4-KLRK1", [("113", "3006")]),            # readthrough
         ("9009", "KLRC4", []), ("9010", "KLRK1", []),
         ("9011", "POTEA", [("113", "3007")])]

FAKE_DATASETS = r'''#!/usr/bin/env python3
import json, os, sys, zipfile
fix = json.load(open(os.environ["FAKE_DS_FIX"]))
a = sys.argv[1:]
if a[:1] == ["--version"]:
    print("datasets version: 16.35.0"); sys.exit(0)
if a[:2] == ["summary", "genome"]:
    if os.environ.get("FAKE_DS_FAIL"):
        print("panic: runtime error: invalid memory address or nil pointer dereference"); sys.exit(2)
    if a[2] == "taxon" and a[3] in fix.get("no_data", []):
        print(f"Error: The taxonomy ID '{a[3]}' is valid for 'X', but no genome data is currently available "
              "for this taxon."); sys.exit(1)
    accs = fix["summary"].get(a[3], []) if a[2] == "taxon" else [a[3]]
    for acc in accs:
        org, taxid, rel, genes = fix["genomes"][acc]
        print(json.dumps({"accession": acc, "organism": {"organism_name": org, "tax_id": taxid},
                          "annotation_info": {"name": acc + "-RS", "release_date": rel,
                                              "stats": {"gene_counts": {"protein_coding": len(genes)}}}}))
    sys.exit(0)
if a[:2] == ["download", "genome"]:
    acc = a[3]; inc = a[a.index("--include") + 1].split(","); fn = a[a.index("--filename") + 1]
    with open(os.environ["FAKE_DS_LOG"], "a") as log: log.write(acc + " " + ",".join(inc) + "\n")
    with zipfile.ZipFile(fn, "w") as z:
        if "gff3" in inc: z.writestr(f"ncbi_dataset/data/{acc}/genomic.gff", fix["gff"][acc])
        if "protein" in inc: z.writestr(f"ncbi_dataset/data/{acc}/protein.faa", fix["faa"][acc])
    sys.exit(0)
sys.exit(1)
'''

FAKE_ORTHOFINDER = r'''#!/usr/bin/env python3
import os, shutil, sys
a = sys.argv[1:]
if a[:1] == ["-h"]:
    print("SIMPLE USAGE: Run full OrthoFinder analysis on FASTA format proteomes"); sys.exit(0)
if a[:1] == ["-v"]:
    print("\x1b[38;5;136;49mOrthoFinder\x1b[0m\x1b[39;49m:v\x1b[0m\x1b[38;5;38;49m3.\x1b[0m\x1b[1;38;5;38;49m2\x1b[0m")
    sys.exit(0)
with open(os.environ["FAKE_OF_ARGS"], "w") as fh: fh.write(" ".join(a))
indir = a[a.index("-f") + 1]; out = a[a.index("-o") + 1]; name = a[a.index("-n") + 1]
species = sorted(f[:-4] for f in os.listdir(indir) if f.endswith(".faa"))
with open(os.environ["FAKE_OF_INPUT"], "w") as fh:
    for sp in species:
        ids = [l[1:].strip() for l in open(os.path.join(indir, sp + ".faa")) if l.startswith(">")]
        fh.write(sp + "\t" + ",".join(ids) + "\n")
res = os.path.join(out, "Results_" + name)
src = os.environ["FAKE_OF_FIX"]
os.makedirs(os.path.join(res, "Orthogroups"))
shutil.copy(os.path.join(src, "Orthogroups.tsv"), os.path.join(res, "Orthogroups"))
if os.environ.get("FAKE_OF_FORMAT") == "2":           # OrthoFinder 2: one table per relative
    os.makedirs(os.path.join(res, "Orthologues", "Orthologues_query"))
    for f in os.listdir(src):
        if f.startswith("query__v__"):
            shutil.copy(os.path.join(src, f), os.path.join(res, "Orthologues", "Orthologues_query"))
else:                                                 # OrthoFinder 3: one gzipped long table, one row per ortholog
    import gzip
    os.makedirs(os.path.join(res, "Orthologs"))
    with gzip.open(os.path.join(res, "Orthologs", "query.tsv.gz"), "wt") as out_fh:
        out_fh.write("Orthogroup\tSpecies\tquery\tOrthologs\n")
        for f in sorted(os.listdir(src)):
            if f.startswith("query__v__"):
                sp = f[len("query__v__"):-4]
                for line in list(open(os.path.join(src, f)))[1:]:
                    og, q, r = line.rstrip("\n").split("\t")
                    for gene in r.split(", "):
                        out_fh.write(f"{og}\t{sp}\t{q}\t{gene}\n")
'''


def setup(tmp: Path):
    fix = {"genomes": {k: [v[0], v[1], v[2], [g[0] for g in v[3]]] for k, v in GENOMES.items()},
           "summary": SUMMARY,
           "no_data": ["999"],
           "gff": {k: gff(v[3], k) for k, v in GENOMES.items()},
           "faa": {k: faa(v[3]) for k, v in GENOMES.items()}}
    (tmp / "ds_fix.json").write_text(json.dumps(fix))
    for name, text in (("datasets", FAKE_DATASETS), ("orthofinder", FAKE_ORTHOFINDER)):
        p = tmp / name
        p.write_text(text)
        p.chmod(p.stat().st_mode | stat.S_IXUSR)
    of = tmp / "of_fix"
    of.mkdir()
    (of / "Orthogroups.tsv").write_text(
        "Orthogroup\tMicrotus_b\tMyodes_c\tquery\n"
        "OG1\tXP_B1\tXP_C1\tg1\nOG2\tXP_B3\tXP_C2\tg2, g5b\nOG3\t\tXP_C3\tg3\nOG5\tXP_B1\tXP_C2\tg5\n"
        "OG6\t\tXP_C1, XP_C3\tg6\nOG7\tXP_B5\tXP_C4\tg7\nOG8\t\tXP_C5\tg8\nOG9\tXP_B7\tXP_C8\tg9\nOG11\t\tXP_C9\tg11\n")
    (of / "query__v__Microtus_b.tsv").write_text(
        "Orthogroup\tquery\tMicrotus_b\nOG1\tg1\tXP_B1\nOG2\tg2\tXP_B3\nOG5\tg5\tXP_B1\nOG7\tg7\tXP_B5\n"
        "OG9\tg9\tXP_B7\n")
    (of / "query__v__Myodes_c.tsv").write_text(
        "Orthogroup\tquery\tMyodes_c\nOG1\tg1\tXP_C1\nOG2\tg2\tXP_C2\nOG3\tg3\tXP_C3\nOG5\tg5\tXP_C2\n"
        "OG6\tg6\tXP_C1, XP_C3\nOG7\tg7\tXP_C4\nOG8\tg8\tXP_C5\nOG9\tg9\tXP_C8\nOG11\tg11\tXP_C9\n")
    (of / "query__v__Mus_e.tsv").write_text(
        "Orthogroup\tquery\tMus_e\nOG1\tg1\tXP_E3\nOG9\tg9\tXP_E1\nOG10\tg10\tXP_E2\n")
    (of / "query__v__Microtus_ochrogaster.tsv").write_text(
        "Orthogroup\tquery\tMicrotus_ochrogaster\nOG1\tg1\tXP_Q1\nOG3\tg3\tXP_Q2\nOG5\tg5\tXP_Q3\n"
        "OG6\tg6\tXP_Q4\nOG7\tg7\tXP_Q5\n")
    db = tmp / "tax.sqlite"
    conn = sqlite3.connect(db)
    conn.execute("CREATE TABLE taxonomy_nodes (tax_id INTEGER, parent_tax_id INTEGER, rank TEXT, name TEXT)")
    conn.executemany("INSERT INTO taxonomy_nodes VALUES (?, ?, ?, ?)",
                     [(79684, 10053, "species", "Microtus ochrogaster"), (10053, 39087, "genus", "Microtus"),
                      (39087, 337677, "subfamily", "Arvicolinae"), (337677, 337687, "family", "Cricetidae"),
                      (337687, 9989, "no rank", "Eumuroida"), (9989, 1, "order", "Rodentia"), (1, 1, "no rank", "root"),
                      (111, 10053, "species", "Microtus a"), (112, 10053, "species", "Microtus b"),
                      (10060, 39087, "genus", "Myodes"), (113, 10060, "species", "Myodes c"),
                      (114, 10060, "species", "Myodes d"),
                      (10088, 39107, "genus", "Mus"), (39107, 10066, "subfamily", "Murinae"),
                      (10066, 337687, "family", "Muridae"), (115, 10088, "species", "Mus e")])
    conn.execute("CREATE TABLE ncbi_orthologs (tax_id TEXT, gene_id TEXT, relationship TEXT, other_tax_id TEXT, "
                 "other_gene_id TEXT)")
    conn.execute("CREATE TABLE ncbi_geneinfo (gene_id TEXT, symbol TEXT, type TEXT)")
    conn.execute("CREATE TABLE ncbi_gene_symbols (tax_id TEXT, symbol TEXT, gene_id TEXT, kind TEXT)")
    conn.execute("INSERT INTO ncbi_gene_symbols VALUES ('9606', 'tnfaip8l3', '9007', 'synonym')")
    for hid, sym, others in HUMAN:
        conn.execute("INSERT INTO ncbi_geneinfo VALUES (?, ?, ?)",
                     (hid, sym, "pseudo" if sym == "CDHR18P" else "protein-coding"))
        conn.execute("INSERT INTO ncbi_gene_symbols VALUES ('9606', ?, ?, 'symbol')", (sym.lower(), hid))
        conn.executemany("INSERT INTO ncbi_orthologs VALUES ('9606', ?, 'Ortholog', ?, ?)",
                         [(hid, t, g) for t, g in others])
    conn.commit()
    conn.close()
    os.environ.update({"FAKE_DS_FIX": str(tmp / "ds_fix.json"), "FAKE_DS_LOG": str(tmp / "ds_log.txt"),
                       "FAKE_OF_ARGS": str(tmp / "of_args.txt"), "FAKE_OF_INPUT": str(tmp / "of_input.txt"),
                       "FAKE_OF_FIX": str(of)})
    return db


def test_helpers():
    check(species_label("Peromyscus maniculatus bairdii") == "Peromyscus_maniculatus_bairdii", "species label")
    check(symbol_base("foxo3b") == "FOXO3" and symbol_base("dhdh.2") == "DHDH" and symbol_base("AVP") == "AVP",
          "suffix-free symbol")
    check(agreed_symbol(["Avp", "AVP"]) == ("Avp", "agree"), "agreement is case-insensitive")
    check(agreed_symbol(["foxo3b", "foxo3b", "foxo3"]) == ("foxo3b", "agree_base"), "most frequent form")
    check(agreed_symbol(["foxo3a", "foxo3b"]) == ("FOXO3", "agree_base"), "tie gives the suffix-free form")
    check(agreed_symbol(["Avp", "Kcnq5"]) == ("", "disagree") and agreed_symbol([]) == ("", "none"),
          "disagreement and no symbol")


def test_orthofinder3_table(tmp: Path):
    import gzip
    path = tmp / "query.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("Orthogroup\tSpecies\tquery\tOrthologs\n"
                 "OG1\tSp_b\tq1\tb1\n"                    # 1:1
                 "OG2\tSp_b\tq2\tb2\nOG2\tSp_b\tq2\tb3\n"   # 1:many split over two rows
                 "OG3\tSp_b\tq3, q4\tb4\n"                # many:1
                 "OG1\tSp_c\tq1\tc1, c2\n")               # 1:many in one row
    rows = parse_orthologs_by_species(path)
    b = {r[1][0]: r for r in rows["Sp_b"]}
    check(b["q1"] == ("OG1", ["q1"], ["b1"]), "OrthoFinder 3: 1:1")
    check(b["q2"] == ("OG2", ["q2"], ["b2", "b3"]), "OrthoFinder 3: orthologs split over rows regrouped")
    check(b["q3"] == ("OG3", ["q3", "q4"], ["b4"]) and b["q4"] == ("OG3", ["q4", "q3"], ["b4"]),
          "OrthoFinder 3: co-orthologs in the query, the gene itself first")
    check(rows["Sp_c"] == [("OG1", ["q1"], ["c1", "c2"])], "OrthoFinder 3: species kept apart")


def test_progress_line(tmp: Path):
    log = tmp / "orthofinder.log"
    log.write_text("\x1b[1;36;49m2026\x1b[0m-11-01 17:40:14 : Scheduling species tree\n"
                   "Processing... \u2501\u2501\u2501   16/16 0:03:08\r\n"
                   "-------------------------------\n=====\n", encoding="utf-8")
    check(_last_line(log) == "Processing... 16/16 0:03:08", "progress line without colour codes and rules")
    log.write_text("\x1b[1;36;49m2026\x1b[0m-11-01 17:40:14 : Scheduling species tree\n------\n",
                   encoding="utf-8")
    check(_last_line(log) == "2026-11-01 17:40:14 : Scheduling species tree", "colour codes removed")
    check(_last_line(tmp / "missing.log") == "", "missing log")


def test_gff(tmp: Path):
    p = tmp / "b.gff"
    p.write_text(gff(GENOMES["GCF_B"][3], "GCF_B"))
    info = parse_refseq_gff(p)
    check(info["XP_B3"] == ("gene-Kcnq5", "Kcnq5", "potassium channel KCNQ5 isoform X1", "protein_coding", "2002"),
          f"CDS traced to its gene and GeneID: {info['XP_B3']}")


def test_module(tmp: Path, db: Path):
    proteins = tmp / "q.faa"
    names = ["g1", "g2", "g3", "g5", "g6", "g7", "g8", "g9", "g10", "g11"]
    proteins.write_text("".join(f">{g}_t1\nM{'A' * 99}\n" for g in names) + ">g2_t2\nM" + "A" * 299 + "\n")
    tx2g = {f"{g}_t1": g for g in names}
    tx2g["g2_t2"] = "g2"
    g2t = {}
    for t, g in tx2g.items():
        g2t.setdefault(g, []).append(t)
    universe = GeneUniverse(gene_ids=names, transcript_to_gene=tx2g, gene_to_transcripts=g2t,
                            coding={g: True for g in names})
    ann = {g: QueryAnnotation(query_id=g) for g in names}
    ann["g2"].consensus_symbol, ann["g2"].symbol_status = "kcnq5b", "foreign_kept"
    ann["g2"].consensus_name = "potassium voltage-gated channel subfamily KQT member 5"
    ann["g3"].consensus_symbol, ann["g3"].symbol_status = "GRM1", "vote_supported"
    ann["g8"].consensus_symbol, ann["g8"].symbol_status = "TNFAIP8L3", "rbh"
    bf = BadNameFilter(str(default_bad_names_path()))
    os.chdir(tmp)
    work, reports = Path("intermediate"), Path("reports")
    work.mkdir()
    reports.mkdir()
    args = Namespace(modules="orthology", orthofinder_bin=str(tmp / "orthofinder"), datasets_bin=str(tmp / "datasets"),
                     relatives=2, relative=[], relative_min_named=0.2, relative_named_ratio=0.75,
                     relative_named_margin=0.1, relatives_dir=None, orthofinder_threads=48,
                     ortholog_symbol="fill", taxid=79684, db=str(db), force_modules=False, same_species=False,
                     reference_taxid="9606", ortholog_anchor_taxid="9606", no_resolve_orthologs=False)
    mods = select_modules(args)
    check(len(mods) == 1 and isinstance(mods[0], OrthologyModule), "orthology module selected")
    ctx = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                        prefix="p", threads=60, universe=universe, annotations=ann, bad_name_filter=bf,
                        symbol_case_fn=str.upper)
    ran = run_modules(mods, ctx)
    check(len(ran) == 1, "module ran with the stand-in tools")
    m = ran[0]
    check(m.status.version == "OrthoFinder 3.2, datasets 16.35.0", f"versions: {m.status.version}")
    check([(r.accession, r.status) for r in m.candidates] == [("GCF_A", "rejected_ratio"), ("GCF_B", "chosen"),
                                                            ("GCF_C", "chosen"), ("GCF_D", "rejected_ratio")],
          f"candidates: query species skipped, A and D below the ratio floor: "
          f"{[(r.accession, r.status) for r in m.candidates]}")
    check(abs(m.candidates[2].named - 7 / 9) < 1e-9 and abs(m.candidates[2].named_symbol - 6 / 9) < 1e-9,
          f"informative fraction counts LOC78 through its reference ortholog POTEA: {m.candidates[2].named}")
    check([r.level for r in m.relatives] == ["Microtus", "Arvicolinae"], "levels recorded")
    check([(r.close, r.lca_rank, r.lca_name) for r in m.relatives] == [(True, "genus", "Microtus"),
                                                                       (True, "subfamily", "Arvicolinae")],
          f"closeness from the taxonomy: {[(r.close, r.lca_rank, r.lca_name) for r in m.relatives]}")
    check(abs(m.candidates[0].named - 1 / 3) < 1e-9 and abs(m.relatives[0].named - 5 / 6) < 1e-9, "named fractions")
    check([r.level_index for r in m.relatives] == [0, 1], "level indices")
    dl = (tmp / "ds_log.txt").read_text().split("\n")
    check("GCF_D gff3" in dl and "GCF_D gff3,protein" not in dl and "GCF_B gff3,protein" in dl
          and "GCF_A gff3,protein" not in dl, f"every candidate's GFF3, only chosen relatives' proteins: {dl}")
    of_args = (tmp / "of_args.txt").read_text().split()
    check(of_args[of_args.index("-t") + 1] == "48" and of_args[of_args.index("-a") + 1] == "12"
          and of_args[of_args.index("-f") + 1].startswith("/"), f"OrthoFinder threads and absolute paths: {of_args}")
    inputs = dict(l.split("\t") for l in (tmp / "of_input.txt").read_text().splitlines())
    check(inputs["Microtus_b"] == "XP_B1,XP_B3,XP_B4,XP_B5,XP_B6,XP_B7" and inputs["query"] == "g1,g2,g3,g5,g6,g7,g8,g9,g10,g11",
          f"one protein per gene, the longest: {inputs}")
    # Symbols
    check(ann["g1"].consensus_symbol == "AVP" and ann["g1"].symbol_status == "ortholog_fill"
          and ann["g1"].consensus_name == "vasopressin-neurophysin 2-copeptin" and ann["g1"].product_note == "ortholog",
          f"unnamed gene filled: {ann['g1'].consensus_symbol} / {ann['g1'].consensus_name}")
    check(ann["g2"].consensus_symbol == "KCNQ5" and "kcnq5b" in ann["g2"].symbol_aliases
          and ann["g2"].consensus_name == "potassium voltage-gated channel subfamily KQT member 5",
          "lineage symbol replaced, kept as alias; product kept")
    check(ann["g3"].consensus_symbol == "GRM1" and m.check["compared"] == 2 and m.check["same"] == 1,
          "supported symbol checked, not changed")
    check(ann["g5"].consensus_symbol == "" and m.check["relatives_disagree"] == 1, "disagreeing relatives fill nothing")
    check(ann["g6"].consensus_symbol == "" and m.gene_symbol["g6"] == ("", "none"), "1:many gives no symbol")
    check(ann["g7"].consensus_symbol == "ZNF558" and m.gene_symbol["g7"] == ("ZNF558", "agree")
          and ann["g7"].consensus_name == "zinc finger protein 558",
          f"mouse-style Zfp558 and Znf558 agree in the reference nomenclature: {m.gene_symbol['g7']}")
    check(m.normalized == {"same": 7, "changed": 2, "no_reference_ortholog": 1, "reference_pseudogene": 1,
                           "identifier_replaced": 2}, f"normalization counts {m.normalized}")
    check(m.prior() == {} and m.gene_symbol["g9"] == ("Oxt", "agree") and m.gene_close["g9"], "fill mode gives no prior")
    pc = m.relatives[1].proteins
    check(pc["XP_C5"].ref_symbol == "TIPE3", "a former name is updated to the reference gene's symbol")
    check(m.gene_symbol["g8"] == ("TIPE3", "agree") and ann["g8"].consensus_symbol == "TNFAIP8L3"
          and m.check["synonym"] == 1,
          f"the core's former name of the same gene agrees (TNFAIP8L3 / TIPE3): {m.check}")
    check(pc["XP_C6"].ref_symbol == "LOC77", "a readthrough reference ortholog is not used")
    check(pc["XP_C7"].ref_symbol == "POTEA", "a LOC identifier takes its reference ortholog's symbol")
    check(m.relatives[0].proteins["XP_B6"].ref_symbol == "Cdhr18", "a pseudogene reference ortholog is not used")
    check(ann["g9"].consensus_symbol == "OXT" and ann["g9"].symbol_status == "ortholog_fill", "g9 filled")
    check(ann["g11"].consensus_symbol == "C19ORF81" and m.gene_symbol["g11"] == ("C19orf81", "agree"),
          f"an official C#orf symbol is accepted from a close relative's ortholog: {m.gene_symbol.get('g11')}")
    check(m.counts["symbols_filled"] == 5 and m.counts["filled_unnamed"] == 4 and m.counts["filled_foreign_kept"] == 1
          and m.counts["products_filled"] == 3 and m.counts["symbols_taken"] == 0, f"counts {m.counts}")
    rel_rows = (reports / "p_orthology_relatives.tsv").read_text().splitlines()
    check(len(rel_rows) == 5 and rel_rows[2].split("\t")[7:13] == ["0.833", "0.833", "chosen", "genus", "Microtus", "close"]
          and rel_rows[3].split("\t")[7:10] == ["0.778", "0.667", "chosen"] and rel_rows[4].split("\t")[9] == "rejected_ratio",
          f"relatives report: {rel_rows[2].split(chr(9))[7:]} / {rel_rows[3].split(chr(9))[7:10]}")
    rep = {r.split("\t")[0]: r.split("\t") for r in (reports / "p_orthology.tsv").read_text().splitlines()[1:]}
    check(rep["g2"][4:9] == ["Kcnq5", "agree", "close", "filled", "kcnq5b"] and rep["g1"][8] == "(none)"
          and rep["g3"][7] == "confirmed", f"per-gene report: {rep['g2'][4:9]}")
    check(rep["g7"][11] == "XP_B5:Zfp558>ZNF558" and rep["g7"][13] == "XP_C4:Znf558",
          f"report shows the relative's symbol and its reference form: {rep['g7'][11:]}")
    check(any("Relatives' symbols" in l and "2 put in the reference nomenclature" in l for l in m.summary_lines()),
          "summary normalization line")
    check(m._orthologs_text("g6") == "Myodes_c:Avp,Grm1(1:many)", f"orthologs column: {m._orthologs_text('g6')}")

    # A relative's gene shared by several query genes (many:1): used with transcript input only
    saved = dict(m.orthologs["g9"])
    m.orthologs["g9"] = {label: ("many:1", named) for label, (_, named) in saved.items()}
    m.gene_symbol.pop("g9")
    m._compute_calls(lambda s: s.startswith("LOC"))
    check(m.gene_symbol["g9"] == ("", "none") and not m._usable("many:1"),
          f"gene models: a many:1 ortholog gives no symbol: {m.gene_symbol['g9']}")
    m.fragment_orthologs = True
    m._compute_calls(lambda s: s.startswith("LOC"))
    check(m.gene_symbol["g9"] == ("Oxt", "agree") and m.shared_calls == {"g9"} and not m._usable("1:many"),
          f"transcript input: a many:1 ortholog is used like 1:1: {m.gene_symbol['g9']}")
    bare = QueryAnnotation(query_id="g9")
    voted = QueryAnnotation(query_id="g9")
    voted.voted_symbol = "OXT"
    ctx_s = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                          prefix="p", threads=60, reference_symbols={"AVP", "OXT"},
                          reference_scores={"g9": {"AVP": 300.0, "OXT": 280.0}})
    check(not m._own_evidence(ctx, bare, "g9", "Oxt") and m._own_evidence(ctx, voted, "g9", "Oxt")
          and m._own_evidence(ctx_s, bare, "g9", "Oxt") and not m._own_evidence(ctx_s, bare, "g1", "Oxt"),
          "a shared ortholog needs the transcript's own vote, hits or reference alignment")
    m.fragment_orthologs = False
    m.shared_calls.clear()
    m.orthologs["g9"] = saved
    check(OrthologyModule(Namespace(**{**vars(args), "transcripts": "t.fasta"})).fragment_orthologs
          and not OrthologyModule(args).fragment_orthologs, "many:1 orthologs follow the input type")

    # A fill is refused when the gene's sequence scores the symbol's reference protein below the margin
    ctx_f = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                          prefix="p", threads=60, reference_symbols={"AVP", "OXT"},
                          reference_scores={"g1": {"AVP": 300.0, "OXT": 200.0}, "g9": {"AVP": 300.0, "OXT": 280.0}})
    check(m._sequence_allows_fill(ctx_f, "g1", "AVP") and not m._sequence_allows_fill(ctx_f, "g1", "Oxt")
          and m._sequence_allows_fill(ctx_f, "g9", "Oxt"), "fill allowed at 90% of the best reference hit or more")
    check(m._sequence_allows_fill(ctx_f, "g1", "GAS8") and m._sequence_allows_fill(ctx_f, "g2", "AVP")
          and m._sequence_allows_fill(ctx_f, "g1", "Kcnq5"),
          "a lineage symbol and a gene without reference hits cannot be judged and are filled")
    ctx_f.reference_symbols.add("GAS8")
    check(not m._sequence_allows_fill(ctx_f, "g1", "Gas8"), "a reference symbol the gene does not hit is refused")
    check(not m._sequence_allows_fill(ctx_f, "g1", "DRC4", "Gas8") and m._sequence_allows_fill(ctx_f, "g1", "DRC4"),
          "the relatives' own form is judged when the reference proteome names the gene by it")
    ann_f = {g: QueryAnnotation(query_id=g) for g in names}
    ann_f["g2"].consensus_symbol, ann_f["g2"].symbol_status = "kcnq5b", "foreign_kept"
    ctx_r = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                          prefix="r", threads=60, universe=universe, annotations=ann_f, bad_name_filter=bf,
                          symbol_case_fn=str.upper, reference_symbols={"AVP", "OXT", "KCNQ5"},
                          reference_scores={"g1": {"OXT": 300.0}, "g2": {"AVP": 300.0, "KCNQ5": 100.0},
                                            "g9": {"OXT": 250.0}})
    m_r = OrthologyModule(args)
    run_modules([m_r], ctx_r)
    check(ann_f["g1"].consensus_symbol == "" and m_r.decision["g1"] == "fill_refused"
          and ann_f["g2"].consensus_symbol == "kcnq5b" and m_r.decision["g2"] == "fill_refused"
          and ann_f["g9"].consensus_symbol == "OXT" and m_r.counts["check_fill_refused"] == 2,
          f"fills refused by the sequence: g1 {m_r.decision['g1']}, g2 {m_r.decision['g2']}, g9 {m_r.decision['g9']}")
    check(any("Fills refused:               2" in l for l in m_r.summary_lines()), "summary line of refused fills")
    ann_v = {g: QueryAnnotation(query_id=g) for g in names}
    ann_v["g1"].voted_symbol, ann_v["g1"].symbol_status = "AVP", "unsupported_removed"
    ctx_r.annotations, ctx_r.prefix = ann_v, "v"
    m_v = OrthologyModule(args)
    run_modules([m_v], ctx_r)
    check(ann_v["g1"].consensus_symbol == "AVP" and m_v.decision["g1"] != "fill_refused",
          f"a symbol the gene's own vote carried is not refused: {m_v.decision['g1']}")
    check(any("Check, supported symbols" in l and "2 of 2 agree" in l for l in m.summary_lines()), "summary check line")
    check(any("2 close (common ancestor within the family)" in l for l in m.summary_lines()), "summary closeness")

    # Reuse: same relatives and query, OrthoFinder not rerun, downloads reused
    (tmp / "of_args.txt").unlink()
    (tmp / "ds_log.txt").write_text("")
    m2 = OrthologyModule(args)
    m2.detect()
    m2.run(ctx)
    check(not (tmp / "of_args.txt").exists() and (tmp / "ds_log.txt").read_text() == "",
          "OrthoFinder results and downloads reused")
    check(m2.orthologs == m.orthologs, "reused results read the same")

    # OrthoFinder 2 layout gives the same orthologs as OrthoFinder 3
    os.environ["FAKE_OF_FORMAT"] = "2"
    v2 = tmp / "of2"
    subprocess.run([str(tmp / "orthofinder"), "-f", str(work / "p_orthology" / "input"), "-o", str(v2),
                    "-n", "goannotate"], check=True)
    del os.environ["FAKE_OF_FORMAT"]
    m5 = OrthologyModule(args)
    m5.relatives = m.relatives
    m5._read_results(v2 / "Results_goannotate")
    check(m5.orthologs == m.orthologs and m5.orthologs, "OrthoFinder 2 and 3 tables agree")

    # Manual relative by accession
    m3 = OrthologyModule(Namespace(**{**vars(args), "relative": ["GCF_D"], "taxid": None}))
    s = m3.detect()
    check(s.available, "manual relative needs no --taxid")
    m3._choose_relatives(ctx, work / "p_orthology" / "relatives", bf.is_bad_symbol)
    check([r.accession for r in m3.relatives] == ["GCF_D"] and m3.relatives[0].level == "--relative",
          "manual relative used")
    # No taxid and no relative
    m4 = OrthologyModule(Namespace(**{**vars(args), "taxid": None}))
    check(not m4.detect().available, "automatic choice needs --taxid")
    m4b = OrthologyModule(Namespace(**{**vars(args), "taxid": None, "relative": ["GCF_D"], "same_species": True}))
    check(not m4b.detect().available, "--same-species needs --taxid")

    # Without the reference tables the relatives' own symbols are used
    m7 = OrthologyModule(Namespace(**{**vars(args), "no_resolve_orthologs": True}))
    m7.relatives = m.relatives
    for r in m7.relatives:
        for g in r.proteins.values():
            g.ref_symbol = g.symbol
    m7._normalize_symbols()
    check(not m7.normalized and m.relatives[0].proteins["XP_B5"].ref_symbol == "Zfp558",
          "--no-resolve-orthologs keeps the relatives' symbols")

    # --same-species: the query species' own annotation is one of the --relatives, chosen first
    ann2 = {g: QueryAnnotation(query_id=g) for g in names}
    ann2["g3"].consensus_symbol, ann2["g3"].symbol_status = "GRM1", "vote_supported"
    ctx2 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="s", threads=60, universe=universe, annotations=ann2, bad_name_filter=bf,
                         symbol_case_fn=str.upper)
    m5 = OrthologyModule(Namespace(**{**vars(args), "same_species": True}))
    m5._choose_relatives(ctx2, work / "rel_cache_m5", lambda s: s.startswith("LOC"))
    check([r.accession for r in m5.relatives] == ["GCF_Q", "GCF_B"],
          f"same species takes one of --relatives 2: {[r.accession for r in m5.relatives]}")
    m6 = OrthologyModule(Namespace(**{**vars(args), "same_species": True, "relatives": 3}))
    ran = run_modules([m6], ctx2)
    check(len(ran) == 1 and m6.same is not None and m6.same.accession == "GCF_Q"
          and [r.accession for r in m6.relatives] == ["GCF_Q", "GCF_B", "GCF_C"],
          f"same-species annotation first among --relatives 3: {[r.accession for r in m6.relatives]}")
    inputs = dict(l.split("\t") for l in (tmp / "of_input.txt").read_text().splitlines())
    check("Microtus_ochrogaster" in inputs, "same-species proteome given to OrthoFinder")
    check(m6.gene_symbol["g1"] == ("Avp", "agree_same_species"), f"agrees with the relatives: {m6.gene_symbol['g1']}")
    check(m6.gene_symbol["g5"] == ("Oldname", "same_species") and ann2["g5"].consensus_symbol == "OLDNAME"
          and ann2["g5"].consensus_name == "old protein",
          f"used when the other relatives disagree: {m6.gene_symbol['g5']}")
    check(m6.gene_symbol["g6"] == ("Sstr2", "same_species") and ann2["g6"].consensus_symbol == "SSTR2",
          f"used when the other relatives have no 1:1 ortholog: {m6.gene_symbol['g6']}")
    check(m6.gene_symbol["g7"] == ("", "same_species_disagree") and ann2["g7"].consensus_symbol == "",
          f"not used when the other relatives agree on another symbol: {m6.gene_symbol['g7']}")
    check(m6.same_check == {"compared": 1, "same": 1} and m6.check["compared"] == 1,
          f"continuity check of supported symbols: {m6.same_check}")
    check(any("Check, previous annotation" in l and "1 of 1" in l for l in m6.summary_lines())
          and any("(same species)" in l for l in m6.summary_lines()), "same-species summary lines")
    failures_before = m6.summary_failures
    check(m6._summary("999") == [] and m6.summary_failures == failures_before,
          "a taxon NCBI reports without genome data is an empty answer, not a failed query")
    rel_rows = (reports / "s_orthology_relatives.tsv").read_text().splitlines()
    check(rel_rows[1].split("\t")[3] == "same species" and rel_rows[1].split("\t")[9:13]
          == ["chosen", "species", "Microtus ochrogaster", "close"], f"same species first in the relatives report: {rel_rows[1]}")

    # Weighted mode with a close (GCF_B, genus) and a distant (GCF_E, another family) relative
    from goannotate.modules import integrate_modules, start_modules
    ann3 = {g: QueryAnnotation(query_id=g) for g in names}
    ann3["g2"].consensus_symbol, ann3["g2"].symbol_status = "kcnq5b", "foreign_kept"
    ann3["g9"].consensus_symbol, ann3["g9"].symbol_status = "AVP", "rbh"
    ann3["g9"].consensus_name = "Vasopressin-neurophysin 2-copeptin"
    ann3["g9"].symbol_aliases = ["Oxt"]          # a hit carried the relatives' symbol
    ann3["g10"].consensus_symbol, ann3["g10"].symbol_status = "SSTR2", "rbh"
    ctx3 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="w", threads=60, universe=universe, annotations={}, bad_name_filter=bf,
                         symbol_case_fn=str.upper)
    m8 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted", "relative": ["GCF_B", "GCF_E"],
                                      "ortholog_close_rank": "family"}))
    started = start_modules([m8], ctx3)
    check(started == [m8] and m8.runs_before_vote, "weighted: OrthoFinder run before the vote")
    check([(r.accession, r.close, r.lca_rank) for r in m8.relatives] == [("GCF_B", True, "genus"), ("GCF_E", False, "order")],
          f"close and distant relatives: {[(r.accession, r.close, r.lca_rank, r.lca_name) for r in m8.relatives]}")
    prior = m8.prior()
    check(prior.get("g9") == ("Oxt", True) and prior.get("g10") == ("Sstr1", False) and prior.get("g1") == ("Avp", True),
          f"prior for the vote: {prior}")
    ctx3.annotations = ann3
    ctx3.reference_descriptions = {"OXT": "Oxytocin-neurophysin 1", "AVP": "Vasopressin-neurophysin 2-copeptin"}
    done = integrate_modules(started, ctx3)
    check(done == [m8], "weighted: integrated after the symbol assignment")
    check(ann3["g9"].consensus_symbol == "OXT" and ann3["g9"].symbol_status == "ortholog"
          and ann3["g9"].voted_symbol == "AVP" and "AVP" not in ann3["g9"].symbol_aliases
          and ann3["g9"].consensus_name == "Oxytocin-neurophysin 1" and ann3["g9"].product_note == "ortholog",
          f"close ortholog symbol taken over a supported core symbol: {ann3['g9'].consensus_symbol}, "
          f"{ann3['g9'].symbol_status}, {ann3['g9'].consensus_name}")
    check(ann3["g10"].consensus_symbol == "SSTR2" and m8.decision["g10"] == "kept" and m8.check["different"] == 2,
          f"distant ortholog symbol does not replace a supported symbol: {m8.decision['g10']}")
    check(ann3["g1"].consensus_symbol == "AVP" and ann3["g1"].symbol_status == "ortholog_fill"
          and ann3["g2"].consensus_symbol == "KCNQ5" and ann3["g2"].symbol_status == "ortholog_fill",
          "unnamed and lineage-symbol genes are filled, not marked as taken")
    check(ann3["g9"].confidence == "medium" and 0.55 <= ann3["g9"].confidence_score <= 0.69
          and ann3["g1"].confidence == "medium" and ann3["g1"].confidence_score == 0.55,
          f"confidence from orthology: taken {ann3['g9'].confidence_score}, close fill without hits "
          f"{ann3['g1'].confidence_score}")
    check(m8.counts["symbols_taken"] == 1 and m8.counts["symbols_taken_conflicts"] == 1
          and m8.counts["symbols_filled"] == 4 and m8.counts["close_relatives"] == 1, f"weighted counts: {m8.counts}")
    rep = {r.split("\t")[0]: r.split("\t") for r in (reports / "w_orthology.tsv").read_text().splitlines()[1:]}
    check(rep["g9"][4:9] == ["Oxt", "agree", "close", "taken", "AVP"] and rep["g10"][6:8] == ["distant", "kept"],
          f"decisions in the report: {rep['g9'][4:9]} {rep['g10'][4:9]}")
    check(any("Symbols taken (close):       1, of which 1 replace" in l for l in m8.summary_lines()),
          "summary line for taken symbols")
    m9 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "off", "relative": ["GCF_B", "GCF_E"]}))
    check(m9.prior() == {}, "off: no prior")

    # The outvote guard: one close relative alone, a reciprocal best hit, and no hit carrying the symbol
    ann4 = {g: QueryAnnotation(query_id=g) for g in names}
    ann4["g9"].consensus_symbol, ann4["g9"].symbol_status = "AVP", "rbh"
    ann4["g1"].consensus_symbol, ann4["g1"].symbol_status = "OXT", "rbh"     # Avp from B and E (E distant): 1 close
    ann4["g1"].hits = [type("H", (), {"gene_symbol": "Avp"})()]            # but a hit carries Avp
    ctx4 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="k", threads=60, universe=universe, annotations={}, bad_name_filter=bf,
                         symbol_case_fn=str.upper)
    m10 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted", "relative": ["GCF_B", "GCF_E"]}))
    start_modules([m10], ctx4)
    ctx4.annotations = ann4
    integrate_modules([m10], ctx4)
    check(ann4["g9"].consensus_symbol == "AVP" and m10.decision["g9"] == "conflict_kept"
          and m10.check["conflict_kept"] == 1, f"strongly outvoted: core symbol kept: {m10.decision['g9']}")
    check(ann4["g1"].consensus_symbol == "AVP" and m10.decision["g1"] == "taken",
          f"a hit carrying the relatives' symbol lifts the guard: {m10.decision['g1']}")
    check(any("Conflicts kept:              1" in l for l in m10.summary_lines()), "summary line for conflicts kept")

    # The gene's own sequence decides a taken symbol: the relatives' reference symbol replaces a supported
    # core symbol only when its reference protein scores at least 90% of the core symbol's
    for scores, expect, label in (({"g1": {"OXT": 300.0, "AVP": 280.0}}, "taken", "within 90%"),
                                  ({"g1": {"OXT": 300.0, "AVP": 200.0}}, "conflict_kept", "below 90%"),
                                  ({"g1": {"OXT": 300.0}}, "conflict_kept", "absent from the reference hits"),
                                  ({}, "taken", "no reference hit")):
        ann7 = {g: QueryAnnotation(query_id=g) for g in names}
        ann7["g1"].consensus_symbol, ann7["g1"].symbol_status = "OXT", "rbh"
        ann7["g1"].consensus_name = "Oxytocin-neurophysin 1"
        ann7["g1"].hits = [type("H", (), {"gene_symbol": "Avp"})()]
        ctx7 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                             prefix="s", threads=60, universe=universe, annotations={}, bad_name_filter=bf,
                             symbol_case_fn=str.upper, reference_symbols={"OXT", "AVP"}, reference_scores=scores,
                             reference_descriptions={"AVP": "vasopressin-neurophysin 2-copeptin"})
        m15 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted", "relative": ["GCF_B", "GCF_E"]}))
        start_modules([m15], ctx7)
        ctx7.annotations = ann7
        integrate_modules([m15], ctx7)
        a = ann7["g1"]
        product_ok = (a.consensus_name == "vasopressin-neurophysin 2-copeptin") if expect == "taken" \
            else (a.consensus_name == "Oxytocin-neurophysin 1")
        check(m15.decision["g1"] == expect and a.consensus_symbol == ("AVP" if expect == "taken" else "OXT") and product_ok,
              f"relatives' reference symbol {label}: {m15.decision['g1']}, {a.consensus_symbol}, {a.consensus_name}")
    # a lineage symbol of the relatives has no reference protein and is taken
    ann8 = {g: QueryAnnotation(query_id=g) for g in names}
    ann8["g1"].consensus_symbol, ann8["g1"].symbol_status = "OXT", "rbh"
    ann8["g1"].hits = [type("H", (), {"gene_symbol": "Avp"})()]
    ctx8 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="s", threads=60, universe=universe, annotations={}, bad_name_filter=bf,
                         symbol_case_fn=str.upper, reference_symbols={"OXT"}, reference_scores={"g1": {"OXT": 300.0}})
    m16 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted", "relative": ["GCF_B", "GCF_E"]}))
    start_modules([m16], ctx8)
    ctx8.annotations = ann8
    integrate_modules([m16], ctx8)
    check(m16.decision["g1"] == "taken", f"a symbol outside the reference nomenclature is taken: {m16.decision['g1']}")

    # The close relatives' own form of the same gene replaces the core's reference-style symbol
    from goannotate.modules.orthology import RelGene
    m12 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted"}))
    rel = type("R", (), {})()
    rel.label, rel.close, rel.same_species = "Poecilia_r", True, False
    rel.proteins = {"XP_P1": RelGene("g-kdm2aa", "kdm2aa", "lysine-specific demethylase 2Aa", "", "kdm2aa")}
    m12.relatives = [rel]
    m12.orthologs = {"q1": {"Poecilia_r": ("1:1", [("XP_P1", "kdm2aa")])}}
    m12._compute_calls(bf.is_bad_symbol)
    annq = {"q1": QueryAnnotation(query_id="q1", consensus_symbol="kdm2a", symbol_status="rbh",
                                  consensus_name="Lysine-specific demethylase 2A", hits=[])}
    ctxq = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="q", threads=1, universe=GeneUniverse(gene_ids=["q1"], transcript_to_gene={},
                                                                      gene_to_transcripts={}, coding={"q1": True}),
                         annotations=annq, bad_name_filter=bf, symbol_case_fn=str.lower)
    m12.integrate(ctxq)
    check(annq["q1"].consensus_symbol == "kdm2aa" and m12.decision["q1"] == "refined"
          and annq["q1"].voted_symbol == "kdm2a" and "kdm2a" in annq["q1"].symbol_aliases
          and annq["q1"].consensus_name == "Lysine-specific demethylase 2A" and m12.check["same_base"] == 1,
          f"ohnolog form from close relatives used, product kept: {annq['q1'].consensus_symbol} {m12.decision['q1']}")

    # A C#orf ortholog symbol fills only; it never replaces a supported symbol
    ann5 = {g: QueryAnnotation(query_id=g) for g in names}
    ann5["g11"].consensus_symbol, ann5["g11"].symbol_status = "ISOC2", "rbh"
    ctx5 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                         prefix="o", threads=60, universe=universe, annotations={}, bad_name_filter=bf,
                         symbol_case_fn=str.upper)
    m11 = OrthologyModule(Namespace(**{**vars(args), "ortholog_symbol": "weighted"}))
    start_modules([m11], ctx5)
    ctx5.annotations = ann5
    integrate_modules([m11], ctx5)
    check(ann5["g11"].consensus_symbol == "ISOC2" and m11.decision["g11"] == "kept",
          f"C#orf symbol does not replace a supported symbol: {m11.decision.get('g11')}")

    # NCBI unreachable: the relatives chosen in the previous run in this output directory are used
    os.environ["FAKE_DS_FAIL"] = "1"
    try:
        ann6 = {g: QueryAnnotation(query_id=g) for g in names}
        ctx6 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                             prefix="p", threads=60, universe=universe, annotations=ann6, bad_name_filter=bf,
                             symbol_case_fn=str.upper)
        m13 = OrthologyModule(args)
        ran6 = run_modules([m13], ctx6)
        check(len(ran6) == 1 and m13.previous_selection and m13.summary_failures > 0
              and [r.accession for r in m13.relatives] == ["GCF_B", "GCF_C"]
              and [r.level for r in m13.relatives] == ["Microtus", "Arvicolinae"]
              and all(r.close for r in m13.relatives),
              f"previous selection used when NCBI fails: {[r.accession for r in m13.relatives]}")
        ctx7 = ModuleContext(query_fasta=Path("q.faa"), search_mode="blastp", work_dir=work, reports_dir=reports,
                             prefix="none", threads=60, universe=universe, annotations=ann6, bad_name_filter=bf,
                             symbol_case_fn=str.upper)
        m14 = OrthologyModule(args)
        ran7 = run_modules([m14], ctx7)
        check(not ran7 and not m14.previous_selection, "without an earlier selection the module is skipped")
    finally:
        del os.environ["FAKE_DS_FAIL"]


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        db = setup(tmp)
        test_helpers()
        test_progress_line(tmp)
        test_orthofinder3_table(tmp)
        test_gff(tmp)
        test_module(tmp, db)
    print()
    if failures:
        print(f"{failures} check(s) failed")
        sys.exit(1)
    print("All checks passed")
