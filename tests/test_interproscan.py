"""
Checks of the InterProScan module on small inputs, without Nextflow or
InterProScan: a stand-in nextflow script writes an InterProScan TSV, and a
miniature InterPro data directory and GO hierarchy are built here.

Run: python3 tests/test_interproscan.py
"""

import json
import os
import stat
import sys
import tempfile
from argparse import Namespace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import output  # noqa: E402
from goannotate.go import GOHierarchy  # noqa: E402
from goannotate.gff3 import Gene, Transcript  # noqa: E402
from goannotate.models import GeneUniverse, QueryAnnotation  # noqa: E402
from goannotate.modules import ModuleContext, run_modules, select_modules  # noqa: E402
from goannotate.modules.interproscan import (InterProScanModule, domain_label, domain_product,  # noqa: E402
                                             parse_applications_config, parse_tsv,
                                             resources_config, RESOURCES_CONFIG)

failures = 0


def check(cond, label):
    global failures
    print(("ok  " if cond else "FAIL") + " " + label)
    if not cond:
        failures += 1


OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0007154
name: cell communication
namespace: biological_process
is_a: GO:0008150

[Term]
id: GO:0007165
name: signal transduction
namespace: biological_process
is_a: GO:0007154

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005509
name: calcium ion binding
namespace: molecular_function
is_a: GO:0003674

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0016020
name: membrane
namespace: cellular_component
is_a: GO:0005575
"""

APPS_CONFIG = """params {
    appsConfig {
        antifam {
            name = "AntiFam"
            dir = "antifam"
            hmm = "AntiFam.hmm"
            has_data = true
        }
        coils {
            name = "COILS"
            has_data = false
        }
        pfam {
            name = "Pfam"
            dir = "pfam"
            hmm = "pfam_a.hmm"
            has_data = true
        }
        prositeprofiles {
            name = "PROSITE profiles"
            dir = "prosite"
            has_data = true
        }
    }
}
"""

# Columns: protein md5 length analysis signature description start end score status date IPR IPR-desc GO pathways
TSV = "\n".join([
    "t1\tx\t300\tSMART\tSM00179\tEGF_CA\t40\t80\t1e-10\tT\t01-11-2026\tIPR001881\tEGF-like calcium-binding domain\tGO:0005509(InterPro)\t-",
    "t1\tx\t300\tSMART\tSM00051\tDSL\t5\t30\t1e-9\tT\t01-11-2026\tIPR001774\tDSL protein\tGO:0007154(InterPro)|GO:0016020(InterPro)\t-",
    "t1\tx\t300\tPfam\tPF00008\tEGF\t40\t75\t1e-8\tT\t01-11-2026\tIPR000742\tEGF-like domain\t-\t-",
    "t2\tx\t400\tPfam\tPF00069\tPkinase\t150\t380\t1e-50\tT\t01-11-2026\tIPR000719\tProtein kinase domain\tGO:0005509(InterPro)\t-",
    "t2\tx\t400\tCDD\tcd00054\tEGF_CA\t20\t60\t1e-5\tT\t01-11-2026\tIPR000742\tEGF-like domain\t-\t-",
    "t2\tx\t400\tMobiDB-lite\tmobidb-lite\tdisorder\t300\t340\t-\tT\t01-11-2026\t-\t-\t-\t-",
]) + "\n"

FAKE_NEXTFLOW = """#!/bin/sh
# Stand-in for nextflow: '-version' prints a version; 'run' copies $FAKE_IPS_TSV to <outdir>/<outprefix>.tsv
# and records its arguments in $FAKE_IPS_ARGS.
if [ "$1" = "-version" ]; then echo "      version 26.10.4 build 14435"; exit 0; fi
echo "$@" > "$FAKE_IPS_ARGS"
outdir=""; prefix=""; input=""
while [ $# -gt 0 ]; do
  case "$1" in
    --input) input="$2"; shift ;;
    --outdir) outdir="$2"; shift ;;
    --outprefix) prefix="$2"; shift ;;
  esac
  shift
done
echo "Submitted process > INTERPROSCAN:X" > .nextflow.log
echo "Task completed > X" >> .nextflow.log
# Nextflow resolves paths against its working directory, as the real run does
[ -f "$input" ] || { echo "Input file '$input' does not exist"; exit 1; }
cp "$FAKE_IPS_TSV" "$outdir/$prefix.tsv"
"""


def setup(tmp: Path):
    data = tmp / "ipsdata"
    ipr = data / "interpro" / "110.0"
    ipr.mkdir(parents=True)
    (ipr / "databases.json").write_text(json.dumps({"AntiFam": "8.0", "Pfam": "38.2", "PROSITE profiles": "2026_01",
                                                    "InterPro": "110.0"}))
    (ipr / "entries.json").write_text(json.dumps({
        "IPR001881": {"description": "EGF-like calcium-binding domain", "type": "Domain"},
        "IPR001774": {"description": "Delta/Serrate/lag-2 (DSL) protein", "type": "Family"},
        "IPR000742": {"description": "EGF-like domain", "type": "Domain"},
        "IPR000719": {"description": "Protein kinase domain", "type": "Domain"},
    }))
    (ipr / "goterms.ipr.json").write_text(json.dumps({
        "IPR001881": ["GO:0005509"], "IPR001774": ["GO:0007154", "GO:0016020"], "IPR000719": ["GO:0005509"]}))
    for d in ("antifam/8.0", "pfam/38.2", "prosite/2026_01"):
        (data / d).mkdir(parents=True)
        (data / d / "file").write_text("x")
    wf = tmp / "InterProScan" / "conf"
    wf.mkdir(parents=True)
    (wf / "applications.config").write_text(APPS_CONFIG)
    fake = tmp / "nextflow"
    fake.write_text(FAKE_NEXTFLOW)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    tsv = tmp / "fake.tsv"
    tsv.write_text(TSV)
    obo = tmp / "go.obo"
    obo.write_text(OBO)
    return data, tmp / "InterProScan", fake, tsv, obo


def test_helpers(tmp: Path):
    apps = parse_applications_config(tmp / "InterProScan" / "conf" / "applications.config")
    check(apps == {"antifam": "antifam", "pfam": "pfam", "prositeprofiles": "prosite"},
          f"applications with data parsed: {apps}")
    check(domain_label("Protein kinase domain", "Domain") == "protein kinase domain", "capitalized word lowered")
    check(domain_label("EGF-like domain", "Domain") == "EGF-like domain", "acronym kept")
    check(domain_label("WD40 repeat", "Repeat") == "WD40 repeat", "repeat label")
    check(domain_label("Zinc finger, C2H2 type", "Domain") == "zinc finger C2H2 type domain", "comma removed, word added")
    check(domain_label("Domain of unknown function DUF4515", "Domain") == "DUF4515 domain", "DUF form")
    check(domain_label("ABC transporter type 1, transmembrane domain MetI-like", "Domain")
          == "ABC transporter type 1 transmembrane domain", "text after 'domain' removed")
    check(domain_product([(1, "Cadherin-like domain 1", "Domain"), (90, "Cadherin-like domain 2", "Domain")])
          == "uncharacterized Cadherin-like domain-containing protein", "near-identical domains named once")
    check(domain_label("UPF0606 protein KIAA1549, alpha/beta domain 1", "Domain") == ""
          and domain_label("C2orf74 domain", "Domain") == "", "entries named after identifiers give no label")
    check(domain_product([(1, "UPF0606 protein KIAA1549, alpha/beta domain 1", "Domain"),
                          (90, "UPF0606 protein KIAA1549, alpha/beta domain 2", "Domain")]) == "",
          "identifier-named domains alone give no product")
    check(domain_product([(1, "UPF0606 protein KIAA1549, alpha/beta domain 1", "Domain"), (300, "SH2 domain", "Domain")])
          == "uncharacterized SH2 domain-containing protein", "identifier-named domains skipped")
    check(domain_product([(20, "EGF-like domain", "Domain")]) == "uncharacterized EGF-like domain-containing protein",
          "one domain")
    check(domain_product([(150, "Protein kinase domain", "Domain"), (20, "EGF-like domain", "Domain")])
          == "uncharacterized EGF-like and protein kinase domain-containing protein", "two domains, by position")
    check(domain_product([(5, "WD40 repeat", "Repeat"), (50, "SH2 domain", "Domain")])
          == "uncharacterized SH2 domain and WD40 repeat-containing protein", "domain before repeat")
    check(domain_product([(5, "Some family", "Family")]) == "", "families give no name")
    check("--cpus" not in resources_config("singularity") and "cpus = 4" in resources_config("singularity"),
          "no container CPU limit for Singularity; step CPUs kept")
    m = parse_tsv(tmp / "fake.tsv")
    check(len(m["t1"]) == 3 and m["t1"][1]["go"] == ["GO:0007154", "GO:0016020"] and m["t2"][2]["interpro"] == "",
          "TSV parsed")


def test_module(tmp: Path, data, workflow, fake, tsv, obo):
    os.environ["FAKE_IPS_TSV"] = str(tsv)
    os.environ["FAKE_IPS_ARGS"] = str(tmp / "args.txt")
    proteins = tmp / "proteins.faa"
    proteins.write_text(">t1 gene=g1\n" + "M" + "A" * 299 + "*\n>t2 gene=g2\n" + "M" + "K" * 399 + "\n"
                        ">t3 gene=g3\nMSTE\n")
    tx2g = {"t1": "g1", "t2": "g2", "t3": "g3"}
    universe = GeneUniverse(gene_ids=["g1", "g2", "g3"], transcript_to_gene=tx2g,
                            gene_to_transcripts={"g1": ["t1"], "g2": ["t2"], "g3": ["t3"]},
                            coding={"g1": True, "g2": True, "g3": True})
    a1 = QueryAnnotation(query_id="g1", consensus_symbol="DLL1", consensus_name="delta-like protein 1",
                         go_terms={"GO:0007165"}, specific_go_terms={"GO:0007165"},
                         go_evidence={"GO:0007165": "IDA"}, go_support={"GO:0007165": 3})
    a2 = QueryAnnotation(query_id="g2")
    a3 = QueryAnnotation(query_id="g3")
    annotations = {"g1": a1, "g2": a2, "g3": a3}
    # Relative paths, as the pipeline passes them (output directory relative to the current directory)
    os.chdir(tmp)
    work = Path("intermediate")
    reports = Path("reports")
    work.mkdir()
    reports.mkdir()
    proteins = Path("proteins.faa")
    gh = GOHierarchy(str(obo))
    args = Namespace(modules="interproscan", interproscan_nextflow=str(fake), interproscan_workflow=str(workflow),
                     interproscan_revision=None, interproscan_datadir=str(data), interproscan_interpro=None,
                     interproscan_profile="docker", interproscan_apps=None, interproscan_local=False,
                     interproscan_keep_work=False, interpro_go="add", interpro_products="fill", nxf_home=None,
                     force_modules=False)
    mods = select_modules(args)
    check(len(mods) == 1 and isinstance(mods[0], InterProScanModule), "interproscan module selected")
    ctx = ModuleContext(query_fasta=proteins, search_mode="blastp", work_dir=work, reports_dir=reports, prefix="p",
                        threads=4, universe=universe, annotations=annotations, go_hierarchy=gh,
                        forbidden_go=set(), namespaces=["BP", "MF", "CC"])
    ran = run_modules(mods, ctx)
    check(len(ran) == 1, "module ran with the stand-in nextflow")
    m = ran[0]
    check(m.status.available and "InterPro 110.0" in m.status.version and "Nextflow 26.10.4" in m.status.version,
          f"status: {m.status.version}")
    argv = (tmp / "args.txt").read_text().split()
    check("--goterms" in argv and "-resume" in argv and "--max-workers" in argv
          and argv[argv.index("--max-workers") + 1] == "4"
          and "--no-matches-api" not in argv and argv[argv.index("--interpro") + 1] == "110.0"
          and argv[argv.index("-profile") + 1] == "docker", f"command line: {' '.join(argv)}")
    conf = Path(argv[argv.index("-c") + 1])
    check(conf.name == RESOURCES_CONFIG and "withName: 'RUN_MOBIDBLITE' { cpus = 4 }" in conf.read_text()
          and '--cpus ${task.cpus}' in conf.read_text(), "resources configuration passed with -c")

    # GO: g1 keeps signal transduction; cell communication (its ancestor) adds nothing;
    # calcium ion binding and membrane are added
    check(a1.specific_go_terms == {"GO:0007165", "GO:0005509", "GO:0016020"},
          f"g1 specific terms without redundancy: {sorted(a1.specific_go_terms)}")
    check(m.gene_go_added["g1"] == ["GO:0005509", "GO:0016020"], "g1 added terms exclude the ancestor")
    check(a1.go_evidence["GO:0007165"] == "IDA" and a1.go_evidence["GO:0005509"] == "IEA", "evidence codes")
    check(a1.consensus_name == "delta-like protein 1", "existing product kept")
    # g2: no product and no GO before
    check(a2.consensus_name == "uncharacterized EGF-like and protein kinase domain-containing protein"
          and a2.product_note == "interpro_domain", f"g2 domain product: {a2.consensus_name!r}")
    check(a2.specific_go_terms == {"GO:0005509"}, "g2 GO from InterPro")
    check(m.counts == {"proteins_with_matches": 2, "genes_with_matches": 2, "genes_with_interpro": 2,
                       "genes_with_pfam": 2, "genes_go_added": 2, "go_terms_added": 3, "genes_go_only_interpro": 1,
                       "products_filled": 1}, f"counts {m.counts}")
    check(m._entries_text(m.gene_entries["g2"]) == "IPR000742:EGF-like domain;IPR000719:Protein kinase domain",
          "gene entries by position")
    check(m.gene_attributes("g1") == {"Dbxref": "InterPro:IPR001774,InterPro:IPR000742,InterPro:IPR001881,Pfam:PF00008"},
          f"gene Dbxref: {m.gene_attributes('g1')}")
    check(not (work / "p_interproscan" / "work").exists(), "work directory removed")
    check(argv[argv.index("--input") + 1].startswith("/"), "absolute input path passed to Nextflow")
    rep = (reports / "p_interproscan.tsv").read_text().splitlines()
    check(len(rep) == 7 and rep[1].startswith("g1\tt1\tSMART\tSM00051"), "report sorted by position")

    # GO evidence file: source column
    ev = tmp / "go_evidence.tsv"
    output.write_go_evidence(annotations, gh, ev)
    rows = {(r.split("\t")[0], r.split("\t")[1]): r.split("\t") for r in ev.read_text().splitlines()[1:]}
    check(rows[("g1", "GO:0007165")][8] == "hits" and rows[("g1", "GO:0005509")][8] == "interpro"
          and rows[("g1", "GO:0005509")][9] == "IPR001881", "GO evidence source column")

    # GFF3: Dbxref merged with the existing value
    gff_in = tmp / "in.gff3"
    gff_in.write_text("##gff-version 3\nchr1\tsrc\tgene\t1\t900\t.\t+\t.\tID=g1;Dbxref=GeneID:1\n"
                      "chr1\tsrc\tmRNA\t1\t900\t.\t+\t.\tID=t1;Parent=g1\n")
    t1 = Transcript(id="t1", ftype="mRNA", seqid="chr1", strand="+", start=1, end=900, key="t1")
    g1 = Gene(id="g1", ftype="gene", seqid="chr1", strand="+", start=1, end=900, transcripts=[t1], key="g1")
    out = tmp / "out.gff3"
    output.write_annotated_gff3(str(gff_in), out, annotations, [g1], mode="add", hypothetical=False, modules=ran)
    gene_line = [l for l in out.read_text().splitlines() if "\tgene\t" in l][0]
    check("Dbxref=GeneID:1,InterPro:IPR001774" in gene_line and "Pfam:PF00008" in gene_line,
          f"Dbxref merged: {gene_line.split(chr(9))[8]}")

    # Reuse, local mode and fill mode
    (tmp / "args.txt").unlink()
    m2 = InterProScanModule(args)
    m2.detect()
    m2.run(ctx)
    check(not (tmp / "args.txt").exists() and len(m2.matches) == 2, "results reused")
    args_local = Namespace(**{**vars(args), "interproscan_local": True})
    m3 = InterProScanModule(args_local)
    m3.detect()
    m3.run(ctx)
    argv = (tmp / "args.txt").read_text().split()
    check("--no-matches-api" in argv and "--skip-interpro-version-check" in argv, "local mode reruns without lookup")
    m7 = InterProScanModule(Namespace(**{**vars(args), "interproscan_workers": 3}))
    m7.detect()
    ctx.force = True
    m7.run(ctx)
    ctx.force = False
    argv = (tmp / "args.txt").read_text().split()
    check(argv[argv.index("--max-workers") + 1] == "3", "--interproscan-workers sets the parallel tasks")
    b1 = QueryAnnotation(query_id="g1", go_terms={"GO:0007165"}, specific_go_terms={"GO:0007165"},
                         go_evidence={"GO:0007165": "IDA"}, go_support={"GO:0007165": 1}, consensus_name="x")
    ctx_fill = ModuleContext(**{**vars(ctx), "annotations": {"g1": b1, "g2": QueryAnnotation(query_id="g2"),
                                                             "g3": QueryAnnotation(query_id="g3")}})
    m4 = InterProScanModule(Namespace(**{**vars(args), "interpro_go": "fill"}))
    m4.detect()
    m4.run(ctx_fill)
    m4.integrate(ctx_fill)
    check(b1.specific_go_terms == {"GO:0007165"} and m4.counts["genes_go_added"] == 1,
          "fill mode leaves genes with GO from hits alone")

    # Missing data and missing workflow
    (data / "pfam" / "38.2" / "file").unlink()
    m5 = InterProScanModule(args)
    s = m5.detect()
    check(not s.available and "pfam 38.2" in s.reason, f"missing data reported: {s.reason}")
    m6 = InterProScanModule(Namespace(**{**vars(args), "interproscan_workflow": None}))
    check(not m6.detect().available, "workflow required")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        data, workflow, fake, tsv, obo = setup(tmp)
        test_helpers(tmp)
        test_module(tmp, data, workflow, fake, tsv, obo)
    print()
    if failures:
        print(f"{failures} check(s) failed")
        sys.exit(1)
    print("All checks passed")
