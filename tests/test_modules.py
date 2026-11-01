"""
Checks of the optional-module framework and the SignalP module on small
inputs, without SignalP installed: a stand-in signalp6 script writes a
prediction_results.txt in the SignalP 6 layout.

Run: python3 tests/test_modules.py
"""

import os
import stat
import sys
import tempfile
from argparse import Namespace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.gff3 import Transcript, Gene  # noqa: E402
from goannotate.models import GeneUniverse, QueryAnnotation  # noqa: E402
from goannotate.modules import parse_selection, run_modules, select_modules, ModuleContext  # noqa: E402
from goannotate.modules.base import protein_to_genomic, write_clean_proteins  # noqa: E402
from goannotate.modules.signalp import SignalPModule, parse_prediction_results  # noqa: E402
from goannotate import output  # noqa: E402

failures = 0


def check(cond, label):
    global failures
    print(("ok  " if cond else "FAIL") + " " + label)
    if not cond:
        failures += 1


FAKE_SIGNALP = """#!/bin/sh
# Stand-in for signalp6: every sequence whose ID contains 'sec' gets a
# signal peptide with cleavage after residue 22; the others are OTHER.
fasta=""; outdir=""
while [ $# -gt 0 ]; do
  case "$1" in
    --version) echo "SignalP 6.0 Signal peptide prediction tool 6.0h"; exit 0 ;;
    --fastafile) fasta="$2"; shift ;;
    --output_dir) outdir="$2"; shift ;;
  esac
  shift
done
mkdir -p "$outdir"
out="$outdir/prediction_results.txt"
printf '# SignalP-6.0\\tOrganism: Eukarya\\tTimestamp: test\\n' > "$out"
printf '# ID\\tPrediction\\tOTHER\\tSP(Sec/SPI)\\tCS Position\\n' >> "$out"
grep '^>' "$fasta" | sed 's/^>//' | while read -r id; do
  case "$id" in
    *sec*) printf '%s\\tSP\\t0.012345\\t0.987655\\tCS pos: 22-23. Pr: 0.9512\\n' "$id" >> "$out" ;;
    *)     printf '%s\\tOTHER\\t0.999900\\t0.000100\\t\\n' "$id" >> "$out" ;;
  esac
done
"""


def test_parse_selection():
    check(parse_selection(None) == [], "no --modules: nothing selected")
    check(parse_selection("signalp") == ["signalp"], "one module")
    check(parse_selection("signalp, signalp") == ["signalp"], "duplicates collapsed")
    check(parse_selection("auto") == ["auto"], "auto")
    try:
        parse_selection("interpro")
        check(False, "unknown module rejected")
    except ValueError:
        check(True, "unknown module rejected")


def test_parse_results(tmp):
    p = tmp / "prediction_results.txt"
    p.write_text("# SignalP-6.0\tOrganism: Eukarya\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
                 "t1\tSP\t0.0123\t0.9877\tCS pos: 22-23. Pr: 0.9512\n"
                 "t2\tOTHER\t0.9999\t0.0001\t\n"
                 "t3 some description\tSP\t0.2\t0.8\tCS pos: 19-20. Pr: 0.61\n")
    r = parse_prediction_results(p)
    check(set(r) == {"t1", "t2", "t3"}, "three results parsed, ID up to first space")
    check(r["t1"].has_sp and abs(r["t1"].probability - 0.9877) < 1e-6, "SP probability from its column")
    check(r["t1"].cs_start == 22 and r["t1"].cs_end == 23 and abs(r["t1"].cs_probability - 0.9512) < 1e-6,
          "cleavage site parsed")
    check(not r["t2"].has_sp and abs(r["t2"].probability - 0.9999) < 1e-6, "OTHER probability")
    check(r["t3"].cs_start == 19, "second SP result")


def test_clean_proteins(tmp):
    p = tmp / "clean.faa"
    n = write_clean_proteins([("a", "MKT*LLU*"), ("b", "m" * 300), ("c", "*")], p, max_length=200)
    recs = dict((h, s) for h, s in _read(p))
    check(n == 2 and set(recs) == {"a", "b"}, "empty record skipped")
    check(recs["a"] == "MKTXLLX", "internal stop and U replaced by X, trailing stop removed")
    check(len(recs["b"]) == 200, "truncated to 200 residues")


def _read(path):
    from goannotate.sequences import read_fasta
    return read_fasta(str(path))


def test_protein_to_genomic():
    # Plus strand, two CDS segments 101-130 (30 nt = 10 aa) and 201-260, phase 0
    tx = Transcript(id="t", ftype="mRNA", seqid="chr1", strand="+", start=101, end=260,
                    cds=[(101, 130, 0), (201, 260, 0)])
    check(protein_to_genomic(tx, 1, 5) == [(101, 115)], "plus strand within first segment")
    check(protein_to_genomic(tx, 9, 12) == [(125, 130), (201, 206)], "plus strand across the intron")
    # Minus strand: translation starts at the high coordinate
    tx2 = Transcript(id="t2", ftype="mRNA", seqid="chr1", strand="-", start=101, end=260,
                     cds=[(101, 160, 0), (231, 260, 0)])
    check(protein_to_genomic(tx2, 1, 5) == [(246, 260)], "minus strand first residues at the high end")
    check(protein_to_genomic(tx2, 9, 12) == [(155, 160), (231, 236)], "minus strand across the intron")
    # Phase 1: the first nucleotide is skipped
    tx3 = Transcript(id="t3", ftype="mRNA", seqid="chr1", strand="+", start=101, end=160,
                     cds=[(101, 160, 1)])
    check(protein_to_genomic(tx3, 1, 2) == [(102, 107)], "phase 1 skips one nucleotide")
    nc = Transcript(id="nc", ftype="transcript", seqid="chr1", strand="+", start=1, end=10)
    check(protein_to_genomic(nc, 1, 5) == [], "non-coding transcript gives nothing")


def test_module_end_to_end(tmp):
    bin_dir = tmp / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "signalp6"
    fake.write_text(FAKE_SIGNALP)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)

    proteins = tmp / "proteins.faa"
    proteins.write_text(">g1_sec_t1 gene=g1\nMKTLLLLAVLLAVAGASAQAQSSDEEKTT*\n"
                        ">g1_t2 gene=g1\nMAAAAKKKDDDEEE\n"
                        ">g2_sec_t1 gene=g2\nMRLLPLLALLAVSSAAEKKPEL\n"
                        ">g3_t1 gene=g3\nMSSSSSSSSTTTT\n")
    universe = GeneUniverse(
        gene_ids=["g1", "g2", "g3"],
        transcript_to_gene={"g1_sec_t1": "g1", "g1_t2": "g1", "g2_sec_t1": "g2", "g3_t1": "g3"},
        gene_to_transcripts={"g1": ["g1_sec_t1", "g1_t2"], "g2": ["g2_sec_t1"], "g3": ["g3_t1"]},
        coding={"g1": True, "g2": True, "g3": True})
    annotations = {g: QueryAnnotation(query_id=g, consensus_symbol=g.upper(), consensus_name=f"protein {g}")
                   for g in universe.gene_ids}
    work = tmp / "intermediate"
    reports = tmp / "reports"
    work.mkdir()
    reports.mkdir()

    args = Namespace(modules="signalp", signalp_bin=str(fake), signalp_mode="fast", signalp_features=True,
                     force_modules=False)
    selected = select_modules(args)
    check(len(selected) == 1 and isinstance(selected[0], SignalPModule), "signalp module selected")
    ctx = ModuleContext(query_fasta=proteins, search_mode="blastp", work_dir=work, reports_dir=reports,
                        prefix="test", threads=2, universe=universe, annotations=annotations)
    ran = run_modules(selected, ctx)
    check(len(ran) == 1, "module ran with the stand-in executable")
    m = ran[0]
    check(m.status.available and m.status.version == "6.0h", "version parsed from --version")
    check(m.counts == {"proteins_scored": 4, "transcripts_with_sp": 2, "genes_with_sp": 2, "genes_tested": 3},
          f"counts {m.counts}")
    check(m._gene_value("g1") == "yes (1/2 isoforms)" and m._gene_value("g2") == "yes"
          and m._gene_value("g3") == "no" and m._gene_value("g4") == "", "gene column values")
    check(m._tx_value("g1_sec_t1") == "yes" and m._tx_probability("g1_sec_t1") == "0.988"
          and m._tx_cleavage("g1_sec_t1") == "22-23" and m._tx_probability("g1_t2") == "", "transcript columns")
    secretome = (reports / "test_secretome.tsv").read_text().splitlines()
    check(len(secretome) == 3 and secretome[1].startswith("g1\tg1_sec_t1\tG1\tprotein g1\t0.988\t22-23\t0.951"),
          "secretome report")

    # Reuse of existing results
    m2 = SignalPModule(args)
    m2.detect()
    m2.run(ctx)
    check(m2.runtime == 0.0 and len(m2.results) == 4, "existing results reused without rerunning")

    # Blastx input skips the module
    ctx_nt = ModuleContext(query_fasta=proteins, search_mode="blastx", work_dir=work, reports_dir=reports,
                           prefix="test", threads=2, universe=universe, annotations=annotations)
    skipped = run_modules(select_modules(args), ctx_nt)
    check(skipped == [], "module skipped for nucleotide input")

    # Missing executable
    args_missing = Namespace(modules="signalp", signalp_bin=str(tmp / "nowhere"), signalp_mode="fast",
                             signalp_features=False, force_modules=False)
    mods = select_modules(args_missing)
    check(run_modules(mods, ctx) == [] and "nowhere" in mods[0].skipped_reason, "missing tool reported")

    # GFF3 output: attributes and sig_peptide features
    gff_in = tmp / "in.gff3"
    gff_in.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t101\t260\t.\t+\t.\tID=g1\n"
        "chr1\tsrc\tmRNA\t101\t260\t.\t+\t.\tID=g1_sec_t1;Parent=g1\n"
        "chr1\tsrc\tCDS\t101\t130\t.\t+\t0\tID=cds1;Parent=g1_sec_t1\n"
        "chr1\tsrc\tCDS\t201\t260\t.\t+\t0\tID=cds1;Parent=g1_sec_t1\n"
        "chr1\tsrc\tmRNA\t101\t260\t.\t+\t.\tID=g1_t2;Parent=g1\n"
        "chr1\tsrc\tCDS\t101\t160\t.\t+\t0\tID=cds2;Parent=g1_t2\n")
    t1 = Transcript(id="g1_sec_t1", ftype="mRNA", seqid="chr1", strand="+", start=101, end=260,
                    cds=[(101, 130, 0), (201, 260, 0)], key="g1_sec_t1", exons=[(101, 130), (201, 260)])
    t2 = Transcript(id="g1_t2", ftype="mRNA", seqid="chr1", strand="+", start=101, end=260,
                    cds=[(101, 160, 0)], key="g1_t2", exons=[(101, 160)])
    g1 = Gene(id="g1", ftype="gene", seqid="chr1", strand="+", start=101, end=260, transcripts=[t1, t2], key="g1")
    gff_out = tmp / "out.gff3"
    output.write_annotated_gff3(str(gff_in), gff_out, annotations, [g1], mode="add", hypothetical=False, modules=ran)
    lines = gff_out.read_text().splitlines()
    gene_line = [l for l in lines if "\tgene\t" in l][0]
    mrna1 = [l for l in lines if "ID=g1_sec_t1;" in l][0]
    mrna2 = [l for l in lines if "ID=g1_t2;" in l][0]
    sig = [l for l in lines if "\tsig_peptide\t" in l]
    check("goannotate_signal_peptide=yes" in gene_line, "gene attribute written")
    check("goannotate_signal_peptide=yes" in mrna1 and "goannotate_sp_cleavage=22-23" in mrna1,
          "transcript attributes written")
    check("goannotate_signal_peptide" not in mrna2, "no attribute on the isoform without SP")
    # residues 1-22 = nt 1-66: 101-130 (30 nt) then 201-236 (36 nt)
    check(len(sig) == 2 and sig[0].split("\t")[3:5] == ["101", "130"] and sig[1].split("\t")[3:5] == ["201", "236"],
          f"sig_peptide features across the intron: {[l.split(chr(9))[3:5] for l in sig]}")
    check(all("Parent=g1_sec_t1" in l and "ID=sig_peptide-g1_sec_t1" in l for l in sig), "sig_peptide parent and ID")
    idx_m = lines.index(mrna1)
    check(lines[idx_m + 1] in sig and lines[idx_m + 2] in sig, "sig_peptide lines follow the mRNA line")

    # Tables
    gene_file, tx_file = output.write_annotation_tables(tmp, "test", annotations, universe, genes=None,
                                                        fasta_headers=str(proteins), mode="add", modules=ran)
    gh = gene_file.read_text().splitlines()
    th = tx_file.read_text().splitlines()
    check(gh[0].endswith("\tgo_terms\tsignal_peptide"), "gene table column appended")
    check(th[0].endswith("\tgo_terms\tsignal_peptide\tsp_probability\tsp_cleavage"), "transcript table columns appended")
    g1_row = [l for l in gh if l.startswith("g1\t")][0]
    check(g1_row.endswith("\tyes (1/2 isoforms)"), "gene table value")
    sec_row = [l for l in th if "\tg1_sec_t1\t" in l][0]
    check(sec_row.endswith("\tyes\t0.988\t22-23"), "transcript table values")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_parse_selection()
        test_parse_results(tmp)
        test_clean_proteins(tmp)
        test_protein_to_genomic()
        test_module_end_to_end(tmp)
    print()
    if failures:
        print(f"{failures} check(s) failed")
        sys.exit(1)
    print("All checks passed")
