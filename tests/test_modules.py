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
from goannotate.modules.signalp import (SignalPModule, alternative_start, alternative_starts, hydrophobic_core,  # noqa: E402
                                        parse_prediction_results, plan_processes)
from goannotate.models import BlastHit  # noqa: E402
from goannotate import output  # noqa: E402

failures = 0


def check(cond, label):
    global failures
    print(("ok  " if cond else "FAIL") + " " + label)
    if not cond:
        failures += 1


FAKE_SIGNALP = """#!/bin/sh
# Stand-in for signalp6, deciding from the sequence (one line per record):
# LLLL within the first 30 residues: SP, probability 0.987655, cleavage after 22;
# WLW within the first 30 residues: SP, probability 0.6, cleavage after 18;
# otherwise OTHER. Writes a tqdm-like progress line and records every input
# sequence in $SIGNALP_FAKE_LOG when that is set.
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
n=$(grep -c '^>' "$fasta")
awk '/^>/{id=substr($0,2); next} {
  head=substr($0,1,30)
  if (index(head,"LLLL")) printf "%s\\tSP\\t0.012345\\t0.987655\\tCS pos: 22-23. Pr: 0.9512\\n", id
  else if (index(head,"WLW")) printf "%s\\tSP\\t0.400000\\t0.600000\\tCS pos: 18-19. Pr: 0.5000\\n", id
  else printf "%s\\tOTHER\\t0.999900\\t0.000100\\t\\n", id
}' "$fasta" >> "$out"
if [ -n "${SIGNALP_FAKE_LOG:-}" ]; then grep -v '^>' "$fasta" >> "$SIGNALP_FAKE_LOG"; fi
printf 'Predicting: 100%%| %s/%s [00:01<00:00]\\n' "$n" "$n" >&2
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

    # g2's isoforms share their first 70 residues (one SignalP input); g4's
    # annotated start has no signal peptide but an internal Met at residue 34
    # starts one; g5's signal peptide has probability 0.6.
    g2_head = "MRLLLLALLAVSSAAEKKPEL" + "A" * 49
    g4_seq = "MSTWWPNTEKNHVHQPYENRSTDEKQPNGSRDE" + "MFQLLLLAVLLASAQAQSSDEEK"
    # g6: the first gated Met (residue 11) has no signal peptide, the second (residue 42) has one
    g6_seq = "MSTDEKPNGS" + "MAVIAVFAGTDEKPNGSRDEQKPNTSRDEQK" + "MFLLLLAVLASAQAQSSDEEK"
    proteins = tmp / "proteins.faa"
    proteins.write_text(">g1_sec_t1 gene=g1\nMKTLLLLAVLLAVAGASAQAQSSDEEKTT*\n"
                        ">g1_t2 gene=g1\nMAAAAKKKDDDEEE\n"
                        f">g2_sec_t1 gene=g2\n{g2_head}KKK\n"
                        f">g2_sec_t2 gene=g2\n{g2_head}DDD\n"
                        ">g3_t1 gene=g3\nMSSSSSSSSTTTT\n"
                        f">g4_t1 gene=g4\n{g4_seq}\n"
                        ">g5_t1 gene=g5\nMAWLWAAQSSDEEK\n"
                        f">g6_t1 gene=g6\n{g6_seq}\n")
    tx2g = {"g1_sec_t1": "g1", "g1_t2": "g1", "g2_sec_t1": "g2", "g2_sec_t2": "g2", "g3_t1": "g3",
            "g4_t1": "g4", "g5_t1": "g5", "g6_t1": "g6"}
    g2t = {}
    for t, g in tx2g.items():
        g2t.setdefault(g, []).append(t)
    universe = GeneUniverse(gene_ids=list(g2t), transcript_to_gene=tx2g, gene_to_transcripts=g2t,
                            coding={g: True for g in g2t})
    annotations = {g: QueryAnnotation(query_id=g, consensus_symbol=g.upper(), consensus_name=f"protein {g}")
                   for g in universe.gene_ids}
    # g4's best hit aligns from the homolog's residue 1 to the query's residue 34 (the internal Met)
    annotations["g4"].hits = [BlastHit("g4_t1", "sp|P99999|TEST_HUMAN", 80.0, 50, 0, 0, 34, 83, 1, 50, 1e-20, 200.0),
                              BlastHit("g4_t1", "tr|A0A000|X_MOUSE", 60.0, 40, 0, 0, 40, 79, 10, 49, 1e-10, 90.0)]
    work = tmp / "intermediate"
    reports = tmp / "reports"
    work.mkdir()
    reports.mkdir()
    fake_log = tmp / "fake_inputs.txt"
    os.environ["SIGNALP_FAKE_LOG"] = str(fake_log)

    args = Namespace(modules="signalp", signalp_bin=str(fake), signalp_mode="fast", signalp_features=True,
                     force_modules=False, signalp_alt_starts=True, signalp_min_probability=0.0, signalp_processes=0)
    selected = select_modules(args)
    check(len(selected) == 1 and isinstance(selected[0], SignalPModule), "signalp module selected")
    ctx = ModuleContext(query_fasta=proteins, search_mode="blastp", work_dir=work, reports_dir=reports,
                        prefix="test", threads=2, universe=universe, annotations=annotations)
    ran = run_modules(selected, ctx)
    check(len(ran) == 1, "module ran with the stand-in executable")
    m = ran[0]
    check(m.status.available and m.status.version == "6.0h", "version parsed from --version")
    check(m.processes == 2 and len(list((work / "test_signalp" / "annotated_starts").glob("chunk_*"))) == 2,
          "two processes for --threads 2, one chunk each")
    inputs = fake_log.read_text().split()
    check(len(inputs) == 10 and all(len(s) <= 70 for s in inputs),
          f"7 distinct annotated starts + 3 alternative starts sent, each at most 70 residues ({len(inputs)})")
    check(m.counts == {"proteins_scored": 8, "distinct_starts": 7, "transcripts_with_sp": 4, "below_cutoff": 0,
                       "genes_with_sp": 3, "genes_tested": 6, "alt_candidates": 2, "alt_distinct": 3,
                       "alt_sp": 2, "alt_supported": 1}, f"counts {m.counts}")
    check(m._tx_value("g2_sec_t2") == "yes" and m._tx_cleavage("g2_sec_t2") == "22-23",
          "result copied to the isoform with the same N-terminus")
    check(m._gene_value("g1") == "yes (1/2 isoforms)" and m._gene_value("g2") == "yes"
          and m._gene_value("g3") == "no" and m._gene_value("g9") == "", "gene column values")
    check(m._tx_value("g1_sec_t1") == "yes" and m._tx_probability("g1_sec_t1") == "0.988"
          and m._tx_cleavage("g1_sec_t1") == "22-23" and m._tx_probability("g1_t2") == "", "transcript columns")
    check(m._tx_value("g4_t1") == "no" and
          m._tx_alt("g4_t1") == "M34: SP 0.99, cleavage 55-56; homolog P99999 aligns from residue 34",
          f"alternative start in protein coordinates with homology note: {m._tx_alt('g4_t1')!r}")
    secretome = (reports / "test_secretome.tsv").read_text().splitlines()
    check(len(secretome) == 5 and secretome[1].startswith("g1\tg1_sec_t1\tG1\tprotein g1\t0.988\t22-23\t0.951"),
          "secretome report")
    alt_rows = (reports / "test_signalp_alternative_starts.tsv").read_text().splitlines()
    alt_g4 = [r.split("\t") for r in alt_rows[1:] if r.startswith("g4\t")]
    check(len(alt_rows) == 3 and alt_g4 and alt_g4[0][5:8] == ["M34", "0.988", "55-56"], "alternative-start report")
    merged = parse_prediction_results(work / "test_signalp" / "annotated_starts" / "prediction_results.txt")
    check(len(merged) == 7, "chunk outputs merged")

    # Probability cutoff: g5 (0.6) is no longer called; results reused
    args_cut = Namespace(**{**vars(args), "signalp_min_probability": 0.7})
    m3 = SignalPModule(args_cut)
    m3.detect()
    m3.run(ctx)
    m3.integrate(ctx)
    check(m3.runtime == 0.0, "results reused for the cutoff run")
    check(m3.counts["transcripts_with_sp"] == 3 and m3.counts["below_cutoff"] == 1 and m3._gene_value("g5") == "no",
          "cutoff removes the 0.6 call")

    # Reuse of existing results
    m2 = SignalPModule(args)
    m2.detect()
    m2.run(ctx)
    check(m2.runtime == 0.0 and len(m2.results) == 8, "existing results reused without rerunning")
    check(m._tx_alt("g6_t1") == "M42: SP 0.99, cleavage 63-64",
          f"a later Met is tested when the first gated Met has no signal peptide: {m._tx_alt('g6_t1')!r}")

    # Process planning
    check(plan_processes(8, "fast", requested=3) == (3, "--signalp-processes"), "requested process count used")
    n, _ = plan_processes(10_000, "fast")
    check(1 <= n <= (os.cpu_count() or 1), f"automatic process count within the CPUs ({n})")

    # Alternative-start gate
    check(alternative_start(g4_seq) == 33 and alternative_start("MSTWWPNTEK" * 8) == -1, "internal Met gate")
    check(alternative_starts(g6_seq) == [10, 41], f"every gated Met listed: {alternative_starts(g6_seq)}")
    check(hydrophobic_core("MKLLLLAVLLAVAG") and not hydrophobic_core("MKDEKDEKDEKDEKDE"), "hydrophobic core test")

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
    check(th[0].endswith("\tgo_terms\tsignal_peptide\tsp_probability\tsp_cleavage\tsp_alternative_start"),
          "transcript table columns appended")
    g1_row = [l for l in gh if l.startswith("g1\t")][0]
    check(g1_row.endswith("\tyes (1/2 isoforms)"), "gene table value")
    sec_row = [l for l in th if "\tg1_sec_t1\t" in l][0]
    check(sec_row.endswith("\tyes\t0.988\t22-23\t"), "transcript table values")
    g4_row = [l for l in th if "\tg4_t1\t" in l][0]
    check(g4_row.endswith("\tno\t\t\tM34: SP 0.99, cleavage 55-56; homolog P99999 aligns from residue 34"),
          "alternative-start column")


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
