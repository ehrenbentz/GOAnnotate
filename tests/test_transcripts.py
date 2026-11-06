#!/usr/bin/env python3
"""
Checks of the transcript input handling on small inputs (no databases
needed): translation in the frame of the best blastx hit, and the grouping
of Trinity isoforms into genes.

Run: python3 tests/test_transcripts.py
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import orfs, sequences  # noqa: E402

failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


# ---------------------------------------------------------------- translation
print("\n# ORF from a blastx hit")
# 5' UTR with an in-frame stop (TAA), then ATG, 6 codons, TGA, 3' UTR
cds = "ATG" + "GCT" * 6 + "TGA"
utr5 = "CCTAACC"          # 7 nt: frame offset 1 for the CDS; TAA in frame at positions 3-5
utr3 = "GGGGG"
tx = utr5 + cds + utr3
# alignment over codons 3-5 of the CDS (nucleotides 7+7=14 .. 7+15=22 on the transcript)
a_start = len(utr5) + 1 + 6
a_end = a_start + 8
strand, frame, s, e, aa, status, note = orfs.orf_from_hit(tx, a_start, a_end)
check(strand == "+" and aa == "M" + "A" * 6 and status == "complete" and note == "",
      f"plus strand: methionine after the upstream stop to the stop codon ({aa}, {status})")
check(s == len(utr5) + 1 and e == len(utr5) + len(cds), f"ORF coordinates include the stop codon ({s}-{e})")
check(tx[s - 1:e] == cds, "ORF coordinates give the CDS")

rc = sequences.reverse_complement(tx)
n = len(tx)
strand, frame, s, e, aa, status, note = orfs.orf_from_hit(rc, n - a_start + 1, n - a_end + 1)
check(strand == "-" and aa == "M" + "A" * 6 and status == "complete",
      f"minus strand: the same protein ({aa}, {status})")
check(sequences.reverse_complement(rc[s - 1:e]) == cds, f"minus strand: ORF coordinates on the given sequence ({s}-{e})")

# no upstream stop and no methionine: 5' partial from the transcript start
tx2 = "GCT" * 5 + "TGA" + "CC"
_, _, s, e, aa, status, note = orfs.orf_from_hit(tx2, 4, 12)
check(aa == "A" * 5 and status == "5_partial" and s == 1 and e == 18, f"5' partial without a methionine ({aa}, {status})")
# no stop downstream: 3' partial
tx3 = "CC" + "ATG" + "GCT" * 5 + "GC"
_, _, s, e, aa, status, note = orfs.orf_from_hit(tx3, 6, 14)
check(aa == "M" + "A" * 5 and status == "3_partial" and e == 2 + 18, f"3' partial without a stop ({aa}, {status}, end {e})")
# upstream stop but no methionine before the alignment: starts after the stop, noted
tx4 = "TAA" + "GCT" * 4 + "TGA"
_, _, s, e, aa, status, note = orfs.orf_from_hit(tx4, 7, 12)
check(aa == "A" * 4 and status == "5_partial" and note == "no_start_codon" and s == 4,
      f"after an upstream stop without a methionine ({aa}, {note})")
# a stop inside the aligned region ends the protein there
tx5 = "ATG" + "GCT" * 3 + "TAG" + "GCT" * 3 + "TGA"
_, _, s, e, aa, status, note = orfs.orf_from_hit(tx5, 4, 24)
check(aa == "M" + "A" * 3 and note == "internal_stop" and e == 15, f"internal stop ends the protein ({aa}, {note})")

# ---------------------------------------------------------------- translate_from_hits
print("\n# translate_from_hits")
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    fasta = tmp / "tx.fna"
    with open(fasta, "w") as fh:
        fh.write(f">DN1_c0_g1_i1 len={len(tx)}\n{tx}\n>DN1_c0_g1_i2\n{rc}\n>DN2_c0_g1_i1\n{tx2}\n>nohit\n{tx3}\n")
    hits = tmp / "blastx.tsv"
    def row(q, qs, qe, qlen, bits, sid="sp|P1|X_HUMAN", ev=1e-20):
        return f"{q}\t{sid}\t60.0\t3\t0\t0\t{qs}\t{qe}\t1\t3\t{ev}\t{bits}\t{qlen}\t100\t{sid} X\n"
    with open(hits, "w") as fh:
        fh.write(row("DN1_c0_g1_i1", a_start, a_end, len(tx), 50.0))
        fh.write(row("DN1_c0_g1_i1", 1, 6, len(tx), 20.0, sid="sp|P2|Y_HUMAN"))          # other frame, separate region
        fh.write(row("DN1_c0_g1_i2", n - a_start + 1, n - a_end + 1, len(rc), 50.0))
        fh.write(row("DN2_c0_g1_i1", 4, 12, len(tx2), 40.0))
        fh.write(row("DN2_c0_g1_i1", 4, 12, len(tx2), 10.0, ev=1e-3))                   # over the threshold
    result = orfs.translate_from_hits(str(fasta), str(hits), tmp / "AA.faa", tmp / "orf.tsv")
    records = dict((h.split()[0], (h, s)) for h, s in sequences.read_fasta(str(tmp / "AA.faa")))
    check(result.n_transcripts == 4 and result.n_with_hits == 3 and result.n_translated == 3,
          f"three of four transcripts translated ({result.n_translated})")
    check(set(records) == {"DN1_c0_g1_i1", "DN1_c0_g1_i2", "DN2_c0_g1_i1"}, "transcripts without hits get no protein")
    check(records["DN1_c0_g1_i1"][1] == "M" + "A" * 6 and records["DN1_c0_g1_i2"][1] == "M" + "A" * 6,
          "proteins from the best hit's frame on either strand")
    check(result.completeness == {"complete": 2, "5_partial": 1}, f"completeness counts {dict(result.completeness)}")
    check(result.notes.get("multiple_frames") == 1 and result.notes.get("separate_regions") == 1,
          f"hits in another frame and region flagged {dict(result.notes)}")
    _, fields = sequences.parse_header(records["DN1_c0_g1_i1"][0])
    check(fields.get("orf") == f"{len(utr5) + 1}-{len(utr5) + len(cds)}" and fields.get("strand") == "+"
          and fields.get("orf_status") == "complete", f"ORF fields in the protein header {fields}")
    header_fields = orfs.orf_header_fields(tmp / "orf.tsv")
    check(header_fields["DN1_c0_g1_i2"]["strand"] == "-" and header_fields["DN2_c0_g1_i1"]["orf_status"] == "5_partial",
          "ORF fields read back from the report")
    check(len(open(tmp / "orf.tsv").read().splitlines()) == 4, "report has one row per protein")

# ---------------------------------------------------------------- upstream stop, completeness, blastx fallback
print("\n# upstream stop, completeness, blastx fallback")
tx6 = "TAA" + "CC" + "C" + "ATG" + "GCT" * 5 + "TGA"      # in-frame TAA, then ATG
_, _, s, e, aa, status, note = orfs.orf_from_hit(tx6, 10, 18)
check(aa == "M" + "A" * 5 and status == "complete" and note == "upstream_stop",
      f"a methionine after an in-frame stop is recorded ({note})")
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    report = tmp / "orf.tsv"
    hdr = "transcript_id\tstrand\tframe\torf_start\torf_end\tprotein_length\tcompleteness\tdefining_hit\tbitscore\thits\tnote\n"
    rows = [("a_i1", "complete", "upstream_stop"),     # methionine after a stop, stop codon
            ("b_i1", "complete", ""),                  # methionine at the transcript start, reference start reached
            ("c_i1", "complete", ""),                  # methionine, no upstream stop, reference start not reached
            ("d_i1", "5_partial", "no_start_codon"),
            ("e_i1", "3_partial", "upstream_stop"),    # no stop, but the reference end is reached
            ("f_i1", "partial", "")]
    with open(report, "w") as fh:
        fh.write(hdr)
        for t, c, n in rows:
            fh.write(f"{t}\t+\t1\t1\t300\t100\t{c}\tsp|X|Y\t50\t3\t{n}\n")
    from goannotate.models import QueryAnnotation
    annotations = {t: QueryAnnotation(query_id=t, consensus_symbol="ABC") for t, _, _ in rows}
    ref = {"b_i1": {"ABC": (90.0, 3, 100, 300)}, "c_i1": {"ABC": (90.0, 80, 300, 300), "DEF": (95.0, 1, 50, 60)},
           "e_i1": {"ABC": (90.0, 40, 295, 300)}}
    comp = orfs.completeness(report, ref, annotations, {})
    check(comp["a_i1"][0] == "complete" and "5:met+upstream_stop" in comp["a_i1"][1], f"stop-bounded methionine: {comp['a_i1']}")
    check(comp["b_i1"][0] == "complete" and "5:met+reference_start" in comp["b_i1"][1], f"reference start reached: {comp['b_i1']}")
    check(comp["c_i1"][0] == "5_partial" and "5:met_without_support" in comp["c_i1"][1],
          f"the symbol's reference protein decides, not the best other hit: {comp['c_i1']}")
    check(comp["d_i1"][0] == "5_partial" and comp["f_i1"][0] == "partial", "missing methionine and stop")
    check(comp["e_i1"][0] == "3_partial" and "3:reference_end" in comp["e_i1"][1], f"no stop codon: {comp['e_i1']}")
    fields = orfs.orf_header_fields(report, comp)
    check(fields["c_i1"]["completeness"] == "5_partial" and fields["c_i1"]["orf_status"] == "complete",
          "header fields carry the reading frame status and the completeness")
    blastp = tmp / "blastp.tsv"
    blastx = tmp / "blastx.tsv"
    blastp.write_text("a_i1\tsp|P|Q\t50\t30\t0\t0\t1\t30\t1\t30\t1e-10\t60\t100\t100\tx\n")
    blastx.write_text("a_i1\tsp|P|Q\t50\t30\t0\t0\t1\t90\t1\t30\t1e-10\t60\t300\t100\tx\n"
                      "d_i1\tsp|P|Q\t50\t30\t0\t0\t300\t211\t1\t30\t1e-10\t60\t300\t100\tx\n")
    merged, n = orfs.blastx_fallback(blastp, str(blastx), report, tmp / "merged.tsv")
    lines = merged.read_text().splitlines()
    d = [ln.split("\t") for ln in lines if ln.startswith("d_i1")]
    check(n == 1 and len(lines) == 2 and not any(ln.startswith("a_i1") and "\t90\t" in ln for ln in lines),
          "only transcripts without blastp hits take their blastx hits")
    check(d and d[0][6] == "1" and d[0][7] == "30" and d[0][12] == "100",
          f"blastx coordinates converted to residues ({d[0][6:8] + [d[0][12]] if d else ''})")

# ---------------------------------------------------------------- Trinity grouping
print("\n# Trinity isoform grouping")
check(sequences.trinity_gene("TRINITY_DN1000_c0_g1_i1") == "TRINITY_DN1000_c0_g1"
      and sequences.trinity_gene("DN26003_c0_g1_i2") == "DN26003_c0_g1"
      and sequences.trinity_gene("TR1000|c0_g1_i1") == "TR1000|c0_g1"
      and sequences.trinity_gene("comp1000_c0_seq2") == "comp1000_c0"
      and sequences.trinity_gene("XM_012345.1") == "" and sequences.trinity_gene("gene-LOC1_i1") == "",
      "Trinity identifiers of each generation give their gene; others give nothing")
with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "tx.fna"
    with open(fasta, "w") as fh:
        fh.write(">DN1_c0_g1_i1\nACGT\n>DN1_c0_g1_i2\nACGT\n>DN2_c0_g1_i1\nACGT\n>other\nACGT\n>withgene gene=G9\nACGT\n")
    u = sequences.universe_from_fasta(str(fasta))
    check(u.gene_ids == ["DN1_c0_g1", "DN2_c0_g1", "other", "G9"] and u.gene_to_transcripts["DN1_c0_g1"] == ["DN1_c0_g1_i1", "DN1_c0_g1_i2"],
          f"isoforms grouped by Trinity identifier, gene= field kept ({u.gene_ids})")
    u = sequences.universe_from_fasta(str(fasta), isoform_grouping="off")
    check(u.gene_ids == ["DN1_c0_g1_i1", "DN1_c0_g1_i2", "DN2_c0_g1_i1", "other", "G9"], "grouping off: one gene per transcript")

print()
if failures:
    print(f"{len(failures)} check(s) failed:")
    for f in failures:
        print("  - " + f)
    sys.exit(1)
print("All checks passed")
