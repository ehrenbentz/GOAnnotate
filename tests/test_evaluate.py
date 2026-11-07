#!/usr/bin/env python3
"""
Checks of the evaluation report on a small set of annotations: each
source's verdict, the split of the relatives' agreement by whether they
decided the symbol, and the conflict column.

Run: python3 tests/test_evaluate.py
"""
import csv
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import evaluate  # noqa: E402
from goannotate.models import GeneUniverse, QueryAnnotation  # noqa: E402

failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


class FakeOrthology:
    def __init__(self):
        self.decision = {"g1": "confirmed", "g2": "taken", "g3": "none", "g4": "kept"}
        self.ev = {"g1": {"relatives": "ABC1", "close": True, "same_species": "ABC1"},
                   "g2": {"relatives": "XYZ2", "close": True, "same_species": ""},
                   "g3": {"relatives": "", "close": False, "same_species": ""},
                   "g4": {"relatives": "ABC2", "close": True, "same_species": "ABC1"}}

    def ortholog_evidence_for(self, gene):
        return self.ev.get(gene, {})

    def _equivalent(self, a, b):
        return "same" if a.upper() == b.upper() else ""


class FakeInterPro:
    gene_entries = {"g1": [(1, "IPR1"), (50, "IPR2")], "g2": [(1, "IPR3")], "g4": [(1, "IPR1")]}
    entries = {"IPR1": ("Kinesin motor domain", "Domain"), "IPR2": ("P-loop containing nucleoside triphosphate hydrolase", "Homologous_superfamily"),
               "IPR3": ("Zinc finger C2H2-type", "Domain")}


def annot(gene, symbol, name, hit_symbols, sym_conc, hit_names, name_conc, status="rbh", source="uniprot", ref_sym=""):
    a = QueryAnnotation(query_id=gene, consensus_symbol=symbol, consensus_name=name, hits=[1], hit_symbols=hit_symbols,
                        symbol_concordance=sym_conc, hit_names=hit_names, name_concordance=name_conc,
                        symbol_status=status, symbol_source=source)
    a.ortholog_reference_symbol = ref_sym
    return a


annotations = {
    # everything agrees
    "g1": annot("g1", "ABC1", "kinesin-like protein ABC1", {"ABC1": 9, "ABC2": 1}, 0.9, {"kinesin-like protein abc1": 8}, 0.8),
    # taken from relatives: hits disagree, the reference does not hit XYZ2, product names another gene, domains disagree
    "g2": annot("g2", "XYZ2", "kinesin-like protein ABC1", {"ABC1": 10}, 0.0, {"kinesin-like protein abc1": 10}, 1.0,
                status="ortholog", source="relative"),
    # no symbol: hits name only
    "g3": annot("g3", "", "uncharacterized protein", {}, 0.0, {"uncharacterized protein": 3}, 1.0),
    # near-best reference, relatives disagree (independent), same species agrees
    "g4": annot("g4", "ABC1", "kinesin-like protein ABC1", {"ABC1": 5, "ABC2": 5}, 0.5, {"kinesin-like protein abc1": 5}, 0.5),
    "g5": QueryAnnotation(query_id="g5"),      # nothing: not evaluated
    # the prior picked a group the hits scored lower: decided by the relatives; product in the hits' style
    "g6": annot("g6", "ABC2", "Protein OS-9", {"ABC1": 8, "ABC2": 2}, 0.2, {"protein os-9": 8}, 1.0,
                source="relative_override"),
}
annotations["g6"].hit_curated_names = {"protein os-9": ["ABC2"]}
universe = GeneUniverse(gene_ids=list(annotations), transcript_to_gene={f"{g}_t": g for g in annotations},
                        gene_to_transcripts={g: [f"{g}_t"] for g in annotations}, coding={})
ref_alignments = {"g1_t": {"ABC1": (300.0, 1, 100, 100), "ABC2": (200.0, 1, 100, 100)},
                  "g2_t": {"ABC1": (300.0, 1, 100, 100)},
                  "g4_t": {"ABC2": (300.0, 1, 100, 100), "ABC1": (280.0, 1, 100, 100)}}
ref_descriptions = {"ABC1": "kinesin-like protein ABC1", "XYZ2": "zinc finger protein XYZ2", "ABC2": "OS9 endoplasmic reticulum lectin"}

with tempfile.TemporaryDirectory() as tmp:
    report = Path(tmp) / "eval.tsv"
    stats = evaluate.evaluate(annotations, universe, ref_alignments, {"ABC1", "ABC2", "XYZ2"}, ref_descriptions,
                              [FakeOrthology(), FakeInterPro()], report)
    rows = {r["gene_id"]: r for r in csv.DictReader(open(report), delimiter="\t")}
    check(set(rows) == {"g1", "g2", "g3", "g4", "g6"}, f"genes with hits or a symbol are evaluated: {sorted(rows)}")
    r = rows["g6"]
    check(r["decided_by_relatives"] == "yes" and r["hits_symbol"] == "minority" and r["product_symbol"] == "agree",
          f"prior override counted as decided; curated hit names carrying the symbol accept the product: {r['product_symbol']}")
    r = rows["g1"]
    check(r["hits_symbol"] == "agree" and r["reference"] == "best" and r["relatives"] == "agree"
          and r["same_species"] == "agree" and r["domains"] == "agree" and r["product_symbol"] == "agree"
          and r["conflict"] == "", f"all sources agree: {r}")
    r = rows["g2"]
    check(r["hits_symbol"] == "disagree" and r["reference"] == "absent" and r["reference_best_symbol"] == "ABC1"
          and r["relatives"] == "agree" and r["decided_by_relatives"] == "yes" and r["domains"] == "no_shared_word"
          and r["product_symbol"] == "disagree"
          and r["conflict"] == "hits_symbol,reference,product_symbol", f"taken symbol contradicted, domains never contradict: {r}")
    r = rows["g3"]
    check(r["hits_symbol"] == "none" and r["reference"] == "none" and r["hits_name"] == "agree" and r["conflict"] == "",
          f"no symbol: nothing to contradict: {r}")
    r = rows["g4"]
    check(r["reference"] == "near_best" and r["reference_ratio"] == "0.93" and r["relatives"] == "disagree"
          and r["decided_by_relatives"] == "no" and r["same_species"] == "agree" and r["hits_symbol"] == "agree"
          and r["conflict"] == "relatives", f"near-best reference, relatives disagree: {r}")
    check(stats["conflicts"] == 2 and stats["relatives_split"]["decided"] == {"agree": 1}
          and stats["relatives_split"]["independent"] == {"agree": 1, "disagree": 1},
          f"summary counts: {stats['conflicts']} {dict(stats['relatives_split'])}")
    lines = evaluate.summary_lines(stats)
    check(any(l.startswith("  Reference proteome") and "2" in l for l in lines) and any("Conflicts:" in l for l in lines),
          "summary lines")

print()
if failures:
    print(f"{len(failures)} check(s) failed:")
    for f in failures:
        print("  - " + f)
    sys.exit(1)
print("All checks passed")
