#!/usr/bin/env python3
"""
Checks of the annotation rules on small inputs (no databases needed):
hit coverage filter, word weighting, GO evidence policy and rules, taxon
constraints, product naming rules, symbol policy and paralog designators.

Run: python3 tests/test_annotation_rules.py
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import annotate, blast, mappings, output, taxon  # noqa: E402
from goannotate.filters import BadNameFilter, default_bad_names_path  # noqa: E402
from goannotate.go import GOHierarchy  # noqa: E402
from goannotate.models import BlastHit, QueryAnnotation  # noqa: E402

failures = []


def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


def hit(pident=60.0, length=200, qstart=1, qend=200, sstart=1, send=200, qlen=200, slen=200, bits=300.0,
        sid="sp|P1|X_HUMAN", title="sp|P1|X_HUMAN Protein X OS=Homo sapiens OX=9606 GN=X PE=1 SV=1"):
    return BlastHit("q1", sid, pident, length, 0, 0, qstart, qend, sstart, send, 1e-50, bits, title, qlen, slen)


# ---------------------------------------------------------------- coverage
print("\n# Hit coverage filter")
full = hit()
domain = hit(qstart=1, qend=120, sstart=400, send=520, qlen=120, slen=1200)   # whole query, 10% of subject
frag = hit(qstart=1, qend=150, sstart=300, send=450, qlen=150, slen=600)      # fragment: 100% q, 25% s
lowid = hit(pident=30.0)
nolen = hit(qlen=0, slen=0)
strict = blast.HitCoverageFilter("strict")
perm = blast.HitCoverageFilter("permissive")
off = blast.HitCoverageFilter("off")
check(abs(domain.subject_coverage - 121 / 1200) < 1e-9 and domain.query_coverage == 1.0, "coverage properties")
check(strict.passes(full) and not strict.passes(domain) and not strict.passes(frag), "strict mode")
check(perm.passes(full) and perm.passes(domain) and perm.passes(frag), "permissive mode")
check(not strict.passes(lowid) and not perm.passes(lowid) and off.passes(lowid), "identity floor")
check(strict.passes(nolen), "hits without lengths pass on identity alone")

bnf = BadNameFilter(str(default_bad_names_path()))
hits = {"t1": [domain, full], "t2": [domain]}
filtered, status = blast.filter_and_select_hits(hits, bnf, top_n=30, coverage=strict)
check(len(filtered["t1"]) == 1 and filtered["t1"][0][0] is full and status["t1"] == "pass",
      "domain-only hit removed when a full hit exists")
check(len(filtered["t2"]) == 1 and status["t2"] == "fallback", "fallback keeps the unfiltered hits")
filtered, status = blast.filter_and_select_hits(hits, bnf, top_n=30,
                                                coverage=blast.HitCoverageFilter("strict", fallback=False))
check("t2" not in filtered, "no fallback drops the transcript")

# ----------------------------------------------------------- blast parsing
print("\n# Results parsing, both layouts")
with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
    fh.write("q1\tsp|P1|X_HUMAN\t60.0\t200\t0\t0\t1\t200\t1\t200\t1e-50\t300\t200\t250\tsp|P1|X_HUMAN Protein X OS=Homo sapiens OX=9606 GN=X PE=1 SV=1\n")
    fh.write("q2\tsp|P2|Y_HUMAN\t60.0\t200\t0\t0\t1\t200\t1\t200\t1e-50\t300\tsp|P2|Y_HUMAN Protein Y OS=Homo sapiens OX=9606 GN=Y PE=1 SV=1\n")
    path = fh.name
parsed, stats = blast.parse_blast_results(path)
check(parsed["q1"][0].qlen == 200 and parsed["q1"][0].slen == 250
      and bnf.clean_name(parsed["q1"][0].gene_name) == "Protein X" and parsed["q1"][0].gene_symbol == "X",
      "current layout with qlen/slen")
check(parsed["q2"][0].qlen == 0 and bnf.clean_name(parsed["q2"][0].gene_name) == "Protein Y"
      and stats.hits_with_lengths == 1, "1.x layout without lengths")
Path(path).unlink()

# ------------------------------------------------------------ word weights
print("\n# Word weighting")
names = ["protein kinase c alpha", "protein kinase c beta", "olfactory receptor 4k36", "olfactory receptor 2a1",
         "zinc finger protein 91", "forkhead box protein o3", "protein fam83a"]
ww = annotate.WordWeights("idf", annotate.load_stop_words()).fit(names)
check(ww.weight("protein") == 0.0 and ww.weight("4k36") > ww.weight("receptor") > 0, "stop words 0, rare > common")
annotate._worker_word_weights = None
plain = annotate.name_similarity("protein fam83a", "protein fam83b")
annotate._worker_word_weights = ww
weighted = annotate.name_similarity("protein fam83a", "protein fam83b")
check(plain >= 0.33 and weighted < plain or weighted == plain, f"weighted Jaccard does not exceed plain ({plain:.2f} -> {weighted:.2f})")
j_plain = len({"olfactory", "receptor", "4k36"} & {"olfactory", "receptor", "2a1"}) / 4
annotate._worker_word_weights = None
check(annotate.name_similarity("protein kinase c alpha", "protein kinase c alpha") == 1.0, "identical names 1.0")
annotate._worker_word_weights = None

# -------------------------------------------------------------- GO evidence
print("\n# GO evidence")
pairs = mappings.parse_go_field("GO:0005515:IDA,IPI;GO:0005634:IEA;GO:0000001")
check(pairs == [("GO:0005515", "IDA,IPI"), ("GO:0005634", "IEA"), ("GO:0000001", "")], "parse_go_field with codes")
check(mappings.parse_go_field("GO:0000001,GO:0000002") == [("GO:0000001", ""), ("GO:0000002", "")],
      "parse_go_field 1.x comma layout")
for policy, expect in (("all", {"GO:0005515", "GO:0005634"}), ("non-electronic", {"GO:0005515"}),
                       ("experimental", {"GO:0005515"})):
    m = {}
    st = {"excluded": 0, "nd": 0, "no_codes": 0, "with_codes": 0}
    mappings._add_terms(m, "P1", "GO:0005515:IDA,IPI;GO:0005634:IEA;GO:0008150:ND", mappings._allowed_codes(policy), st)
    check(set(m["P1"]) == expect and st["nd"] == 1, f"policy {policy}")
codes = {"GO:1": {"IEA"}, "GO:2": {"IEA", "IDA"}, "GO:3": {"IEA"}, "GO:0008150": {"IDA"}}
support = {"GO:1": 1, "GO:2": 1, "GO:3": 3}
kept = annotate._apply_go_evidence_rules(codes, support, iea_min_hits=2)
check(kept == {"GO:2", "GO:3"}, "IEA-only terms need the minimum support; roots never assigned")
check(annotate._apply_go_evidence_rules(codes, support, 1) == {"GO:1", "GO:2", "GO:3"}, "iea_min_hits 1 keeps all")
tc, ts = mappings.collect_go_evidence([hit()], {"P1": {"GO:1": "IEA,IDA", "GO:2": ""}})
check(tc == {"GO:1": {"IEA", "IDA"}, "GO:2": set()} and ts == {"GO:1": 1, "GO:2": 1}, "collect_go_evidence")

# ----------------------------------------------------------- taxon filter
print("\n# Taxon constraints")
obo = """format-version: 1.2

[Term]
id: GO:0000001
name: root-ish
namespace: biological_process

[Term]
id: GO:0000002
name: photosynthesis
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000003
name: photosystem part
namespace: biological_process
relationship: part_of GO:0000002

[Term]
id: GO:0000004
name: lactation
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000005
name: spore wall
namespace: cellular_component
"""
tc_obo = """format-version: 1.2

[Term]
id: GO:0000002
relationship: RO:0002162 NCBITaxon_Union:0000002

[Term]
id: GO:0000004
relationship: RO:0002162 NCBITaxon:40674

[Term]
id: GO:0000005
property_value: RO:0002161 NCBITaxon:40674
"""
tg_obo = """format-version: 1.2

[Term]
id: NCBITaxon_Union:0000002
name: Viridiplantae or Cyanobacteria
union_of: NCBITaxon:33090
union_of: NCBITaxon:1117
"""
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "go.obo").write_text(obo)
    (d / "tc.obo").write_text(tc_obo)
    (d / "tg.obo").write_text(tg_obo)
    hier = GOHierarchy(str(d / "go.obo"))
    cons = taxon.parse_taxon_constraints_obo(d / "tc.obo")
    unions = taxon.parse_taxon_groupings_obo(d / "tg.obo")
    check(cons["GO:0000002"] == [("in_taxon", "NCBITaxon_Union:0000002")]
          and cons["GO:0000005"] == [("never_in_taxon", "NCBITaxon:40674")], "constraint parsing")
    check(taxon.expand_taxon("NCBITaxon_Union:0000002", unions) == {33090, 1117}, "union expansion")
    vole_lineage = [79684, 10053, 9989, 40674, 7742, 33208, 2759, 131567, 1]
    forb = taxon.forbidden_go_terms(cons, unions, vole_lineage, hier)
    check(set(forb) == {"GO:0000002", "GO:0000003", "GO:0000005"},
          f"vole: photosynthesis (and its part) and spore wall excluded, lactation kept: {sorted(forb)}")
    plant_lineage = [3702, 33090, 2759, 131567, 1]
    forb = taxon.forbidden_go_terms(cons, unions, plant_lineage, hier)
    check(set(forb) == {"GO:0000004"}, f"plant: lactation excluded only: {sorted(forb)}")
    (d / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2759\t|\t1\t|\tdomain\t|\n40674\t|\t2759\t|\tclass\t|\n"
                                 "79684\t|\t40674\t|\tspecies\t|\n")
    (d / "merged.dmp").write_text("99999\t|\t79684\t|\n")
    check(taxon.lineage_from_taxdump(d, 99999) == [79684, 40674, 2759, 1], "lineage from taxdump with merged ID")

# ----------------------------------------------------------- product rules
print("\n# Product naming rules")
output.configure_products(bnf, cleanup=True)
cases = {
    "Hexokinase-1 (EC 2.7.1.1)": "hexokinase-1",
    "Heat shock 70 kDa protein 1A": "heat shock 70 kDa protein 1A",
    "Glutathione S-transferase (45 kDa)": "glutathione S-transferase",
    "Protein wntless homolog (Drosophila homolog)": "protein wntless homolog",
    "Disks large homolog 5 (Homo sapiens)": "disks large homolog 5",
    "FAM83A protein": "protein FAM83A",
    "Kinesin protein": "kinesin protein",
    "Annexin A5 precursor": "annexin A5",
    "Ras-related protein Rab-11b, tandem duplicate 1": "ras-related protein Rab-11b, tandem duplicate 1",
    "Cytochrome b": "cytochrome b",
    "Uncharacterized protein protein": "uncharacterized protein",
}
for raw, expect in cases.items():
    got = output.sanitize_product(raw)
    check(got == expect, f"{raw!r} -> {got!r}")
check(bnf.clean_name("Ras-related protein Rab-11b, tandem duplicate 1") == "Ras-related protein Rab-11b",
      "tandem duplicate stripped from hit names (qualifier)")
check(bnf.process_name("Not available") is None and bnf.process_name("NA") is None, "'not available' placeholder")

# ------------------------------------------------------ paralog designators
print("\n# Paralog designators")
norm = bnf.normalize_for_comparison
ref = {"FOXO3": "Forkhead box protein O3", "PRKCB": "Protein kinase C beta type", "CYB": "Cytochrome b"}
annots = {}
for i, (name, status, pool) in enumerate([
        ("Forkhead box protein O3-B", "vote", ["Forkhead box protein O3-B"]),
        ("Protein kinase C beta type a", "vote", []),
        ("Cytochrome b", "vote", ["Cytochrome b"]),
        ("Annexin A", "vote", ["Annexin A"]),
        ("Forkhead box protein O3-B", "rbh", []),
        ("Olfactory receptor 2b", "vote", ["Olfactory receptor 2b", "Olfactory receptor 2"])]):
    a = QueryAnnotation(query_id=f"g{i}")
    a.consensus_name = name
    a.symbol_status = status
    a.winning_cluster_names = pool
    annots[a.query_id] = a
assign_mod = __import__("goannotate.assign", fromlist=["x"])
r = assign_mod.strip_paralog_designators(annots, ref, norm)
check(annots["g0"].consensus_name == "Forkhead box protein O3" and annots["g0"].product_note == "designator_stripped",
      "O3-B -> O3 (reference name)")
check(annots["g1"].consensus_name == "Protein kinase C beta type", "trailing ' a' stripped (reference name)")
check(annots["g2"].consensus_name == "Cytochrome b", "Cytochrome b untouched")
check(annots["g3"].consensus_name == "Annexin A", "Annexin A untouched (uppercase)")
check(annots["g4"].consensus_name == "Forkhead box protein O3-B", "reference-derived products untouched")
check(annots["g5"].consensus_name == "Olfactory receptor 2", "2b -> 2 when the base name is among the hits")
check(r["designators_stripped"] == 3, f"3 products changed ({r})")

print()
if failures:
    print(f"{len(failures)} check(s) failed")
    sys.exit(1)
print("All checks passed")
