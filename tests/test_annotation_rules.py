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
g = QueryAnnotation(query_id="w", go_terms={"GO:1", "GO:2", "GO:3", "GO:4", "GO:5", "GO:6"},
                    go_evidence={"GO:1": "IEA", "GO:2": "IEA,IDA", "GO:3": "IEA", "GO:4": "ISS", "GO:5": "ISS",
                                 "GO:6": "TAS"},
                    go_support={"GO:1": 1, "GO:2": 1, "GO:3": 2, "GO:4": 2, "GO:5": 3, "GO:6": 1},
                    go_reference={"GO:6": ["P9"]})
st = annotate.filter_weak_go({"w": g}, None, min_hits=3, iea_min_hits=2)
check(g.go_terms == {"GO:2", "GO:5", "GO:6"} and set(g.go_weak_removed) == {"GO:1", "GO:3", "GO:4"}
      and st["terms"] == 3,
      f"weak terms from the hits alone removed; experimental and reference-confirmed kept: {sorted(g.go_terms)}")
g3 = QueryAnnotation(query_id="w3", go_terms={"GO:7", "GO:8"}, go_evidence={"GO:7": "ISS", "GO:8": "IBA"},
                    go_support={"GO:7": 1, "GO:8": 1}, go_curated={"GO:8"}, go_fallback={"GO:7": {"GO:9": 4}})
annotate.filter_weak_go({"w3": g3}, None, min_hits=3, iea_min_hits=2)
check(g3.go_terms == {"GO:8", "GO:9"} and g3.go_support["GO:9"] == 4,
      f"curated SwissProt terms of the gene kept; a weak term replaced by its supported ancestor: {sorted(g3.go_terms)}")
g4 = QueryAnnotation(query_id="w4", go_terms={"GO:10", "GO:11"}, go_evidence={"GO:10": "ISS", "GO:11": "IEA"},
                    go_support={"GO:10": 2, "GO:11": 1}, go_need=2)
g5 = QueryAnnotation(query_id="w5", go_terms={"GO:12", "GO:13"}, go_evidence={"GO:12": "IBA", "GO:13": "IEA"},
                    go_support={"GO:12": 1, "GO:13": 1}, go_need=1, go_single_strong=True)
annotate.filter_weak_go({"w4": g4, "w5": g5}, None, min_hits=5, iea_min_hits=2)
check(g4.go_terms == {"GO:10"} and g5.go_terms == {"GO:12"},
      f"a majority of few source hits suffices; one strong hit keeps its non-electronic terms: "
      f"{sorted(g4.go_terms)} {sorted(g5.go_terms)}")
g6 = QueryAnnotation(query_id="w6", go_terms={"GO:14", "GO:15"}, go_evidence={"GO:14": "ISS", "GO:15": "ISS"},
                    go_support={"GO:14": 1, "GO:15": 1}, go_need=2, go_ortholog_support={"GO:14": 1})
annotate.filter_weak_go({"w6": g6}, None, min_hits=5, iea_min_hits=2)
check(g6.go_terms == {"GO:14"}, f"a close relative's ortholog counts as one supporting source: {sorted(g6.go_terms)}")
g7 = QueryAnnotation(query_id="w7", go_terms={"GO:20", "GO:21", "GO:22", "GO:23"},
                    go_evidence={"GO:20": "IEA", "GO:21": "TAS", "GO:22": "IEA", "GO:23": "IEA"},
                    go_support={"GO:23": 5},
                    go_reference={"GO:20": ["P1"], "GO:21": ["P1"], "GO:22": ["P1"], "GO:23": ["P1"]},
                    go_ortholog_support={"GO:22": 1})
st7 = annotate.filter_weak_go({"w7": g7}, None, min_hits=5, iea_min_hits=2)
check(g7.go_terms == {"GO:21", "GO:22", "GO:23"} and g7.go_weak_removed["GO:20"][2] == "reference"
      and st7["reference_electronic"] == 1,
      f"electronic-only reference terms need the hits or a close relative's ortholog: {sorted(g7.go_terms)}")
g2 = QueryAnnotation(query_id="w2", go_terms={"GO:1"}, go_evidence={"GO:1": "IEA"}, go_support={"GO:1": 1})
annotate.filter_weak_go({"w2": g2}, None, min_hits=1, iea_min_hits=1)
check(g2.go_terms == {"GO:1"}, "minimum 1 keeps every term")
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
    # GO terms from the reference gene of the final symbol
    print("\n# Reference-gene GO terms")
    import sqlite3
    db = d / "ref.sqlite"
    conn = sqlite3.connect(db)
    conn.execute("CREATE TABLE go_mapping (accession TEXT, go_terms TEXT)")
    conn.executemany("INSERT INTO go_mapping VALUES (?, ?)",
                     [("P_AVP1", "GO:0000004:IDA;GO:0000001:IEA;GO:0000005:IDA"),      # lactation, its root-ish
                      ("P_AVP2", "GO:0000004:IEA;GO:0000003:IMP"),                      # ancestor, forbidden term
                      ("P_OXT", "GO:0000003:IDA")])
    conn.commit()
    conn.close()
    ref_symbols = {"P_AVP1": "AVP", "P_AVP2": "AVP", "P_OXT": "OXT", "P_X": "GRM1"}
    annots = {g: QueryAnnotation(query_id=g) for g in ("a", "b", "c", "e", "f")}
    annots["a"].consensus_symbol, annots["a"].symbol_status = "AVP", "rbh"            # gets lactation, root kept out
    annots["b"].consensus_symbol, annots["b"].symbol_status = "OXT", "ortholog"       # forbidden term only: nothing
    annots["c"].consensus_symbol, annots["c"].symbol_status = "AVP", "vote"           # vote only: not eligible
    annots["e"].consensus_symbol, annots["e"].symbol_status = "GRM1", "vote_supported"  # no GO on its proteins
    annots["f"].consensus_symbol, annots["f"].symbol_status = "AVP", "ortholog_fill"
    annots["f"].go_terms = {"GO:0000004"}
    annots["f"].specific_go_terms = {"GO:0000004"}
    annots["f"].go_evidence = {"GO:0000004": "IEA"}
    annots["f"].go_support = {"GO:0000004": 2}
    from goannotate import refgo
    annots["g"] = QueryAnnotation(query_id="g", consensus_symbol="avpb", symbol_status="ortholog", voted_symbol="AVP")
    forbidden = {"GO:0000003": "x", "GO:0000005": "y"}
    st = refgo.add_reference_go(annots, ref_symbols, str(db), hier, set(forbidden), [], "all", "add")
    check(annots["a"].specific_go_terms == {"GO:0000004"} and annots["a"].go_evidence["GO:0000004"] == "IDA,IEA"
          and annots["a"].go_reference["GO:0000004"] == ["P_AVP1", "P_AVP2"],
          f"terms of the reference gene's proteins, ancestors and forbidden terms left out: {annots['a'].specific_go_terms}")
    check(not annots["b"].go_terms and not annots["c"].go_terms and not annots["e"].go_terms,
          "forbidden-only, vote-only and GO-less symbols receive nothing")
    check(annots["f"].go_evidence["GO:0000004"] == "IEA" and annots["f"].go_reference["GO:0000004"]
          and annots["f"].specific_go_terms == {"GO:0000004"}, "a term the hits gave keeps its evidence, gains the source")
    check(st["genes_eligible"] == 5 and st["genes_added"] == 2 and annots["g"].specific_go_terms == {"GO:0000004"},
          f"counts {st}; a relatives' own form (avpb) reaches the reference gene through the replaced symbol")
    annots2 = {"a": QueryAnnotation(query_id="a", consensus_symbol="AVP", symbol_status="rbh",
                                    go_terms={"GO:0000002"}, specific_go_terms={"GO:0000002"})}
    st = refgo.add_reference_go(annots2, ref_symbols, str(db), hier, set(), [], "all", "fill")
    check(annots2["a"].specific_go_terms == {"GO:0000002"} and st["genes_added"] == 0, "fill leaves genes with GO alone")
    st = refgo.add_reference_go(annots2, ref_symbols, str(db), hier, set(), [], "all", "off")
    check(st["mode"] == "off" and st["genes_eligible"] == 0, "off")

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

# ---------------------------------------------------------- ortholog prior
print("\n# Ortholog symbols in the vote")


def hwn(symbol, bits, name=None, tx="t1"):
    """(hit, name, symbol, origin, source transcripts) as annotate_single_gene receives them."""
    h = hit(bits=bits, sid=f"sp|{symbol}{bits:.0f}|X_HUMAN")
    return (h, name or f"{symbol} protein", symbol, "uniprot", (tx,))


# Three hits say KCNQ4 (mean 300), two say KCNQ5 (mean 250): KCNQ4 wins the plain vote
pool = [hwn("KCNQ4", 320), hwn("KCNQ4", 300), hwn("KCNQ4", 280), hwn("KCNQ5", 260), hwn("KCNQ5", 240)]
idx = list(range(len(pool)))
name, sym, origin, _ = annotate.select_name_and_symbol(idx, pool, bnf)
check(sym == "KCNQ4" and origin == "uniprot", f"plain vote: {sym}")
name, sym, origin, _ = annotate.select_name_and_symbol(idx, pool, bnf, prior=("Kcnq5", True))
check(sym == "KCNQ5" and origin == "relative" and name == "KCNQ5 protein", f"close ortholog decides: {sym} ({origin})")
name, sym, origin, _ = annotate.select_name_and_symbol(idx, pool, bnf, prior=("Kcnq5", False))
check(sym == "KCNQ5" and origin == "uniprot", f"distant ortholog wins within the margin (250/300): {sym}")
weak = pool + [hwn("KCNQ3", 100)]
name, sym, origin, _ = annotate.select_name_and_symbol(list(range(len(weak))), weak, bnf, prior=("KCNQ3", False))
check(sym == "KCNQ4", f"distant ortholog below the margin does not win: {sym}")
name, sym, origin, _ = annotate.select_name_and_symbol(idx, pool, bnf, prior=("OXT", True))
check(sym == "KCNQ4", f"a prior symbol absent from the hits changes nothing here: {sym}")

# The close prior also picks the cluster: two clusters, the prior's symbol only in the weaker one
annotate._init_worker({}, None, bnf, mappings.SymbolResolver())
far = [hwn("ACTB", 500, "Actin, cytoplasmic 1"), hwn("ACTB", 480, "Actin, cytoplasmic 1"),
       hwn("OXT", 200, "Oxytocin-neurophysin 1")]
gid, res = annotate.annotate_single_gene("g", far, 0.5, [], "upper", "mean", prior=("Oxt", True))
check(res["consensus_symbol"] == "OXT" and res["symbol_source"] == "relative"
      and res["consensus_name"] == "Oxytocin-neurophysin 1", f"close prior selects the cluster: {res['consensus_symbol']}")
gid, res = annotate.annotate_single_gene("g", far, 0.5, [], "upper", "mean", prior=("Oxt", False))
check(res["consensus_symbol"] == "ACTB", f"distant prior does not change the cluster: {res['consensus_symbol']}")
gid, res = annotate.annotate_single_gene("g", far, 0.5, [], "upper", "mean")
check(res["consensus_symbol"] == "ACTB", "no prior")

# ------------------------------------------------- excluded and weak reference symbols
print("\n# Pseudogenes, readthroughs and weakly evidenced paralogs")
from goannotate import assign as assign_ref  # noqa: E402
nom = mappings.ReferenceNomenclature()
for sym in ("NANOG", "NANOGP8", "IQCJ", "SCHIP1", "IQCJ-SCHIP1", "HLA-A", "CDHR18P"):
    nom.official[sym.lower()] = sym
nom.pseudogenes.add("cdhr18p")
check(nom.is_readthrough("IQCJ-SCHIP1") and not nom.is_readthrough("HLA-A") and nom.is_excluded("CDHR18P")
      and not nom.is_excluded("NANOGP8"), "readthrough and pseudogene recognition")
res = mappings.SymbolResolver()
res.reference_symbols = frozenset(nom.official)
res.reference_pseudogenes = frozenset(nom.pseudogenes)
pool = [hwn("IQCJ-SCHIP1", 330), hwn("IQCJ-SCHIP1", 320), hwn("SCHIP1", 300), hwn("SCHIP1", 290)]
name, sym, origin, _ = annotate.select_name_and_symbol(list(range(4)), pool, bnf, resolver=res)
check(sym == "SCHIP1", f"a readthrough group does not win the vote: {sym}")
nom.official["uox"] = "UOX"
nom.pseudogenes.add("uox")
res.reference_symbols = frozenset(nom.official)
res.reference_pseudogenes = frozenset(nom.pseudogenes)
mouse = hit(bits=400, sid="sp|P_UOX|URIC_MOUSE", title="sp|P_UOX|URIC_MOUSE Uricase OS=Mus musculus OX=10090 GN=Uox PE=1 SV=1")
human = hit(bits=400, sid="sp|P_U2|X_HUMAN", title="sp|P_U2|X_HUMAN Uricase OS=Homo sapiens OX=9606 GN=UOX PE=5 SV=1")
pool = [(mouse, "Uricase", "Uox", "uniprot", ("t1",)), (mouse, "Uricase", "Uox", "uniprot", ("t1",)),
        hwn("ABC1", 200), hwn("ABC1", 190)]
name, sym, origin, _ = annotate.select_name_and_symbol(list(range(4)), pool, bnf, resolver=res)
check(sym == "Uox", f"a pseudogene symbol carried by other species' hits names their functional gene: {sym}")
pool = [(human, "Uricase", "UOX", "uniprot", ("t1",)), (human, "Uricase", "UOX", "uniprot", ("t1",)),
        hwn("ABC1", 200), hwn("ABC1", 190)]
name, sym, origin, _ = annotate.select_name_and_symbol(list(range(4)), pool, bnf, resolver=res)
check(sym == "ABC1", f"a pseudogene symbol carried only by reference-species hits does not win: {sym}")
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "fwd.tsv").write_text("t1\tsp|P_NP8|X\t95\t300\t0\t400\t300\t300\n"
                               "t1\tsp|P_NANOG|X\t94\t300\t0\t390\t300\t305\n"
                               "t2\tsp|P_IQ|X\t90\t400\t0\t500\t400\t420\n"
                               "t2\tsp|P_SCH|X\t90\t390\t0\t480\t400\t400\n")
    (d / "rev.tsv").write_text("sp|P_NP8|X\tt1\t95\t300\t0\t400\t300\t300\n"
                               "sp|P_NANOG|X\tt1\t94\t300\t0\t390\t305\t300\n"
                               "sp|P_IQ|X\tt2\t90\t400\t0\t500\t420\t400\n"
                               "sp|P_SCH|X\tt2\t90\t390\t0\t480\t400\t400\n")
    ref = {"P_NP8": "NANOGP8", "P_NANOG": "NANOG", "P_IQ": "IQCJ-SCHIP1", "P_SCH": "SCHIP1"}
    uni = type("U", (), {"transcript_to_gene": {"t1": "g1", "t2": "g2"}})()
    ann = {"g1": QueryAnnotation(query_id="g1", consensus_symbol="NANOGP8", hits=[hit()]),
           "g2": QueryAnnotation(query_id="g2", hits=[hit()])}
    ann["g1"].symbol_status = "vote"
    st = assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", ref, str.upper,
                                   ref_descriptions={"NANOG": "Homeobox protein NANOG"}, nomenclature=nom,
                                   existence={"NANOGP8": 2, "NANOG": 1, "IQCJ-SCHIP1": 1, "SCHIP1": 1})
    check(ann["g1"].consensus_symbol == "NANOG" and ann["g1"].voted_symbol == "NANOGP8"
          and ann["g1"].consensus_name == "Homeobox protein NANOG" and st["established_preferred"] == 1,
          f"established gene preferred over a weakly evidenced near-identical paralog: {ann['g1'].consensus_symbol}")
    check(ann["g2"].consensus_symbol == "SCHIP1" and ann["g2"].symbol_status == "rbh_added",
          f"a readthrough is never the reciprocal best hit: {ann['g2'].consensus_symbol}")
    ann["g1"].consensus_symbol, ann["g1"].voted_symbol, ann["g1"].symbol_status = "NANOGP8", "", "vote"
    assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", ref, str.upper, nomenclature=nom,
                              existence={"NANOGP8": 2, "NANOG": 1}, established_margin=0.99)
    check(ann["g1"].consensus_symbol == "NANOGP8", "not replaced when the established gene scores below the margin")

# ------------------------------------------------- copy forms and propagation
print("\n# Copy forms of the RBH symbol and symbol propagation")
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "fwd.tsv").write_text("t1\tsp|P_ADG|X\t90\t400\t0\t500\t400\t400\n"
                               "t2\tsp|P_AG1|X\t90\t400\t0\t500\t400\t400\n")
    (d / "rev.tsv").write_text("sp|P_ADG|X\tt1\t90\t400\t0\t500\t400\t400\n"
                               "sp|P_AG1|X\tt2\t90\t400\t0\t500\t400\t400\n")
    ref = {"P_ADG": "ADGRL1", "P_AG1": "AGPAT1", "P_AG2": "AGPAT2"}
    uni = type("U", (), {"transcript_to_gene": {"t1": "g1", "t2": "g2"}})()
    ann = {"g1": QueryAnnotation(query_id="g1", consensus_symbol="adgrl1a", hits=[hit()]),
           "g2": QueryAnnotation(query_id="g2", consensus_symbol="AGPAT2", hits=[hit()])}
    ann["g2"].winning_cluster_symbols = ["AGPAT2", "AGPAT1"]
    st = assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", ref, str.upper,
                                   ref_descriptions={"ADGRL1": "Adhesion G protein-coupled receptor L1"})
    check(ann["g1"].symbol_status == "rbh" and ann["g1"].consensus_symbol == "ADGRL1"
          and ann["g1"].voted_symbol == "adgrl1a" and st["rbh_same_gene_form"] == 1 and st["rbh_override"] == 1,
          f"a copy form of the RBH symbol is a confirmation, another gene an override: "
          f"{ann['g1'].symbol_status} {st['rbh_same_gene_form']} {st['rbh_override']}")
sib = lambda sym: (hit(), "Glutamine synthetase", sym, "uniprot", ("t1",))
ann = {"a": QueryAnnotation(query_id="a", consensus_name="Glutamine synthetase", consensus_symbol="glula"),
       "b": QueryAnnotation(query_id="b", consensus_name="Glutamine synthetase", consensus_symbol="glula"),
       "c": QueryAnnotation(query_id="c", consensus_name="Glutamine synthetase"),
       "d": QueryAnnotation(query_id="d", consensus_name="Glutamine synthetase")}
n, ids = annotate.propagate_symbols(ann, {"c": [sib("GLUL"), sib("gs")], "d": [sib("CA2")]}, "lower")
check(n == 1 and ann["c"].consensus_symbol == "glul" and not ann["d"].consensus_symbol,
      f"propagation gives the hits' reference form of a relative's copy form, never a symbol the hits lack: "
      f"{ann['c'].consensus_symbol!r} {ann['d'].consensus_symbol!r}")
ann["c"].consensus_symbol = ""
ann["a"].consensus_symbol = "glul"
n, ids = annotate.propagate_symbols(ann, {"c": [sib("GLUL")]}, "lower")
check(ann["c"].consensus_symbol == "glul", "a sibling symbol the hits carry is taken as it is")

# ------------------------------------------------------ confidence with orthology
print("\n# Confidence with orthology evidence")
a = QueryAnnotation(query_id="c", hits=[hit()], consensus_symbol="X", name_concordance=0.6, n_cluster_hits=5,
                    symbol_concordance=0.5)
base = annotate.compute_confidence(a)[1]
check(abs(annotate.compute_confidence(a, "close_confirmed")[1] - min(1.0, base + 0.15)) < 1e-9
      and abs(annotate.compute_confidence(a, "conflict")[1] - (base - 0.15)) < 1e-9,
      f"confirmation adds, a kept conflict subtracts (base {base})")
check(annotate.compute_confidence(a, "close_taken")[0] == "medium", "taken symbols are medium")
check(annotate.compute_confidence(QueryAnnotation(query_id="n"), "close_filled") == ("medium", 0.55)
      and annotate.compute_confidence(QueryAnnotation(query_id="n")) == ("none", 0.0), "a fill gives a gene without hits a score")

# ------------------------------------------------------ reference names of confirmed symbols
print("\n# Reference names and clone names")
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "fwd.tsv").write_text("t1\tsp|P_MND|X\t90\t400\t0\t500\t400\t400\n"
                               "t2\tsp|P_PAK|X\t90\t400\t0\t500\t400\t400\n")
    (d / "rev.tsv").write_text("sp|P_MND|X\tt1\t90\t400\t0\t500\t400\t400\n"
                               "sp|P_PAK|X\tt2\t90\t400\t0\t500\t400\t400\n")
    ref = {"P_MND": "MND", "P_PAK": "PAK1"}
    uni = type("U", (), {"transcript_to_gene": {"t1": "g1", "t2": "g2"}})()
    ann = {"g1": QueryAnnotation(query_id="g1", consensus_symbol="MND", hits=[hit()],
                                 consensus_name="Y+L amino acid transporter 2",
                                 hit_names={"y+l amino acid transporter 2": 9, "minidiscs": 1}),
           "g2": QueryAnnotation(query_id="g2", consensus_symbol="PAK1", hits=[hit()],
                                 consensus_name="Serine/threonine-protein kinase",
                                 hit_names={"serine/threonine-protein kinase pak 1": 6, "serine/threonine-protein kinase": 4})}
    st = assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", ref, str.upper,
                                   ref_descriptions={"MND": "Minidiscs", "PAK1": "Serine/threonine-protein kinase PAK 1"})
    check(ann["g1"].consensus_name == "Minidiscs" and ann["g2"].consensus_name == "Serine/threonine-protein kinase PAK 1",
          f"without curated names of other species the reference name is used: {ann['g1'].consensus_name}")
    ann = {"g1": QueryAnnotation(query_id="g1", consensus_symbol="MND", hits=[hit()],
                                 consensus_name="Y+L amino acid transporter 2",
                                 hit_names={"y+l amino acid transporter 2": 9, "minidiscs": 1},
                                 hit_curated_names={"y+l amino acid transporter 2": ["SLC7A6"]}),
           "g2": QueryAnnotation(query_id="g2", consensus_symbol="PAK1", hits=[hit()], consensus_name="RNA helicase",
                                 hit_names={"rna helicase": 30, "eukaryotic initiation factor 4a-iii": 2},
                                 hit_curated_names={"eukaryotic initiation factor 4a-iii": ["EIF4A3"]}),
           "g3": QueryAnnotation(query_id="g3", consensus_symbol="KLK13", hits=[hit()], consensus_name="Chymotrypsin",
                                 hit_names={"chymotrypsin": 30}, hit_curated_names={"kallikrein 1-related peptidase b3": ["KLK1B3"]})}
    (d / "fwd.tsv").write_text((d / "fwd.tsv").read_text() + "t3\tsp|P_KLK|X\t90\t400\t0\t500\t400\t400\n")
    (d / "rev.tsv").write_text((d / "rev.tsv").read_text() + "sp|P_KLK|X\tt3\t90\t400\t0\t500\t400\t400\n")
    ref["P_KLK"] = "KLK13"
    uni.transcript_to_gene["t3"] = "g3"
    assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", ref, str.upper,
                              ref_descriptions={"MND": "Minidiscs", "PAK1": "Eukaryotic initiation factor 4A-III",
                                                "KLK13": "Kallikrein-13"})
    check(ann["g1"].consensus_name == "Y+L amino acid transporter 2"
          and ann["g2"].consensus_name == "Eukaryotic initiation factor 4A-III"
          and ann["g3"].consensus_name == "Kallikrein-13",
          f"a curated name of other species decides: {ann['g1'].consensus_name} / {ann['g2'].consensus_name}")
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "fwd.tsv").write_text("t1\tsp|P_U3|X\t90\t400\t0\t500\t400\t400\nt1\tsp|P_U2|X\t80\t400\t0\t450\t400\t400\n")
    (d / "rev.tsv").write_text("sp|P_U3|X\tt1\t90\t400\t0\t500\t400\t400\n")
    uni = type("U", (), {"transcript_to_gene": {"t1": "g1"}})()
    ann = {"g1": QueryAnnotation(query_id="g1", consensus_symbol="UCP3", hits=[hit()], consensus_name="Dicarboxylate carrier UCP2",
                                 hit_names={"dicarboxylate carrier ucp2": 9},
                                 hit_curated_names={"dicarboxylate carrier ucp2": ["UCP2"]})}
    assign_ref.assign_symbols(ann, uni, d / "fwd.tsv", d / "rev.tsv", {"P_U3": "UCP3", "P_U2": "UCP2"}, str.upper,
                              ref_descriptions={"UCP3": "Mitochondrial uncoupling protein 3"})
    check(ann["g1"].consensus_name == "Mitochondrial uncoupling protein 3",
          f"curated entries of a paralog do not keep the hits' name: {ann['g1'].consensus_name}")
g = QueryAnnotation(query_id="s", consensus_symbol="TRET1", consensus_name="Nebulosa", secondary_symbol="SLC2A8",
                    hit_names={"facilitated trehalose transporter tret1": 5, "solute carrier family 2, facilitated glucose transporter member 8": 2})
g2 = QueryAnnotation(query_id="s2", consensus_symbol="MND", consensus_name="Minidiscs", secondary_symbol="SLC7A6",
                     hit_names={"minidiscs": 1, "y+l amino acid transporter 2": 5})
r = assign_ref.secondary_names({"s": g, "s2": g2},
                               {"SLC2A8": "Solute carrier family 2, facilitated glucose transporter member 8",
                                "SLC7A6": "Y+L amino acid transporter 2"},
                               {"TRET1": "Nebulosa", "MND": "Minidiscs"})
check(g.consensus_name.startswith("Solute carrier family 2") and g.product_note == "secondary_reference"
      and g2.consensus_name == "Minidiscs" and r["secondary_names"] == 1,
      f"secondary reference name when no hit uses the primary one: {g.consensus_name} / {g2.consensus_name}")
check(not bnf.clean_name("Protein CG4768-like protein") or bnf.is_bad_name(bnf.clean_name("Protein CG4768-like protein")),
      "an identifier-only -like name leaves a placeholder")
check(bnf.is_bad_name("protein -like protein") and bnf.is_bad_name("Protein, protein")
      and not bnf.is_bad_name("Ras-like protein"), "'protein like protein' is a placeholder")
check(all(bnf.is_bad_name(n) for n in ("Chromosome 9", "Chromosome LG6 open reading frame, human", "CD033",
                                       "CQ067 protein precursor", "protein (Mus musculus)", "Linkage group 7"))
      and not bnf.is_bad_name("Chromosome partition protein Smc"), "chromosome-only and species-only names are bad")
check(bnf.clean_name("LD33793p") in (None, "") and bnf.clean_name("FI18136p1") in (None, "")
      and bnf.clean_name("GEO02646p1") in (None, "") and bnf.clean_name("bcDNA.GH11023") in (None, "") and bnf.clean_name("Ubiquitin-like protein RH60221p") == "Ubiquitin-like protein",
      "BDGP cDNA clone names are identifiers")

# ------------------------------------------------------ final concordance
print("\n# Concordance with the final annotation")
fh = lambda: hit()
a = QueryAnnotation(query_id="k", hits=[fh()], consensus_symbol="kdm2aa", voted_symbol="KDM2A",
                    consensus_name="Lysine-specific demethylase 2A", n_cluster_hits=10,
                    hit_symbols={"KDM2A": 6, "KDM2AA": 2, "KDM2B": 2}, hit_names={"lysine specific demethylase 2a": 8,
                                                                                "f-box protein": 2})
t = QueryAnnotation(query_id="t", hits=[fh()], consensus_symbol="agpat2", voted_symbol="AGPAT1",
                    symbol_aliases=["AGPAT1"], consensus_name="1-acyl-sn-glycerol-3-phosphate acyltransferase beta",
                    hit_symbols={"AGPAT1": 9, "AGPAT2": 1}, ortholog_evidence="close_taken")
n = QueryAnnotation(query_id="n", consensus_symbol="abc1", ortholog_evidence="close_filled")
st = annotate.final_concordance({"k": a, "t": t, "n": n}, consensus_threshold=0.5)
check(abs(a.symbol_concordance - 0.8) < 1e-9 and abs(a.name_concordance - 0.8) < 1e-9,
      f"a refined copy form counts the hits of its gene: {a.symbol_concordance} {a.name_concordance}")
check(abs(t.symbol_concordance - 0.1) < 1e-9, f"a taken symbol counts only hits of its own gene: {t.symbol_concordance}")
check(st["symbol"]["n"] == 2 and n.confidence == "medium", f"genes without hits are left out: {st['symbol']}")

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
check(annots["g4"].consensus_name == "Forkhead box protein O3", "a hits' name kept on a confirmed gene is stripped too")
check(annots["g5"].consensus_name == "Olfactory receptor 2", "2b -> 2 when the base name is among the hits")
check(r["designators_stripped"] == 4, f"4 products changed ({r})")

print()
if failures:
    print(f"{len(failures)} check(s) failed")
    sys.exit(1)
print("All checks passed")
