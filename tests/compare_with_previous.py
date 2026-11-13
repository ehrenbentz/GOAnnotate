#!/usr/bin/env python3
"""
Compare a GOAnnotate run with an earlier functional annotation of the same gene models.

Previous symbols come from a GFF3 attribute (default Name on gene features);
previous products and GO terms from an optional table with a gene ID column,
a description column (UniProt-style titles are cleaned of OS=/OX=/GN= fields)
and a GO column (IDs separated by , ; or |).

Reports:
  - previous symbols by class: none ("NONE", empty), gene-model IDs, genome-
    specific identifiers (by the GOAnnotate symbol patterns), informative
  - symbol agreement (case-insensitive; "related" also counts fish paralog
    suffixes a/b/l and species case differences such as foxo3a vs FOXO3)
  - product agreement (GOAnnotate name similarity >= 0.8 match, >= 0.5 partial)
  - GO coverage and terms per gene; for genes annotated in both, the overlap of
    the two term sets after adding all ancestors (Jaccard)
Writes {out}/reports/{prefix}_previous_comparison.tsv with one row per gene.

Usage:
    compare_with_previous.py OUTPUT_DIR --previous-gff3 OLD.gff3 [--symbol-attr Name]
        [--previous-table OLD.tsv --id-col Name --desc-col Protein_Description --go-col Annotation_GO_ID]
        [--go-obo go.obo] [--prefix PREFIX]
"""
import argparse
import csv
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.annotate import name_similarity, normalize_gene_name  # noqa: E402
from goannotate.filters import BadNameFilter, default_bad_patterns_path  # noqa: E402
from goannotate.gff3 import parse_attributes  # noqa: E402

csv.field_size_limit(sys.maxsize)

MODEL_ID = re.compile(r"^(?:anno\d+\.)?g\d+(?:\.t\d+)?$|^gene\d+$|^mRNA\d+$|^jg\d+", re.I)


def clean_title(desc):
    if not desc:
        return ""
    s = re.sub(r"\s+(?:OS|OX|GN|PE|SV)=.*$", "", desc.strip())
    return s


def symbol_class(sym, bnf):
    if not sym or sym.upper() in ("NONE", "NA", "-", "."):
        return "none"
    if MODEL_ID.match(sym):
        return "gene_model_id"
    if bnf.is_bad_symbol(sym):
        return "genome_specific_identifier"
    return "informative"


def _stem(sym):
    """Symbol without teleost paralog suffixes (a, b, l, la, lb, l2) and punctuation."""
    s = sym.upper().replace(".", "").replace("-", "")
    return re.sub(r"(?:L?[AB]|L\d*)$", "", s) or s


def related(old, new):
    """
    Same gene allowing teleost paralog suffixes on either side (myh9 / MYH9A,
    PCDH7 / PCDH7A), species casing, and old symbols given as comma-separated
    alternatives (ARHGEF1,ARHGEF12).
    """
    n = new.upper()
    for o in old.split(","):
        o = o.strip().upper()
        if not o:
            continue
        if o == n or _stem(o) == _stem(n):
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("output_dir")
    ap.add_argument("--prefix")
    ap.add_argument("--previous-gff3", required=True)
    ap.add_argument("--symbol-attr", default="Name")
    ap.add_argument("--previous-table")
    ap.add_argument("--id-col", default="Name")
    ap.add_argument("--desc-col", default="Protein_Description")
    ap.add_argument("--go-col", default="Annotation_GO_ID")
    ap.add_argument("--go-obo")
    args = ap.parse_args()

    out = Path(args.output_dir)
    prefix = args.prefix or out.resolve().name
    bnf = BadNameFilter(str(default_bad_patterns_path()))

    old_sym = {}
    with open(args.previous_gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) == 9 and p[2] == "gene":
                a = parse_attributes(p[8])
                old_sym[a.get("ID", "")] = a.get(args.symbol_attr, "")

    old_desc, old_go = {}, {}
    if args.previous_table:
        with open(args.previous_table) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gid = r.get(args.id_col, "")
                old_desc[gid] = clean_title(r.get(args.desc_col, ""))
                gos = r.get(args.go_col, "") or ""
                old_go[gid] = {t for t in re.split(r"[,;| ]+", gos) if t.startswith("GO:")}

    new = {}
    with open(out / f"{prefix}_gene_annotations.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            new[r["gene_id"]] = r

    hier = None
    if args.go_obo:
        from goannotate.go import GOHierarchy
        hier = GOHierarchy(args.go_obo)

    def closure(terms):
        s = set()
        for t in terms:
            t2 = hier.normalize_id(t) if hier else t
            if not t2:
                continue
            s.add(t2)
            if hier:
                s |= hier.get_ancestors(t2)
        return s

    rows = []
    sym_classes = Counter()
    sym_outcome = Counter()
    by_class_outcome = Counter()
    prod = Counter()
    jaccards = []
    go_old_n = go_new_n = 0
    old_terms = []
    new_terms = []
    for gid, r in new.items():
        o = old_sym.get(gid, "")
        n = r["symbol"]
        cls = symbol_class(o, bnf)
        sym_classes[cls] += 1
        if cls == "informative":
            if n:
                if o.upper() == n.upper():
                    outcome = "match"
                elif related(o, n):
                    outcome = "related"
                else:
                    outcome = "different"
            else:
                outcome = "previous_only"
        else:
            outcome = "new_only" if n else "neither"
        sym_outcome[outcome] += 1
        by_class_outcome[(cls, outcome)] += 1

        od, nd = old_desc.get(gid, ""), r["product"]
        if od and nd:
            sim = name_similarity(normalize_gene_name(od), normalize_gene_name(nd))
            pst = "match" if sim >= 0.8 else ("partial" if sim >= 0.5 else "different")
        else:
            sim = ""
            pst = "previous_only" if od else ("new_only" if nd else "neither")
        prod[pst] += 1

        og = old_go.get(gid, set())
        ng = set(t for t in r["go_terms"].split(";") if t)
        if og:
            go_old_n += 1
            old_terms.append(len(og))
        if ng:
            go_new_n += 1
            new_terms.append(len(ng))
        jac = ""
        if og and ng and hier:
            a, b = closure(og), closure(ng)
            jac = len(a & b) / len(a | b) if (a | b) else 0.0
            jaccards.append(jac)
        rows.append([gid, o, cls, n, outcome, od, nd, f"{sim:.2f}" if sim != "" else "", pst,
                     str(len(og)), str(len(ng)), f"{jac:.2f}" if jac != "" else ""])

    total = len(new)
    print(f"Genes: {total:,}")
    print("\nPrevious symbols by class")
    for k, v in sym_classes.most_common():
        print(f"  {k:<30}{v:>8,}  ({100 * v / total:.1f}%)")
    print("\nSymbols, previous vs GOAnnotate")
    for k in ("match", "related", "different", "previous_only", "new_only", "neither"):
        print(f"  {k:<30}{sym_outcome[k]:>8,}")
    print("  new symbols on genes whose previous symbol was:")
    for cls in ("none", "gene_model_id", "genome_specific_identifier"):
        print(f"    {cls:<28}{by_class_outcome[(cls, 'new_only')]:>8,} of {sym_classes[cls]:,}")
    n_sym_new = sum(1 for r in new.values() if r["symbol"])
    print(f"  genes with an informative symbol: previous {sym_classes['informative']:,} "
          f"({100 * sym_classes['informative'] / total:.1f}%), new {n_sym_new:,} ({100 * n_sym_new / total:.1f}%)")
    if old_desc:
        print("\nProducts (previous description vs GOAnnotate product)")
        for k in ("match", "partial", "different", "previous_only", "new_only", "neither"):
            print(f"  {k:<30}{prod[k]:>8,}")
    if old_go:
        print("\nGO terms")
        print(f"  genes with GO: previous {go_old_n:,} ({100 * go_old_n / total:.1f}%), new {go_new_n:,} "
              f"({100 * go_new_n / total:.1f}%)")
        if old_terms and new_terms:
            print(f"  terms per annotated gene (median): previous {statistics.median(old_terms):.0f} (all terms), "
                  f"new {statistics.median(new_terms):.0f} (most specific terms only)")
        if jaccards:
            print(f"  genes with GO in both: {len(jaccards):,}; overlap of the two sets with ancestors added "
                  f"(Jaccard): median {statistics.median(jaccards):.2f}, "
                  f">= 0.5 for {100 * sum(j >= 0.5 for j in jaccards) / len(jaccards):.1f}%")

    path = out / "reports" / f"{prefix}_previous_comparison.tsv"
    path.parent.mkdir(exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(["gene_id", "previous_symbol", "previous_symbol_class", "symbol", "symbol_outcome",
                            "previous_product", "product", "product_similarity", "product_outcome",
                            "previous_n_go", "n_go", "go_jaccard_with_ancestors"]) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
