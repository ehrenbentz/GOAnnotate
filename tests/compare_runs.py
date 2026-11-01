#!/usr/bin/env python3
"""
Compare two GOAnnotate runs of the same gene models (an earlier run and a
later one, e.g. before and after a change of options or databases).

Reads {prefix}_gene_annotations.tsv and reports/{prefix}_annotation_evidence.tsv
of both runs. Reports, for coding genes:
  - coverage of products, symbols and GO terms in each run
  - symbols: unchanged, changed, gained, lost (case-insensitive); the most
    frequent changes; whether lost or changed symbols were outside the
    reference nomenclature (reference proteome symbols of the new run)
  - products: identical, similar (name similarity >= 0.8), different, gained, lost
  - GO: terms per gene; overlap of the two sets with ancestors added (Jaccard)
  - confidence labels; symbol_status and coverage_status of the new run
Writes {new}/reports/{new_prefix}_vs_{old_prefix}.tsv with one row per gene.

Usage:
    compare_runs.py OLD_RUN_DIR NEW_RUN_DIR [--go-obo go.obo] [--old-prefix P] [--new-prefix P]
"""
import argparse
import csv
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.annotate import name_similarity  # noqa: E402

csv.field_size_limit(sys.maxsize)


def read_table(path, key="gene_id"):
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def find(run, prefix, name, sub=None):
    p = run / sub / f"{prefix}_{name}" if sub else run / f"{prefix}_{name}"
    return p if p.exists() else None


def norm(s):
    return " ".join(s.lower().split())


def reference_symbols(run, prefix):
    p = run / "intermediate" / f"{prefix}_reference_proteome.faa"
    syms = set()
    if p.exists():
        with open(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    m = re.search(r"GN=(\S+)", line)
                    if m:
                        syms.add(m.group(1).upper())
    return syms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old_run")
    ap.add_argument("new_run")
    ap.add_argument("--old-prefix")
    ap.add_argument("--new-prefix")
    ap.add_argument("--go-obo")
    ap.add_argument("--examples", type=int, default=15)
    args = ap.parse_args()

    old_dir, new_dir = Path(args.old_run), Path(args.new_run)
    op = args.old_prefix or old_dir.resolve().name
    np_ = args.new_prefix or new_dir.resolve().name
    old = read_table(old_dir / f"{op}_gene_annotations.tsv")
    new = read_table(new_dir / f"{np_}_gene_annotations.tsv")
    old_ev = read_table(find(old_dir, op, "annotation_evidence.tsv", "reports")) \
        if find(old_dir, op, "annotation_evidence.tsv", "reports") else {}
    new_ev = read_table(find(new_dir, np_, "annotation_evidence.tsv", "reports")) \
        if find(new_dir, np_, "annotation_evidence.tsv", "reports") else {}
    ref = reference_symbols(new_dir, np_)

    hier = None
    if args.go_obo:
        from goannotate.go import GOHierarchy
        hier = GOHierarchy(args.go_obo)

    def closure(terms):
        s = set()
        for t in terms:
            t2 = hier.normalize_id(t) if hier else t
            if t2:
                s.add(t2)
                if hier:
                    s |= hier.get_ancestors(t2)
        return s

    genes = [g for g in new if g in old and new[g].get("coding", "yes") != "no"]
    n = len(genes)
    sym = Counter()
    prod = Counter()
    sym_changes = Counter()
    sym_lost = Counter()
    lost_foreign = changed_foreign = 0
    jac = []
    old_terms, new_terms = [], []
    rows = []
    cov = {"old": Counter(), "new": Counter()}
    conf = {"old": Counter(), "new": Counter()}
    for g in genes:
        o, w = old[g], new[g]
        for label, r in (("old", o), ("new", w)):
            cov[label]["product"] += bool(r["product"])
            cov[label]["symbol"] += bool(r["symbol"])
            cov[label]["go"] += bool(r["go_terms"])
            conf[label][r["confidence"] or "none"] += 1
        os_, ns = o["symbol"].upper(), w["symbol"].upper()
        if os_ and ns:
            if os_ == ns:
                s_out = "same"
            else:
                s_out = "changed"
                sym_changes[f"{os_} -> {ns}"] += 1
                if ref and os_ not in ref:
                    changed_foreign += 1
        elif os_:
            s_out = "lost"
            sym_lost[os_] += 1
            if ref and os_ not in ref:
                lost_foreign += 1
        elif ns:
            s_out = "gained"
        else:
            s_out = "neither"
        sym[s_out] += 1

        op_, nw = o["product"], w["product"]
        if op_ and nw:
            if norm(op_) == norm(nw):
                p_out, sim = "identical", 1.0
            else:
                sim = name_similarity(norm(op_), norm(nw))
                p_out = "similar" if sim >= 0.8 else "different"
        else:
            sim = None
            p_out = "lost" if op_ else ("gained" if nw else "neither")
        prod[p_out] += 1

        og = {t for t in o["go_terms"].split(";") if t}
        ng = {t for t in w["go_terms"].split(";") if t}
        if og:
            old_terms.append(len(og))
        if ng:
            new_terms.append(len(ng))
        j = ""
        if og and ng:
            a, b = closure(og), closure(ng)
            j = len(a & b) / len(a | b) if a | b else 0.0
            jac.append(j)
        e = new_ev.get(g, {})
        rows.append([g, o["symbol"], w["symbol"], s_out, op_, nw, p_out,
                     "" if sim is None else f"{sim:.2f}", str(len(og)), str(len(ng)),
                     "" if j == "" else f"{j:.2f}", o["confidence"], w["confidence"],
                     e.get("symbol_status", ""), e.get("voted_symbol", ""), e.get("coverage_status", ""),
                     e.get("go_evidence_classes", ""), e.get("n_go_taxon_removed", ""), e.get("product_note", "")])

    def pct(x):
        return f"{100 * x / n:.1f}%" if n else "-"

    print(f"Coding genes in both runs: {n:,}")
    print(f"\n{'':<22}{op:>22}{np_:>22}")
    for k, lab in (("product", "with product"), ("symbol", "with symbol"), ("go", "with GO terms")):
        print(f"  {lab:<20}{cov['old'][k]:>14,} {pct(cov['old'][k]):>7}{cov['new'][k]:>14,} {pct(cov['new'][k]):>7}")
    for k in ("high", "medium", "low", "none"):
        print(f"  confidence {k:<9}{conf['old'][k]:>14,}{'':>8}{conf['new'][k]:>14,}")
    if old_terms and new_terms:
        print(f"  GO terms/gene (median, mean){'':<0}  {statistics.median(old_terms):.0f}, "
              f"{statistics.mean(old_terms):.1f}   ->   {statistics.median(new_terms):.0f}, {statistics.mean(new_terms):.1f}")

    print("\nSymbols (case-insensitive)")
    for k in ("same", "changed", "gained", "lost", "neither"):
        print(f"  {k:<10}{sym[k]:>8,}")
    if ref:
        print(f"  of the lost symbols, outside the reference nomenclature: {lost_foreign:,}")
        print(f"  of the changed symbols, old symbol outside the reference nomenclature: {changed_foreign:,}")
    print(f"  most frequent changes: " + "; ".join(f"{k} ({v})" for k, v in sym_changes.most_common(args.examples)))
    print(f"  most frequent lost: " + "; ".join(f"{k} ({v})" for k, v in sym_lost.most_common(args.examples)))

    print("\nProducts")
    for k in ("identical", "similar", "different", "gained", "lost", "neither"):
        print(f"  {k:<10}{prod[k]:>8,}")

    if jac:
        print(f"\nGO, genes with terms in both: {len(jac):,}; overlap with ancestors (Jaccard) median "
              f"{statistics.median(jac):.2f}, mean {statistics.mean(jac):.2f}, = 1 for "
              f"{100 * sum(1 for x in jac if x == 1.0) / len(jac):.1f}%")

    if new_ev:
        ev = [new_ev[g] for g in genes if g in new_ev]
        for col in ("symbol_status", "coverage_status", "product_note"):
            c = Counter(e.get(col, "") for e in ev if e.get(col, ""))
            if c:
                print(f"\n{np_} {col}: " + ", ".join(f"{k} {v:,}" for k, v in c.most_common()))
        classes = Counter()
        for e in ev:
            for item in e.get("go_evidence_classes", "").split(","):
                if ":" in item:
                    k, v = item.rsplit(":", 1)
                    classes[k] += int(v)
        if classes:
            tot = sum(classes.values())
            print(f"{np_} specific GO terms by strongest evidence: " +
                  ", ".join(f"{k} {v:,} ({100 * v / tot:.1f}%)" for k, v in classes.most_common()))
        removed = [int(e.get("n_go_taxon_removed") or 0) for e in ev]
        if any(removed):
            print(f"{np_} GO terms removed by taxon constraints: {sum(removed):,} on "
                  f"{sum(1 for x in removed if x):,} genes")

    out = new_dir / "reports" / f"{np_}_vs_{op}.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(["gene_id", "old_symbol", "new_symbol", "symbol_outcome", "old_product", "new_product",
                            "product_outcome", "product_similarity", "old_n_go", "new_n_go",
                            "go_jaccard_with_ancestors", "old_confidence", "new_confidence", "symbol_status",
                            "voted_symbol", "coverage_status", "go_evidence_classes", "n_go_taxon_removed",
                            "product_note"]) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
