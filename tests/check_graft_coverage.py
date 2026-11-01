#!/usr/bin/env python3
"""
Coverage report for a GOAnnotate run on a GRAFT-updated GFF3.

Reports, for coding genes, the fraction with a product, a symbol and GO
terms, split by:
  - whether the input (EGAPx) annotation names the gene (gene= symbol)
  - graft_status
  - the source of the final symbol/product (symbol_status)
and the evidence behind GO terms (SwissProt vs TrEMBL top hit).

Usage: check_graft_coverage.py OUTPUT_DIR [PREFIX]
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path



def find_output(out, name):
    """A file of a GOAnnotate output directory: reports/, intermediate/ or the top level (older layout)."""
    for sub in ("reports", "intermediate", "."):
        p = out / sub / name
        if p.exists():
            return p
    return out / name

def pct(n, d):
    return f"{100 * n / d:5.1f}%" if d else "  n/a"


def main():
    out = Path(sys.argv[1])
    prefix = sys.argv[2] if len(sys.argv) > 2 else out.resolve().name
    ev = {r["gene_id"]: r for r in csv.DictReader(open(find_output(out, f"{prefix}_annotation_evidence.tsv")), delimiter="\t")}
    cmp = {r["gene_id"]: r for r in csv.DictReader(open(find_output(out, f"{prefix}_comparison.tsv")), delimiter="\t")}

    coding = [g for g, r in ev.items() if r["coding"] == "yes"]

    def row(label, genes):
        n = len(genes)
        prod = sum(1 for g in genes if ev[g]["protein_name"])
        sym = sum(1 for g in genes if ev[g]["gene_symbol"])
        go = sum(1 for g in genes if int(ev[g]["n_go_terms"]) > 0)
        print(f"  {label:<38}{n:>7,}  product {pct(prod, n)}  symbol {pct(sym, n)}  GO {pct(go, n)}")

    print("Coding genes")
    row("all", coding)
    named = [g for g in coding if cmp[g]["existing_symbol"]]
    unnamed = [g for g in coding if not cmp[g]["existing_symbol"]]
    row("EGAPx symbol present", named)
    row("EGAPx symbol absent (LOC, novel)", unnamed)
    print()

    print("Coding genes by graft_status")
    by_status = defaultdict(list)
    for g in coding:
        by_status[cmp[g].get("graft_status", "NA")].append(g)
    for st, genes in sorted(by_status.items(), key=lambda kv: -len(kv[1])):
        row(st, genes)
    print()

    print("Coding genes without an EGAPx symbol, by graft_status")
    by_status = defaultdict(list)
    for g in unnamed:
        by_status[cmp[g].get("graft_status", "NA")].append(g)
    for st, genes in sorted(by_status.items(), key=lambda kv: -len(kv[1])):
        row(st, genes)
    print()

    print("Symbol status (coding genes)")
    for st, n in Counter(ev[g]["symbol_status"] or "no hits" for g in coding).most_common():
        print(f"  {st:<38}{n:>7,}")
    print()

    print("Top hit source behind the annotation (coding genes with hits)")
    src = Counter(ev[g]["top_hit_source"] or "none" for g in coding if ev[g]["protein_name"])
    for s, n in src.most_common():
        print(f"  {s:<38}{n:>7,}")
    print()

    # Product names that look truncated or uninformative
    bad = re.compile(r"(^|\W)(protein|domain-containing|family)$|^uncharacterized|^putative$", re.I)
    flagged = [g for g in coding if ev[g]["protein_name"] and bad.search(ev[g]["protein_name"].replace("-like", ""))]
    print(f"Products that look truncated or generic: {len(flagged):,}")
    for g in flagged[:15]:
        print(f"  {g}\t{ev[g]['gene_symbol']}\t{ev[g]['protein_name']}\t{ev[g]['symbol_status']}")


if __name__ == "__main__":
    main()
