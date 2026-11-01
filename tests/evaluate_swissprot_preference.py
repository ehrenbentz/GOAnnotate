#!/usr/bin/env python3
"""
Evaluate, on an existing GOAnnotate run, what a nomenclature preference for
SwissProt hits of reference species would do to symbol agreement with the
existing annotation.

For every gene that has both an existing symbol and a GOAnnotate symbol,
the SwissProt hits from the preferred species (default human, mouse, rat)
whose cleaned protein name is similar to the GOAnnotate product (a proxy
for membership in the winning cluster) are collected. The symbol of the
best such hit, in species priority order and then by bitscore, is the
symbol the preference rule would assign. Agreement before and after is
reported.

Usage:
    evaluate_swissprot_preference.py -o OUTPUT_DIR [--prefix PREFIX] [--taxids 9606 10090 10116]
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.annotate import name_similarity, normalize_gene_name  # noqa: E402
from goannotate.filters import BadNameFilter, default_bad_names_path  # noqa: E402
from goannotate.models import BlastHit  # noqa: E402



def find_output(out, name):
    """A file of a GOAnnotate output directory: reports/, intermediate/ or the top level (older layout)."""
    for sub in ("reports", "intermediate", "."):
        p = out / sub / name
        if p.exists():
            return p
    return out / name

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", required=True, help="GOAnnotate output directory")
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--taxids", nargs="+", default=["9606", "10090", "10116"],
                    help="Preferred species in priority order (default: human mouse rat)")
    ap.add_argument("--similarity", type=float, default=0.5)
    ap.add_argument("--bad-names", default=str(default_bad_names_path()))
    args = ap.parse_args()

    out = Path(args.output)
    prefix = args.prefix or out.resolve().name
    priority = {t: i for i, t in enumerate(args.taxids)}
    bad = BadNameFilter(args.bad_names)

    # Gene universe and comparison
    tx2gene = {}
    with open(find_output(out, f"{prefix}_extraction_report.tsv")) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tx2gene[row["transcript_id"]] = row["gene_id"]

    genes = {}
    with open(find_output(out, f"{prefix}_comparison.tsv")) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["existing_symbol"] and row["goannotate_symbol"]:
                genes[row["gene_id"]] = row

    print(f"Genes with both an existing and a GOAnnotate symbol: {len(genes):,}")

    # Preferred SwissProt hits per gene
    candidates = defaultdict(list)
    n_lines = 0
    with open(find_output(out, f"{prefix}_diamond_results.tsv")) as fh:
        for line in fh:
            n_lines += 1
            if not line.startswith("sp|", line.find("\t") + 1):
                continue
            parts = line.rstrip("\n").split("\t")
            gene = tx2gene.get(parts[0])
            if gene is None or gene not in genes:
                continue
            hit = BlastHit(parts[0], parts[1], float(parts[2]), int(parts[3]), int(parts[4]),
                           int(parts[5]), int(parts[6]), int(parts[7]), int(parts[8]), int(parts[9]),
                           float(parts[10]), float(parts[11]), "\t".join(parts[12:]))
            taxid = hit.taxid
            if taxid not in priority:
                continue
            symbol = hit.gene_symbol
            if not symbol or bad.is_bad_symbol(symbol):
                continue
            name = bad.clean_name(hit.gene_name)
            if name is None or bad.is_bad_name(name):
                continue
            candidates[gene].append((priority[taxid], -hit.bitscore, symbol, name, taxid))

    before = Counter()
    after = Counter()
    transitions = Counter()
    with_candidate = 0
    examples = defaultdict(list)

    for gene, row in genes.items():
        existing = row["existing_symbol"].lower()
        current = row["goannotate_symbol"].lower()
        product_norm = normalize_gene_name(row["goannotate_product"])
        before["match" if existing == current else "mismatch"] += 1

        chosen = None
        for cand in sorted(candidates.get(gene, [])):
            _, _, symbol, name, taxid = cand
            if name_similarity(product_norm, normalize_gene_name(name)) >= args.similarity:
                chosen = symbol.lower()
                break
        if chosen is None:
            new = current
        else:
            with_candidate += 1
            new = chosen
        after["match" if existing == new else "mismatch"] += 1
        key = ("match" if existing == current else "mismatch") + "->" + ("match" if existing == new else "mismatch")
        transitions[key] += 1
        if len(examples[key]) < 8 and new != current:
            examples[key].append((gene, row["existing_symbol"], row["goannotate_symbol"], new))

    n = len(genes) or 1
    print(f"Genes with a preferred SwissProt hit in the winning cluster (proxy): {with_candidate:,} "
          f"({100 * with_candidate / n:.1f}%)")
    print()
    print(f"{'':<22}{'before':>10}{'after':>10}")
    for k in ("match", "mismatch"):
        print(f"{k:<22}{before[k]:>10,}{after[k]:>10,}")
    print(f"{'agreement':<22}{100 * before['match'] / n:>9.1f}%{100 * after['match'] / n:>9.1f}%")
    print()
    print("Transitions:")
    for k in ("match->match", "mismatch->match", "match->mismatch", "mismatch->mismatch"):
        print(f"  {k:<22}{transitions[k]:>10,}")
    for k in ("mismatch->match", "match->mismatch"):
        if examples[k]:
            print(f"\nExamples {k} (gene, existing, current, with preference):")
            for ex in examples[k]:
                print("  " + "\t".join(ex))


if __name__ == "__main__":
    main()
