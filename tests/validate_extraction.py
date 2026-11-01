#!/usr/bin/env python3
"""
Validate GOAnnotate's sequence extraction against an independent protein FASTA
(for example the protein set written by the annotation tool that produced the GFF3).

Usage:
    validate_extraction.py --gff3 FILE --genome FILE --reference proteins.fasta -o DIR [--mito-seqid ...]

Proteins are matched by transcript ID (first header token). Differences are
written to DIR/extraction_differences.tsv and summarized on stdout.
"""
import argparse
import logging
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import gff3, sequences  # noqa: E402


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--reference", required=True, help="Protein FASTA to compare against")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--prefix", default="validate")
    ap.add_argument("--mito-seqid", nargs="+", default=[])
    ap.add_argument("--codon-table", type=int, default=1)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    genes, stats = gff3.parse_gff3(args.gff3)
    gff3.assign_keys(genes)
    print(f"GFF3 parsed in {time.time() - t0:.1f}s")

    t0 = time.time()
    genome = sequences.load_genome(args.genome)
    print(f"Genome loaded in {time.time() - t0:.1f}s")

    t0 = time.time()
    result = sequences.extract_sequences(genes, genome, out, args.prefix,
                                         default_table=args.codon_table, mito_seqids=args.mito_seqid)
    print(f"Extraction in {time.time() - t0:.1f}s")

    ours = read_fasta(result.protein_fasta)
    ref = read_fasta(args.reference)

    only_ours = sorted(set(ours) - set(ref))
    only_ref = sorted(set(ref) - set(ours))
    shared = sorted(set(ours) & set(ref))

    categories = Counter()
    diff_path = out / "extraction_differences.tsv"
    with open(diff_path, "w") as fh:
        fh.write("transcript_id\tcategory\tours_len\tref_len\tdetail\n")
        for tid in shared:
            a = ours[tid]
            b = ref[tid].rstrip("*")
            if a == b:
                categories["identical"] += 1
                continue
            if a.replace("*", "X") == b.replace("*", "X"):
                categories["identical_except_stop_symbol"] += 1
                cat = "identical_except_stop_symbol"
            elif a.startswith(b) or b.startswith(a):
                cat = "prefix_difference"
                categories[cat] += 1
            elif a.endswith(b) or b.endswith(a):
                cat = "suffix_difference"
                categories[cat] += 1
            elif len(a) == len(b):
                n = sum(1 for x, y in zip(a, b) if x != y)
                cat = f"substitutions"
                categories[cat] += 1
                fh.write(f"{tid}\t{cat}\t{len(a)}\t{len(b)}\t{n} differing residues\n")
                continue
            else:
                cat = "different"
                categories[cat] += 1
            fh.write(f"{tid}\t{cat}\t{len(a)}\t{len(b)}\t{a[:30]}... vs {b[:30]}...\n")
        for tid in only_ours:
            fh.write(f"{tid}\tonly_in_goannotate\t{len(ours[tid])}\t0\t\n")
        for tid in only_ref:
            fh.write(f"{tid}\tonly_in_reference\t0\t{len(ref[tid])}\t\n")

    print()
    print(f"GOAnnotate proteins: {len(ours):,}; reference proteins: {len(ref):,}; shared IDs: {len(shared):,}")
    print(f"Only in GOAnnotate: {len(only_ours):,}; only in reference: {len(only_ref):,}")
    for cat, n in categories.most_common():
        print(f"  {cat:<32}{n:>8,}")
    print(f"Differences written to {diff_path}")


if __name__ == "__main__":
    main()
