#!/usr/bin/env python3
"""
Compare GOAnnotate's extraction from a GFF3 + genome with a reference
protein FASTA (e.g. GRAFT's <name>_graft_proteins.fasta): IDs, header
fields and sequences.

Usage: compare_headers.py GFF3 GENOME REFERENCE_PROTEINS OUTDIR
"""
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import gff3, sequences  # noqa: E402


def main():
    gff, genome_path, ref, outdir = sys.argv[1:5]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    genes, _ = gff3.parse_gff3(gff)
    gff3.assign_keys(genes)
    genome = sequences.load_genome(genome_path)
    ext = sequences.extract_sequences(genes, genome, outdir, "hdrtest")
    del genome

    mine = {}
    for h, s in sequences.read_fasta(str(ext.protein_fasta)):
        i, f = sequences.parse_header(h)
        mine[i] = (h, f, s)
    theirs = {}
    for h, s in sequences.read_fasta(ref):
        i, f = sequences.parse_header(h)
        theirs[i] = (h, f, s.replace("*", ""))

    print(f"GOAnnotate proteins: {len(mine):,}; reference: {len(theirs):,}")
    print(f"IDs only in GOAnnotate: {len(set(mine) - set(theirs)):,}; only in reference: "
          f"{len(set(theirs) - set(mine)):,}")
    common = sorted(set(mine) & set(theirs))
    same_header = sum(1 for i in common if mine[i][0] == theirs[i][0])
    same_seq = sum(1 for i in common if mine[i][2].replace("*", "") == theirs[i][2])
    print(f"Identical header line: {same_header:,} of {len(common):,}")
    print(f"Identical sequence (internal stops removed): {same_seq:,} of {len(common):,}")
    field_diff = Counter()
    examples = {}
    for i in common:
        a, b = mine[i][1], theirs[i][1]
        for k in set(a) | set(b):
            if a.get(k) != b.get(k):
                field_diff[k] += 1
                examples.setdefault(k, []).append((i, a.get(k), b.get(k)))
    for k, n in field_diff.most_common():
        print(f"  field {k}: {n:,} differ, e.g.")
        for ex in examples[k][:4]:
            print(f"     {ex[0]}: GOAnnotate={ex[1]!r} reference={ex[2]!r}")
    for i in common[:3]:
        print("GOAnnotate:", mine[i][0])
        print("reference: ", theirs[i][0])


if __name__ == "__main__":
    main()
