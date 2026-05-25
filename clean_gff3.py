#!/usr/bin/env python3
"""
clean_gff3_ids.py

Strips all existing ID/Parent/gene_id/transcript_id attributes from a GFF3 file
and replaces them with clean, sequentially numbered IDs. Parent references are
updated to match the new IDs. All other column 9 attributes are removed.

Reads entire input into memory, then performs two passes (ID mapping, then output).

Usage:
    gt gff3 -sort -tidy input.gff3 | python clean_gff3_ids.py > output.gff3
    python clean_gff3_ids.py input.gff3 output.gff3
    python clean_gff3_ids.py input.gff3 > output.gff3
    cat input.gff3 | python clean_gff3_ids.py > output.gff3
"""

import sys
import re


def main():
    if len(sys.argv) > 3:
        print(f"Usage: {sys.argv[0]} [input.gff3] [output.gff3]", file=sys.stderr)
        sys.exit(1)

    # Input: file arg or stdin
    if len(sys.argv) >= 2:
        with open(sys.argv[1]) as fh:
            lines = fh.readlines()
    else:
        lines = sys.stdin.readlines()

    # Output: file arg or stdout
    if len(sys.argv) == 3:
        fh_out = open(sys.argv[2], "w")
    else:
        fh_out = sys.stdout

    # Map old ID -> new ID
    id_map = {}
    # Per-feature-type counters
    type_counters = {}

    # First pass: build the ID mapping
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#") or line.strip() == "":
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue

        feature_type = fields[2]
        attrs = fields[8]

        id_match = re.search(r'ID=([^;]+)', attrs)
        if id_match:
            old_id = id_match.group(1)
            clean_type = feature_type.replace(" ", "_")
            type_counters[clean_type] = type_counters.get(clean_type, 0) + 1
            new_id = f"{clean_type}_{type_counters[clean_type]}"
            id_map[old_id] = new_id

    # Second pass: write output with new IDs and Parents
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#") or line.strip() == "":
            fh_out.write(line + "\n")
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            fh_out.write(line + "\n")
            continue

        attrs = fields[8]
        new_attrs = []

        id_match = re.search(r'ID=([^;]+)', attrs)
        if id_match:
            old_id = id_match.group(1)
            new_attrs.append(f"ID={id_map[old_id]}")

        parent_match = re.search(r'Parent=([^;]+)', attrs)
        if parent_match:
            old_parent = parent_match.group(1)
            new_parents = []
            for p in old_parent.split(","):
                if p in id_map:
                    new_parents.append(id_map[p])
                else:
                    print(f"WARNING: Parent '{p}' not found in ID map", file=sys.stderr)
                    new_parents.append(p)
            new_attrs.append(f"Parent={','.join(new_parents)}")

        fields[8] = ";".join(new_attrs)
        fh_out.write("\t".join(fields) + "\n")

    if fh_out is not sys.stdout:
        fh_out.close()

    print(f"Processed {sum(type_counters.values())} features", file=sys.stderr)
    for ftype, count in sorted(type_counters.items()):
        print(f"  {ftype}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
