#!/usr/bin/env python3
"""
Convert UniProt .dat (flat file) format to FASTA format.

Usage:
    ./dat_to_fasta.py input.dat.gz output.fasta
    ./dat_to_fasta.py input.dat output.fasta
"""

import sys
import gzip
import re

def convert_dat_to_fasta(input_file, output_file):
    """Convert UniProt .dat format to FASTA."""

    open_func = gzip.open if input_file.endswith('.gz') else open
    mode = 'rt' if input_file.endswith('.gz') else 'r'

    count = 0
    with open_func(input_file, mode) as inf, open(output_file, 'w') as outf:
        acc = ""
        entry_name = ""
        protein_name = ""
        org = ""
        ox = ""
        gn = ""
        sequence = ""
        in_seq = False
        is_reviewed = False  # Swiss-Prot = Reviewed, TrEMBL = Unreviewed

        for line in inf:
            if line.startswith('ID   '):
                entry_name = line.split()[1]
                is_reviewed = 'Reviewed' in line and 'Unreviewed' not in line
            elif line.startswith('AC   '):
                if not acc:  # Take first accession only
                    acc = line[5:].strip().split(';')[0]
            elif line.startswith('DE   RecName: Full='):
                protein_name = line.split('Full=')[1].split(';')[0].strip()
                # Remove ECO evidence codes {ECO:...}
                protein_name = re.sub(r'\s*\{ECO:[^}]+\}', '', protein_name)
            elif line.startswith('DE   SubName: Full=') and not protein_name:
                protein_name = line.split('Full=')[1].split(';')[0].strip()
                # Remove ECO evidence codes {ECO:...}
                protein_name = re.sub(r'\s*\{ECO:[^}]+\}', '', protein_name)
            elif line.startswith('OS   '):
                org = line[5:].strip().rstrip('.')
            elif line.startswith('OX   '):
                match = re.search(r'NCBI_TaxID=(\d+)', line)
                if match:
                    ox = match.group(1)
            elif line.startswith('GN   '):
                match = re.search(r'Name=([^;\s]+)', line)
                if match:
                    gn = match.group(1)
            elif line.startswith('SQ   '):
                in_seq = True
                sequence = ""
            elif line.startswith('//'):
                # End of entry - write FASTA
                if acc and sequence:
                    prefix = "sp" if is_reviewed else "tr"
                    header = f">{prefix}|{acc}|{entry_name} {protein_name}"
                    if org:
                        header += f" OS={org}"
                    if ox:
                        header += f" OX={ox}"
                    if gn:
                        header += f" GN={gn}"
                    outf.write(header + "\n")
                    # Write sequence in 60-char lines
                    for i in range(0, len(sequence), 60):
                        outf.write(sequence[i:i+60] + "\n")
                    count += 1
                    if count % 50000 == 0:
                        print(f"  Converted {count:,} entries...", file=sys.stderr)

                # Reset for next entry
                acc = ""
                entry_name = ""
                protein_name = ""
                org = ""
                ox = ""
                gn = ""
                sequence = ""
                in_seq = False
                is_reviewed = False
            elif in_seq:
                # Sequence line - remove spaces and add to sequence
                sequence += line.strip().replace(' ', '')

    print(f"Converted {count:,} entries to {output_file}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: ./dat_to_fasta.py input.dat[.gz] output.fasta")
        sys.exit(1)

    convert_dat_to_fasta(sys.argv[1], sys.argv[2])
