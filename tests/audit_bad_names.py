#!/usr/bin/env python3
"""
Audit a pattern file against a Diamond results file without running the pipeline.

Writes the same report as the pipeline's {prefix}_name_audit.txt, for the
names and GN= symbols of a sample of hits (symbols are tested without
ortholog resolution here). Optionally reports how often discarded and kept
hits carry GO terms in the SQLite database.

Usage:
    audit_bad_names.py --blast-results results.tsv [--bad-names FILE] [--bad-names-extra FILE]
                       [--db GOAnnotate_db.sqlite] [--sample 500000] [--top 40] [-o report.txt]
"""
import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate.filters import BadNameFilter, NameAudit, default_bad_names_path, write_name_audit  # noqa: E402
from goannotate.models import BlastHit  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast-results", required=True)
    ap.add_argument("--bad-names", default=str(default_bad_names_path()))
    ap.add_argument("--bad-names-extra", action="append", default=[])
    ap.add_argument("--db", help="GOAnnotate SQLite database, to test GO coverage of discarded hits")
    ap.add_argument("--sample", type=int, default=500_000, help="Hits to sample (default 500,000)")
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("-o", "--output", default="/dev/stdout")
    args = ap.parse_args()

    bnf = BadNameFilter(args.bad_names, extra_files=args.bad_names_extra)
    audit = NameAudit()

    random.seed(1)
    n_total = sum(1 for _ in open(args.blast_results))
    keep_prob = min(1.0, args.sample / n_total)

    discarded_acc, kept_acc = [], []
    with open(args.blast_results) as fh:
        for line in fh:
            if random.random() > keep_prob:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            hit = BlastHit(p[0], p[1], float(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6]),
                           int(p[7]), int(p[8]), int(p[9]), float(p[10]), float(p[11]), "\t".join(p[12:]))
            name = bnf.process_name(hit.gene_name, audit)
            if name is None:
                discarded_acc.append(hit.accession)
                continue
            kept_acc.append(hit.accession)
            sym = hit.gene_symbol
            if sym:
                reason = bnf.explain_bad_symbol(sym)
                if reason:
                    audit.symbols_raw_bad += 1
                    audit.symbol_discards[reason] += 1
                    if len(audit.symbol_discard_examples[reason]) < 50:
                        audit.symbol_discard_examples[reason][sym] += 1
                else:
                    audit.kept_symbols[sym] += 1

    write_name_audit(audit, bnf, args.output, top=args.top)

    if args.db:
        import sqlite3
        with open(args.output, "a") as out:
            out.write("\nGO terms available for discarded vs kept hits (sample of 20,000 accessions each)\n")
            conn = sqlite3.connect(args.db)
            for label, accs in (("discarded", discarded_acc), ("kept", kept_acc)):
                sample = random.sample(accs, min(20_000, len(accs)))
                conn.execute("DROP TABLE IF EXISTS _acc")
                conn.execute("CREATE TEMP TABLE _acc (acc TEXT PRIMARY KEY)")
                conn.executemany("INSERT OR IGNORE INTO _acc VALUES (?)", [(a,) for a in sample])
                n_go = conn.execute("SELECT count(DISTINCT g.accession) FROM go_mapping g "
                                    "INNER JOIN _acc a ON g.accession = a.acc").fetchone()[0]
                out.write(f"  {label:<10} {n_go:,} of {len(sample):,} accessions have GO terms "
                          f"({100 * n_go / max(1, len(sample)):.1f}%)\n")
            conn.close()


if __name__ == "__main__":
    main()
