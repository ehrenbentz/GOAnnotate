#!/usr/bin/env python3
"""Tests of the GTF parser and the annotated GTF writer."""

import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from goannotate import gff3, gtf, output  # noqa: E402

failures = 0


def check(condition, label):
    global failures
    print(("ok    " if condition else "FAIL  ") + label)
    if not condition:
        failures += 1


def row(seqid, ftype, start, end, strand, frame, attrs):
    return "\t".join([seqid, "src", ftype, str(start), str(end), ".", strand, frame, attrs]) + "\n"


# NCBI style: gene and transcript lines, gene symbol under "gene", stop codon outside the CDS
NCBI = (
    "#gtf-version 2.2\n"
    + row("chr1", "gene", 101, 400, "+", ".", 'gene_id "G1"; transcript_id ""; db_xref "GeneID:11"; '
          'description "alpha; beta hydrolase"; gene "LOC11"; gene_biotype "protein_coding";')
    + row("chr1", "transcript", 101, 400, "+", ".", 'gene_id "G1"; transcript_id "T1"; db_xref "GeneID:11"; '
          'gene "LOC11"; product "uncharacterized protein LOC11"; transcript_biotype "mRNA";')
    + row("chr1", "exon", 101, 200, "+", ".", 'gene_id "G1"; transcript_id "T1"; exon_number "1";')
    + row("chr1", "exon", 301, 400, "+", ".", 'gene_id "G1"; transcript_id "T1"; exon_number "2";')
    + row("chr1", "CDS", 151, 200, "+", "0", 'gene_id "G1"; transcript_id "T1"; protein_id "XP_1"; '
          'product "uncharacterized protein LOC11"; transl_except "(pos:160..162,aa:Sec)"; '
          'transl_except "(pos:170..172,aa:Sec)";')
    + row("chr1", "CDS", 301, 347, "+", "1", 'gene_id "G1"; transcript_id "T1"; protein_id "XP_1";')
    + row("chr1", "start_codon", 151, 153, "+", "0", 'gene_id "G1"; transcript_id "T1";')
    + row("chr1", "stop_codon", 348, 350, "+", "0", 'gene_id "G1"; transcript_id "T1";')
    # minus strand gene, stop codon split from the CDS by an intron
    + row("chr1", "gene", 1001, 1500, "-", ".", 'gene_id "G2"; transcript_id ""; gene "ABC1";')
    + row("chr1", "transcript", 1001, 1500, "-", ".", 'gene_id "G2"; transcript_id "T2"; gene "ABC1"; '
          'product "ABC transporter 1";')
    + row("chr1", "exon", 1001, 1100, "-", ".", 'gene_id "G2"; transcript_id "T2";')
    + row("chr1", "exon", 1301, 1500, "-", ".", 'gene_id "G2"; transcript_id "T2";')
    + row("chr1", "CDS", 1301, 1450, "-", "2", 'gene_id "G2"; transcript_id "T2";')
    + row("chr1", "stop_codon", 1098, 1100, "-", "0", 'gene_id "G2"; transcript_id "T2";')
    + row("chr1", "gene", 2001, 2100, "+", ".", 'gene_id "G3"; transcript_id ""; gene_biotype "lncRNA";')
    + row("chr1", "transcript", 2001, 2100, "+", ".", 'gene_id "G3"; transcript_id "T3"; '
          'transcript_biotype "lnc_RNA";')
    + row("chr1", "exon", 2001, 2100, "+", ".", 'gene_id "G3"; transcript_id "T3";')
)

# Exon and CDS lines only, gene symbol under gene_name
MINIMAL = (
    row("chr2", "exon", 11, 100, "+", ".", 'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "Foo1"; '
        'gene_biotype "protein_coding";')
    + row("chr2", "CDS", 21, 80, "+", "0", 'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "Foo1";')
    + row("chr2", "exon", 11, 60, "+", ".", 'gene_id "ENSG1"; transcript_id "ENST2"; gene_name "Foo1";')
    + row("chr2", "exon", 501, 600, "-", ".", 'gene_id "ENSG2"; transcript_id "ENST3";')
    + row("chr2", "CDS", 511, 570, "-", "0", 'gene_id "ENSG2"; transcript_id "ENST3";')
)

# AUGUSTUS style: bare identifiers on gene and transcript lines, stop codon inside the CDS
AUGUSTUS = (
    row("chr3", "gene", 1, 300, "+", ".", "g1")
    + row("chr3", "transcript", 1, 300, "+", ".", "g1.t1")
    + row("chr3", "start_codon", 1, 3, "+", "0", 'transcript_id "g1.t1"; gene_id "g1";')
    + row("chr3", "CDS", 1, 300, "+", "0", 'transcript_id "g1.t1"; gene_id "g1";')
    + row("chr3", "exon", 1, 300, "+", ".", 'transcript_id "g1.t1"; gene_id "g1";')
    + row("chr3", "stop_codon", 298, 300, "+", "0", 'transcript_id "g1.t1"; gene_id "g1";')
)


def annotation(symbol, product, go=(), hits=True, aliases=()):
    return SimpleNamespace(consensus_symbol=symbol, consensus_name=product, product=product,
                           symbol_aliases=list(aliases), secondary_symbol="", confidence="high",
                           hits=[1] if hits else [], filtered_hits=[], specific_go_terms=set(go))


def test_attributes():
    pairs = gtf.split_gtf_attributes('gene_id "G1"; note "a; b"; exon_number 2; tag "x"; tag "y";')
    check([p[0] for p in pairs] == ["gene_id", "note", "exon_number", "tag", "tag"], "semicolon inside quotes kept")
    attrs = gtf.parse_gtf_attributes("exon", 'gene_id "G1"; note "a; b"; exon_number 2; tag "x"; tag "y";')
    check(attrs["note"] == "a; b" and attrs["exon_number"] == "2", "quoted and unquoted values")
    check(attrs["tag"] == "x,y", "repeated key joined with commas")
    check(gtf.parse_gtf_attributes("gene", "g1") == {"gene_id": "g1"}, "bare gene identifier")
    check(gtf.parse_gtf_attributes("transcript", "g1.t1") == {"transcript_id": "g1.t1"},
          "bare transcript identifier")


def test_join_stop_codons():
    check(gtf._join_stop_codons([(10, 30, 0), (30, 60, 1)], [(61, 63, 0)], "+") == [(10, 30, 0), (30, 63, 1)],
          "overlapping CDS segments (ribosomal slippage) stay separate")
    check(gtf._join_stop_codons([(13, 30, 2)], [(10, 12, 0)], "-") == [(10, 30, 2)], "minus strand stop joined")
    check(gtf._join_stop_codons([(10, 30, 0), (50, 61, 0)], [(62, 62, 0), (80, 81, 1)], "+")
          == [(10, 30, 0), (50, 62, 0), (80, 81, 1)], "stop codon split over two exons")


def test_parse_ncbi(tmp):
    path = tmp / "ncbi.gtf"
    path.write_text(NCBI)
    genes, stats = gtf.parse_gtf(str(path))
    gff3.assign_keys(genes)
    check([g.key for g in genes] == ["G1", "G2", "G3"], "genes in file order")
    g1, g2, g3 = genes
    t1 = g1.transcripts[0]
    check(t1.key == "T1" and t1.ftype == "mRNA", "transcript identifier and biotype")
    check(t1.exons == [(101, 200), (301, 400)], "exons")
    check(t1.cds == [(151, 200, 0), (301, 350, 1)], "adjacent stop codon joined to the CDS")
    check(t1.first_phase == 0, "first phase, plus strand")
    check(t1.cds_attrs.get("transl_except") == "(pos:160..162,aa:Sec),(pos:170..172,aa:Sec)",
          "repeated transl_except joined")
    check("gene_id" not in t1.cds_attrs and t1.cds_attrs.get("protein_id") == "XP_1", "CDS attributes")
    check(g1.attrs["description"] == "alpha; beta hydrolase" and g1.biotype == "protein_coding", "gene attributes")
    t2 = g2.transcripts[0]
    check(t2.cds == [(1098, 1100, 0), (1301, 1450, 2)], "stop codon in another exon kept as a segment")
    check(t2.first_phase == 2, "first phase, minus strand")
    check(not g3.coding and g3.biotype == "lncRNA" and g3.transcripts[0].ftype == "lnc_RNA", "non-coding gene")
    check(stats.n_genes == 3 and stats.n_coding_genes == 2 and stats.n_wrapper_genes == 0
          and stats.n_implicit_transcripts == 0, "statistics")
    universe = gff3.universe_from_genes(genes)
    check(universe.transcript_to_gene == {"T1": "G1", "T2": "G2", "T3": "G3"}, "gene universe")
    return genes


def test_parse_minimal(tmp):
    path = tmp / "minimal.gtf"
    path.write_text(MINIMAL)
    genes, stats = gtf.parse_gtf(str(path))
    gff3.assign_keys(genes)
    check([g.key for g in genes] == ["ENSG1", "ENSG2"], "genes built from exon and CDS lines")
    g1, g2 = genes
    check((g1.start, g1.end, g1.strand) == (11, 100, "+") and g1.wrapper, "gene span from its transcripts")
    check(g1.attrs.get("gene_name") == "Foo1" and g1.biotype == "protein_coding", "gene attributes from child lines")
    check(sorted(t.key for t in g1.transcripts) == ["ENST1", "ENST2"], "two transcripts")
    t1 = [t for t in g1.transcripts if t.key == "ENST1"][0]
    check(t1.cds == [(21, 80, 0)] and (t1.start, t1.end) == (11, 100) and t1.implicit, "transcript from child lines")
    check(g2.transcripts[0].strand == "-", "minus strand transcript")
    check(stats.n_wrapper_genes == 2 and stats.n_implicit_transcripts == 3, "statistics")
    from goannotate.sequences import gene_symbol_of
    check(gene_symbol_of(g1) == "Foo1" and gene_symbol_of(g2) == "", "gene_name read as the symbol")
    return genes


def test_parse_augustus(tmp):
    path = tmp / "augustus.gtf"
    path.write_text(AUGUSTUS)
    genes, stats = gtf.parse_gtf(str(path))
    gff3.assign_keys(genes)
    check(len(genes) == 1 and genes[0].key == "g1" and not genes[0].wrapper, "gene from a bare identifier")
    t = genes[0].transcripts[0]
    check(t.key == "g1.t1" and not t.implicit, "transcript from a bare identifier, gene from its children")
    check(t.cds == [(1, 300, 0)], "stop codon inside the CDS changes nothing")
    return genes


def test_write_ncbi(tmp, genes):
    annotations = {
        "G1": annotation("KIN1", "kinase 1", go=("GO:0004672", "GO:0005524"), aliases=("KINA",)),
        "G2": annotation("ABC1", "ABC transporter 1"),
    }
    output._missing_symbol_check = lambda value: value.startswith("LOC")
    output._missing_product_check = lambda value: "uncharacterized" in value
    out = tmp / "ncbi_out.gtf"
    output.write_annotated_gtf(str(tmp / "ncbi.gtf"), out, annotations, genes, mode="fill", hypothetical=False)
    text = out.read_text()
    lines = text.splitlines()
    original = NCBI.splitlines()
    check(len(lines) == len(original) and lines[0] == original[0], "line count and comment kept")
    gene1 = lines[1]
    check('goannotate_symbol "KIN1";' in gene1 and 'goannotate_product "kinase 1";' in gene1
          and 'goannotate_confidence "high";' in gene1 and 'goannotate_symbol_aliases "KINA";' in gene1,
          "gene results on the gene line")
    check('Ontology_term "GO:0004672"; Ontology_term "GO:0005524";' in gene1, "one Ontology_term per GO term")
    check('db_xref "GeneID:11";' in gene1, "existing db_xref kept")
    check('gene "KIN1";' in gene1 and 'gene "LOC11"' not in text.split("G2")[0], "identifier replaced on every line")
    check('description "kinase 1";' in gene1, "description replaced with the symbol")
    check('description "alpha; beta hydrolase"' not in gene1, "old description gone")
    tx1 = lines[2]
    check('product "kinase 1";' in tx1 and 'gene "KIN1";' in tx1 and "goannotate_symbol" not in tx1,
          "transcript line: product and symbol only")
    check('gene "KIN1";' in lines[3] and 'exon_number "1";' in lines[3], "symbol on exon lines")
    check('product "kinase 1";' in lines[5] and 'transl_except "(pos:160..162,aa:Sec)";' in lines[5],
          "CDS line: product replaced, other attributes kept")
    g2_lines = [l for l in lines if '"G2"' in l]
    g2_original = [l for l in original if '"G2"' in l]
    check(g2_lines[1:] == g2_original[1:], "lines of a gene that keeps its input annotation are unchanged")
    check('goannotate_symbol "ABC1";' in g2_lines[0] and g2_lines[0].startswith(g2_original[0].rstrip()[:-1]),
          "its gene line only gains the GOAnnotate attributes")
    g3_lines = [l for l in lines if '"G3"' in l]
    check(g3_lines == [l for l in original if '"G3"' in l], "gene without results unchanged")
    reparsed, _ = gtf.parse_gtf(str(out))
    check([g.id for g in reparsed] == ["G1", "G2", "G3"]
          and reparsed[0].transcripts[0].cds == genes[0].transcripts[0].cds, "output parses to the same models")

    out_add = tmp / "ncbi_add.gtf"
    output.write_annotated_gtf(str(tmp / "ncbi.gtf"), out_add, annotations, genes, mode="add", hypothetical=False)
    add_lines = out_add.read_text().splitlines()
    check('gene "LOC11";' in add_lines[1] and 'goannotate_symbol "KIN1";' in add_lines[1], "add mode keeps the symbol")
    check(add_lines[2:8] == original[2:8], "add mode changes only the gene line")


def test_write_minimal(tmp, genes):
    annotations = {"ENSG1": annotation("FOO1", "foo protein 1"),
                   "ENSG2": annotation("BAR2", "bar protein 2", go=("GO:0005515",))}
    out = tmp / "minimal_out.gtf"
    output.write_annotated_gtf(str(tmp / "minimal.gtf"), out, annotations, genes, mode="fill", hypothetical=False)
    lines = out.read_text().splitlines()
    check('gene_name "Foo1";' in lines[0] and 'goannotate_symbol "FOO1";' in lines[0],
          "input symbol kept, results on the lines of a transcript without a line")
    check('product "foo protein 1";' in lines[0] and 'product "foo protein 1";' in lines[1],
          "product on exon and CDS lines of a transcript without a line")
    check('gene_name "BAR2";' in lines[3] and 'gene_name "BAR2";' in lines[4], "symbol under gene_name")
    check('Ontology_term "GO:0005515";' in lines[3], "GO terms on the lines of a gene without a gene line")


def test_write_augustus(tmp, genes):
    class Module:
        def gene_attributes(self, key):
            return {"goannotate_signal_peptide": "yes", "Dbxref": "InterPro:IPR000001,Pfam:PF00001"}

        def transcript_attributes(self, key):
            return {"goannotate_sp_cleavage": "22-23"}

        def transcript_features(self, tx, gene):
            return ["\t".join([tx.seqid, "GOAnnotate", "sig_peptide", "1", "66", ".", "+", ".",
                               "ID=sig_peptide-g1.t1;Parent=g1.t1;Note=signal peptide%3B predicted"])]

    annotations = {"g1": annotation("SEC1", "secreted protein 1")}
    out = tmp / "augustus_out.gtf"
    output.write_annotated_gtf(str(tmp / "augustus.gtf"), out, annotations, genes, mode="fill",
                               hypothetical=False, modules=[Module()])
    lines = out.read_text().splitlines()
    check(lines[0].split("\t")[8].startswith('gene_id "g1"; goannotate_symbol "SEC1";'),
          "bare gene line written with attributes")
    check('db_xref "InterPro:IPR000001"; db_xref "Pfam:PF00001";' in lines[0]
          and 'goannotate_signal_peptide "yes";' in lines[0], "module gene attributes")
    check(lines[1].split("\t")[8].startswith('transcript_id "g1.t1";') and 'product "secreted protein 1";' in lines[1]
          and 'goannotate_sp_cleavage "22-23";' in lines[1] and 'gene_name "SEC1";' in lines[1],
          "bare transcript line written with attributes")
    check(lines[2].split("\t")[2] == "sig_peptide"
          and lines[2].split("\t")[8] == 'gene_id "g1"; transcript_id "g1.t1"; Note "signal peptide, predicted";',
          "module feature written as a GTF line")
    check(len(lines) == len(AUGUSTUS.splitlines()) + 1, "one line added")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_attributes()
        test_join_stop_codons()
        ncbi = test_parse_ncbi(tmp)
        minimal = test_parse_minimal(tmp)
        augustus = test_parse_augustus(tmp)
        test_write_ncbi(tmp, ncbi)
        test_write_minimal(tmp, minimal)
        test_write_augustus(tmp, augustus)
    print()
    if failures:
        print(f"{failures} check(s) failed")
        sys.exit(1)
    print("All checks passed")
