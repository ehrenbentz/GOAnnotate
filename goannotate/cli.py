"""Command-line interface for GOAnnotate."""

import argparse
import logging
import os
import shutil
import sys

from goannotate import GOAnnotateError, __version__
from goannotate.filters import default_bad_names_path

logger = logging.getLogger("goannotate")


class VerticalHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Help formatter that keeps each argument's help on its own lines."""

    def __init__(self, prog, indent_increment=2, max_help_position=36, width=None):
        if width is None:
            width = max(shutil.get_terminal_size().columns, 100)
        super().__init__(prog, indent_increment, max_help_position, width)

    def _split_lines(self, text, width):
        return text.splitlines()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="goannotate",
        description=f"GOAnnotate {__version__}: functional annotation of gene models "
                    f"(gene symbols, protein names, GO terms) from Diamond searches against UniProt",
        formatter_class=VerticalHelpFormatter,
        epilog="""\
Input modes:
  Genome:        --gff3 FILE --genome FILE       extract CDS, translate, search proteins (blastp)
  Proteins:      --gff3 FILE --proteins FILE     search pre-extracted proteins (blastp); IDs checked
                                                 against the GFF3 transcripts before use
  Transcriptome: --transcripts FILE              search nucleotide sequences (blastx);
                                                 gene grouping from gene= header fields or --tx2gene
  Precomputed:   --transcripts FILE --blast-results FILE   annotate existing search results

Reference data: --db (SQLite from build_databases.py) or --go-mapping [--ncbi-idmapping --ncbi-geneinfo]
Diamond:        --diamond-swissprot and/or --diamond-trembl (omit when --blast-results is given)

Examples:
  goannotate --gff3 annotation.gff3 --genome genome.fasta \\
      --diamond-swissprot SwissProt.dmnd --diamond-trembl TrEMBL.dmnd \\
      --db GOAnnotate_db.sqlite --threads 32 -o vole_annotation

  goannotate --transcripts CDS.fasta --blast-results UniProt_results.tsv \\
      --db GOAnnotate_db.sqlite -o output_dir
"""
    )

    g = parser.add_argument_group("Input options")
    g.add_argument("--gff3", "--gff", dest="gff3", metavar="FILE",
                   help="GFF3 annotation. With --genome, sequences are extracted from it;\n"
                        "with --transcripts, it is used for the annotated GFF3 output only")
    g.add_argument("--genome", metavar="FILE", help="Genome FASTA matching the GFF3")
    g.add_argument("--proteins", metavar="FILE",
                   help="Pre-extracted protein FASTA for the GFF3 (used with --gff3, instead of\n"
                        "extraction from --genome). IDs must be the GFF3 transcript identifiers;\n"
                        "they are checked before use. With --genome as well, the proteins are\n"
                        "also compared with GOAnnotate's own extraction")
    g.add_argument("--transcripts", metavar="FILE",
                   help="Transcript (CDS or mRNA) FASTA; gene= header field groups isoforms")
    g.add_argument("--tx2gene", metavar="FILE",
                   help="Two-column file (transcript_id, gene_id) grouping transcripts into genes\n"
                        "for a FASTA without gene= header fields")
    g.add_argument("--blast-results", "--blast", dest="blast_results", metavar="FILE",
                   help="Precomputed Diamond/BLAST results, format 6 with stitle (skips the search)")
    g.add_argument("--gene-id-attr", default="ID", metavar="ATTR",
                   help="GFF3 gene attribute used as the gene identifier in outputs (default: ID)")
    g.add_argument("--transcript-id-attr", default="ID", metavar="ATTR",
                   help="GFF3 transcript attribute used as the transcript identifier (default: ID)")
    g.add_argument("--codon-table", type=int, default=1, metavar="N",
                   help="NCBI genetic code for translation (default: 1); transl_table attributes override it")
    g.add_argument("--mito-seqid", nargs="+", default=[], metavar="SEQID",
                   help="Sequence name(s) translated with the vertebrate mitochondrial code (table 2)")

    g = parser.add_argument_group("Diamond options")
    g.add_argument("--diamond-swissprot", metavar="FILE", help="SwissProt Diamond database (.dmnd)")
    g.add_argument("--diamond-trembl", metavar="FILE", help="TrEMBL Diamond database (.dmnd)")
    g.add_argument("--diamond-bin", metavar="PATH", help="Diamond executable (default: diamond on PATH)")
    g.add_argument("--sensitivity", default="fast",
                   choices=["fast", "mid-sensitive", "sensitive", "more-sensitive",
                            "very-sensitive", "ultra-sensitive"],
                   help="Diamond sensitivity mode (default: fast)")
    g.add_argument("--max-target-seqs", type=int, default=100, metavar="N",
                   help="Hits kept per query and per database (default: 100)")
    g.add_argument("--diamond-args", default="", metavar="STR",
                   help="Extra arguments passed to Diamond, quoted (for example \"-b 12 -c 1\")")
    g.add_argument("--force-search", action="store_true",
                   help="Rerun Diamond even when results exist in the output directory")

    g = parser.add_argument_group("Reference data options")
    g.add_argument("--db", metavar="FILE",
                   help="SQLite database built by build_databases.py (GO mapping and NCBI cross-references)")
    g.add_argument("--go-mapping", metavar="FILE",
                   help="Accession-to-GO mapping file (accession<TAB>GO:xxxx;GO:yyyy); alternative to --db")
    g.add_argument("--ncbi-idmapping", metavar="FILE",
                   help="UniProt idmapping_selected.tab for UniProt->NCBI GeneID lookup (with --go-mapping)")
    g.add_argument("--ncbi-geneinfo", metavar="FILE",
                   help="NCBI gene_info for GeneID->Symbol lookup (with --go-mapping)")
    g.add_argument("--go-obo", metavar="FILE",
                   help="GO OBO hierarchy file (auto-downloaded if missing or >30 days old)")
    g.add_argument("--namespace", nargs='+', default=["BP", "MF", "CC"], choices=["BP", "MF", "CC"],
                   help="GO namespaces to include (default: BP MF CC)")

    g = parser.add_argument_group("GO term options")
    g.add_argument("--go-evidence", default="all", choices=["all", "non-electronic", "experimental"],
                   help="GO annotations used, by evidence code (default: all)\n"
                        "  all:            every code except ND (no data)\n"
                        "  non-electronic: not IEA (electronic) annotations\n"
                        "  experimental:   EXP, IDA, IPI, IMP, IGI, IEP and the HTP codes only\n"
                        "Needs a database built by the current build_databases.py (evidence codes stored)")
    g.add_argument("--go-iea-min-hits", type=int, default=1, metavar="N",
                   help="A term whose only evidence is electronic (IEA) must be carried by at least N of\n"
                        "the hits GO terms are taken from (default: 1, no requirement)")
    g.add_argument("--taxid", type=int, metavar="TAXID",
                   help="NCBI taxon ID of the organism being annotated; enables the GO taxon constraints")
    g.add_argument("--taxon-constraints", default="on", choices=["on", "off"],
                   help="Remove GO terms that the Gene Ontology states do not apply to the --taxid organism\n"
                        "(default: on; needs --taxid and the taxonomy tables of the current build_databases.py,\n"
                        "or the two OBO files below)")
    g.add_argument("--taxon-constraints-obo", metavar="FILE",
                   help="go-computed-taxon-constraints.obo, instead of the constraints stored in --db")
    g.add_argument("--taxon-groupings-obo", metavar="FILE",
                   help="go-taxon-groupings.obo (GO taxon union classes), with --taxon-constraints-obo")
    g.add_argument("--taxdump", metavar="DIR",
                   help="Directory with NCBI nodes.dmp (and merged.dmp) for the --taxid lineage, instead of\n"
                        "the taxonomy tables stored in --db")

    g = parser.add_argument_group("Gene symbol options")
    g.add_argument("--reference-taxid", default="auto", metavar="TAXID",
                   help="Species whose nomenclature gene symbols follow; symbols of hits from other\n"
                        "species are replaced by the NCBI ortholog's symbol. auto (default) chooses from the\n"
                        "--taxid lineage: vertebrates human (9606), arthropods Drosophila melanogaster (7227),\n"
                        "fungi Saccharomyces cerevisiae (559292), otherwise human")
    g.add_argument("--ortholog-anchor-taxid", default="auto", metavar="TAXID",
                   help="Anchor species of the NCBI gene_orthologs rows used for resolution\n"
                        "(default auto: the reference species when it is an NCBI anchor, else 9606)")
    g.add_argument("--secondary-reference-taxid", default="auto", metavar="TAXID",
                   help="Second nomenclature, written to the secondary_symbol column: the symbol of the\n"
                        "gene's unambiguous reciprocal best hit in this species' proteome. auto (default):\n"
                        "human when the reference species is not human; none to disable")
    g.add_argument("--secondary-reference-proteome", metavar="FILE",
                   help="Protein FASTA of the secondary reference species with UniProt-style headers\n"
                        "(default: extracted from --diamond-swissprot)")
    g.add_argument("--no-resolve-orthologs", action="store_true",
                   help="Keep each hit's own GN= symbol; do not resolve through NCBI orthologs")
    g.add_argument("--symbol-case", default="upper", choices=["upper", "title", "lower", "source"],
                   help="Casing of output symbols: upper (AVP), title (Avp), lower (avp), or\n"
                        "source (as in the winning hit) (default: upper)")
    g.add_argument("--symbol-score", default="mean_sqrt", choices=["mean", "mean_sqrt"],
                   help="Score of a symbol group within the winning cluster: mean bitscore, or\n"
                        "mean bitscore times the square root of the group size (default: mean_sqrt)")
    g.add_argument("--reference-proteome", metavar="FILE",
                   help="Protein FASTA of the reference species with UniProt-style headers (GN= symbols)\n"
                        "for reciprocal best hits (default: extracted from --diamond-swissprot)")
    g.add_argument("--no-reciprocal", action="store_true",
                   help="Skip the checks of voted symbols against the reference proteome")
    g.add_argument("--symbol-support", type=float, default=0.9, metavar="F",
                   help="A voted symbol is kept when the gene's hit to that symbol's reference protein\n"
                        "scores at least F times the gene's best reference hit; otherwise the closest\n"
                        "reference protein's symbol is used (default: 0.9)")
    g.add_argument("--symbol-min-coverage", type=float, default=0.5, metavar="F",
                   help="The gene's alignment must also cover at least F of the reference protein\n"
                        "(default: 0.5)")
    g.add_argument("--symbol-min-query-coverage", type=float, default=0.8, metavar="F",
                   help="A voted symbol is also supported when the gene's alignment covers at least F of the\n"
                        "gene itself, so that partial gene models keep the symbol of the protein they are\n"
                        "part of (default: 0.8)")
    g.add_argument("--symbol-policy", default="ortholog", choices=["ortholog", "reference", "any"],
                   help="Voted symbols outside the reference nomenclature (zebrafish foxo3a, rodent Vmn1r,\n"
                        "outdated AARS); comparisons are case-insensitive\n"
                        "  ortholog:  replaced by a reference symbol only with evidence that it is the same\n"
                        "             gene: the symbol without a genome suffix (.2, a, b, l), a unique NCBI\n"
                        "             synonym (EIF2C3 -> AGO3), or a reciprocal best hit; the reference\n"
                        "             protein must be covered to --symbol-min-coverage. Others are kept\n"
                        "             (default)\n"
                        "  reference: replaced by the closest reference protein the gene covers, or removed\n"
                        "  any:       kept as voted")
    g.add_argument("--rbh-margin", type=float, default=0.9, metavar="F",
                   help="A reciprocal best hit counts as ortholog evidence only when the second-best\n"
                        "reference protein scores below F of the best (expanded families such as olfactory\n"
                        "receptors fail); 0 accepts every reciprocal best hit (default: 0.9)")
    g.add_argument("--drop-unmatched-symbols", action="store_true",
                   help="With --symbol-policy ortholog, remove non-reference symbols without ortholog\n"
                        "evidence instead of keeping them")
    g.add_argument("--keep-paralog-designators", action="store_true",
                   help="Keep paralog designators of other genomes in voted product names (\"O3-B\",\n"
                        "\"beta a\"); by default they are removed when the name without them is a reference\n"
                        "protein name or a name among the gene's hits")

    g = parser.add_argument_group("Filtering options")
    g.add_argument("--evalue", type=float, default=1e-5,
                   help="E-value threshold for the search and for hit filtering (default: %(default)s)")
    g.add_argument("--max-gene-hits", type=int, default=300,
                   help="Unique UniProt hits kept per gene after pooling its isoforms (each subject once,\n"
                        "best score); 0 for no cap (default: %(default)s)")
    g.add_argument("--top-n", type=int, default=100,
                   help="Informative hits kept per transcript after bad-name filtering (default: %(default)s)")
    g.add_argument("--hit-coverage", default="permissive", choices=["strict", "permissive", "off"],
                   help="Alignment coverage a hit needs to take part in naming and symbol voting\n"
                        "(default: permissive)\n"
                        "  strict:     covers --min-query-cov of the query AND --min-subject-cov of the subject\n"
                        "  permissive: covers either, OR the alignment is at least --min-align-length long\n"
                        "  off:        no requirement\n"
                        "Both modes also require --min-identity. Needs qlen/slen in the search results\n"
                        "(searches by earlier versions are rerun)")
    g.add_argument("--min-query-cov", type=float, default=0.5, metavar="F",
                   help="Fraction of the query the alignment must span (default: %(default)s)")
    g.add_argument("--min-subject-cov", type=float, default=0.5, metavar="F",
                   help="Fraction of the database protein the alignment must span (default: %(default)s)")
    g.add_argument("--min-identity", type=float, default=40.0, metavar="P",
                   help="Percent identity a hit needs under the coverage filter (default: %(default)s)")
    g.add_argument("--min-align-length", type=int, default=100, metavar="N",
                   help="Alignment length that passes the permissive coverage filter (default: %(default)s)")
    g.add_argument("--no-coverage-fallback", action="store_true",
                   help="Leave a transcript without hits when none passes the coverage filter; by default\n"
                        "its unfiltered hits are used and the gene is marked (coverage_status fallback,\n"
                        "confidence lowered by 0.1)")
    g.add_argument("--word-weighting", default="idf", choices=["idf", "stoplist", "off"],
                   help="Token weights in the name similarity (default: idf)\n"
                        "  idf:      stop words count 0; other words by log((N+1)/(df+1)) over the run's names\n"
                        "  stoplist: stop words count 0, other words 1\n"
                        "  off:      every word 1 (1.x behaviour)")
    g.add_argument("--stop-words", metavar="FILE",
                   help="Stop word list, one per line (default: the list shipped with GOAnnotate)")
    g.add_argument("--bad-names", metavar="FILE", default=None,
                   help="Name and symbol pattern file (default: the file shipped with GOAnnotate)")
    g.add_argument("--bad-names-extra", metavar="FILE", action="append", default=[],
                   help="Additional pattern file layered on top of --bad-names, e.g. locus-tag patterns\n"
                        "for your taxon; may be given more than once")
    g.add_argument("--consensus-threshold", type=float, default=0.5,
                   help="Similarity threshold for name clustering, 0-1 (default: %(default)s)")

    g = parser.add_argument_group("Comparison with existing annotation (GFF3 input)")
    g.add_argument("--no-compare", action="store_true",
                   help="Skip the comparison with functional annotation already in the GFF3")
    g.add_argument("--existing-symbol-attr", nargs="+", default=["gene"], metavar="ATTR",
                   help="Attribute(s) holding the existing gene symbol (default: gene)")
    g.add_argument("--existing-product-attr", nargs="+", default=["product", "description"], metavar="ATTR",
                   help="Attribute(s) holding the existing product name (default: product description)")
    g.add_argument("--strata-gene", nargs="*", default=["graft_status"], metavar="ATTR",
                   help="Gene-level attribute(s) to stratify the comparison by (default: graft_status)")
    g.add_argument("--strata-transcript", nargs="*", default=["graft_protein"], metavar="ATTR",
                   help="Transcript-level attribute(s) to stratify by (default: graft_protein)")

    g = parser.add_argument_group("Output options")
    g.add_argument("--output", "-o", metavar="DIR",
                   help="Output directory (required); its name is the file prefix unless --prefix is given")
    g.add_argument("--prefix", default=None, help="Output file prefix")
    g.add_argument("--annotation-mode", default="fill", choices=["add", "fill", "replace"],
                   help="How GOAnnotate results enter the annotated GFF3 (default: fill)\n"
                        "  add:     goannotate_* attributes, Ontology_term and Dbxref only\n"
                        "  fill:    also gene/Name/description/product where absent (existing\n"
                        "           annotation is kept; disagreements go to the review table)\n"
                        "  replace: also overwrite gene/Name/description/product")
    g.add_argument("--keep-existing-identifiers", action="store_true",
                   help="In fill and add modes, keep existing symbols that are identifiers (LOC numbers,\n"
                        "locus tags) and placeholder products (\"uncharacterized protein LOC...\"); by default\n"
                        "they count as missing and are filled with GOAnnotate's symbol and product")
    g.add_argument("--hypothetical", action="store_true",
                   help="Write product=hypothetical protein for coding models that end up without a\n"
                        "product (fill and replace modes). Off by default; table2asn assigns it itself")
    g.add_argument("--no-product-cleanup", action="store_true",
                   help="Do not apply the product_regex naming rules of the pattern file (EC numbers,\n"
                        "molecular weights, species names, \"<SYMBOL> protein\" word order) to products")

    g = parser.add_argument_group("Optional modules (external tools installed separately)")
    g.add_argument("--modules", default=None, metavar="LIST",
                   help="Modules to run, comma-separated (signalp), or auto to run every module whose\n"
                        "tool is found. Default: none. Modules add columns, GFF3 attributes and reports;\n"
                        "they do not change the symbol, product or GO assignment")
    g.add_argument("--check-tools", action="store_true",
                   help="Report which external tools (Diamond and the module tools) are available, then exit")
    g.add_argument("--force-modules", action="store_true",
                   help="Rerun the module tools even when their results exist in the output directory")
    g.add_argument("--signalp-bin", metavar="PATH",
                   help="signalp6 executable (default: signalp6 on PATH); may point into another\n"
                        "conda environment")
    g.add_argument("--signalp-mode", default="fast", choices=["fast", "slow", "slow-sequential"],
                   help="SignalP 6 model: fast (default), slow (full model, needs more than 14 GB RAM)\n"
                        "or slow-sequential (full model, memory of fast, about six times slower)")
    g.add_argument("--signalp-features", action="store_true",
                   help="Also write sig_peptide features (INSDC type) under coding transcripts with a\n"
                        "signal peptide in the annotated GFF3, with genomic coordinates")
    g.add_argument("--signalp-min-probability", type=float, default=0.0, metavar="P",
                   help="Signal peptide calls also need a SignalP probability of at least P\n"
                        "(default: 0, SignalP's own call)")
    g.add_argument("--signalp-processes", type=int, default=0, metavar="N",
                   help="Parallel signalp6 processes, one thread each (default: one per physical core,\n"
                        "limited by --threads and by available memory)")
    g.add_argument("--signalp-alt-starts", action="store_true",
                   help="Check proteins without a signal peptide for one at an internal methionine within\n"
                        "the first 60 residues that is followed by a hydrophobic core (a possible start\n"
                        "codon error of the gene model). Reported as a flag only; adds about a third to the\n"
                        "SignalP runtime")

    g = parser.add_argument_group("Performance options")
    g.add_argument("--threads", type=int, default=None,
                   help="Threads for Diamond and workers for annotation (default: all CPUs)")
    g.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING"],
                   help="Logging level (default: INFO)")

    parser.add_argument("--version", action="version", version=f"GOAnnotate {__version__}")
    return parser


def validate(args, parser) -> None:
    if args.modules is not None:
        from goannotate.modules import parse_selection
        try:
            parse_selection(args.modules)
        except ValueError as e:
            parser.error(f"--modules: {e}")
    if not 0.0 <= args.signalp_min_probability <= 1.0:
        parser.error("--signalp-min-probability must be between 0 and 1")
    if args.signalp_processes < 0:
        parser.error("--signalp-processes must be 0 (automatic) or a positive number")
    if args.check_tools:
        return
    if not args.output:
        parser.error("--output/-o is required")
    if not args.gff3 and not args.transcripts:
        parser.error("one of --gff3 (with --genome or --proteins) or --transcripts is required")
    if args.gff3 and not args.genome and not args.transcripts and not args.proteins:
        parser.error("--gff3 requires --genome (to extract sequences), --proteins or --transcripts")
    if args.genome and not args.gff3:
        parser.error("--genome requires --gff3")
    if args.proteins and not args.gff3:
        parser.error("--proteins requires --gff3 (use --transcripts for a FASTA without annotation)")
    if args.proteins and args.transcripts:
        parser.error("--proteins and --transcripts cannot be used together")
    if not args.db and not args.go_mapping:
        parser.error("--db or --go-mapping is required")
    if not args.blast_results and not (args.diamond_swissprot or args.diamond_trembl):
        parser.error("--diamond-swissprot and/or --diamond-trembl is required unless --blast-results is given")
    if bool(args.ncbi_idmapping) != bool(args.ncbi_geneinfo):
        parser.error("--ncbi-idmapping and --ncbi-geneinfo must be given together")
    if args.taxon_groupings_obo and not args.taxon_constraints_obo:
        parser.error("--taxon-groupings-obo is used with --taxon-constraints-obo")
    for name in ("min_query_cov", "min_subject_cov"):
        if not 0.0 <= getattr(args, name) <= 1.0:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1")
    if not 0.0 <= args.min_identity <= 100.0:
        parser.error("--min-identity is a percentage between 0 and 100")
    if args.go_iea_min_hits < 1:
        parser.error("--go-iea-min-hits must be at least 1")
    if args.bad_names is None:
        args.bad_names = str(default_bad_names_path())
    if args.threads is None:
        args.threads = os.cpu_count()


def check_tools(args) -> int:
    """Print the availability of Diamond and of every module's tool; exit 0 when all are found."""
    from goannotate.modules import detect_all
    from goannotate.search import find_diamond

    print(f"GOAnnotate {__version__}: external tools")
    all_found = True
    try:
        diamond = find_diamond(args.diamond_bin)
        print(f"  Diamond: {diamond}")
    except GOAnnotateError as e:
        print(f"  Diamond: not available ({e})")
        all_found = False
    for status in detect_all(args):
        print("  " + status.describe())
        all_found = all_found and status.available
    return 0 if all_found else 1


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format='%(asctime)s - %(levelname)s - %(message)s')
    validate(args, parser)

    if args.check_tools:
        return check_tools(args)

    logger.info(f"GOAnnotate {__version__}")
    try:
        from goannotate.pipeline import run
        run(args)
    except GOAnnotateError as e:
        logger.error(str(e))
        return 1
    logger.info("GOAnnotate completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
