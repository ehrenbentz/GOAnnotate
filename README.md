# GOAnnotate

Functional annotation of gene models using Diamond searches against UniProt.

GOAnnotate assigns protein names, gene symbols, and Gene Ontology (GO) terms to genes from Diamond hits against SwissProt and TrEMBL. Starting from a GFF3 and a genome, it extracts and translates the coding sequences itself, runs the Diamond searches, aggregates hits across isoforms, clusters them by protein name similarity, and uses a voting/scoring system to pick the best annotation for each gene. Optional modules add orthologs in close relatives (OrthoFinder, run before the vote so that well-supported orthologs decide it), signal peptides (SignalP 6), and protein domains with their GO terms (InterProScan 6). Every symbol is checked against the proteome of a reference species, and the reference gene's curated GO terms are added. Genomes that already carry functional annotation (for example from EGAPx) are compared gene by gene with the GOAnnotate result, and the existing annotation is kept in the output.

**Author:** E. J. Bentz (2025-2026)

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Database Setup](#database-setup)
- [Input Modes](#input-modes)
- [Running GOAnnotate](#running-goannotate)
- [How It Works](#how-it-works)
- [Output Files](#output-files)
- [Optional Modules](#optional-modules)
- [The bad_patterns.txt File](#the-bad_patternstxt-file)
- [Parameters](#parameters)
- [Utility Scripts](#utility-scripts)

## Requirements

- Python 3.9+
- [Diamond](https://github.com/bbuchfink/diamond) 2.x (not needed when precomputed results are supplied with `--blast-results`)
- Optional, for the modules (see Optional Modules):
  - SignalP 6 (`signalp` module)
  - InterProScan 6 with Nextflow 25.10.4 or later, Java 17 or later, a container runtime (Docker, Singularity, Apptainer or Podman) and the InterPro data (`interproscan` module)
  - [OrthoFinder](https://github.com/davidemms/OrthoFinder) and the NCBI [`datasets`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) command-line tool (`orthology` module)

No Python packages beyond the standard library are required.

## Installation

From a checkout, either run the wrapper script directly:

```bash
./GOAnnotate.py --help
```

or install the package, which provides the `goannotate` command:

```bash
pip install -e .
goannotate --help
```

## Quick Start

```bash
# 1. Build databases (one-time setup, several hours)
./build_databases.py -o /path/to/databases --threads 8

# 2. Annotate a genome from its GFF3
goannotate \
    --gff3 Annotations.gff3 \
    --genome Genome.fasta \
    --diamond-swissprot /path/to/databases/SwissProt_Diamond_DB.dmnd \
    --diamond-trembl /path/to/databases/TrEMBL_Diamond_DB.dmnd \
    --db /path/to/databases/GOAnnotate_db.sqlite \
    --threads 32 \
    --output my_annotations
```

GOAnnotate extracts the proteins, searches them against SwissProt and TrEMBL with `diamond blastp`, and writes the annotation. Diamond results are kept in the output directory and reused by later runs with the same output directory (`--force-search` reruns them), so the annotation parameters can be tuned without repeating the search.

## Database Setup

`build_databases.py` downloads, verifies and builds the reference data GOAnnotate needs from GOA, UniProt and NCBI. The databases can be reused across projects and should be rebuilt periodically, all sources together, so that they come from one release.

```bash
./build_databases.py -o /path/to/databases --protein-dbs metazoa --threads 64 --yes
```

`--protein-dbs` takes the lineage of the protein databases as an NCBI taxon name or taxon ID (`metazoa`, `viridiplantae`, `fungi`, `7742`), or `all`. Several lineages are given comma-separated and built from one download, one pair of Diamond databases each (`--protein-dbs metazoa,viridiplantae`); `--skip-sqlite` adds a lineage beside an existing SQLite database. Each lineage FASTA records the UniProt release it was built from, and a later run that downloads another release rebuilds it with its Diamond database.

**SQLite database** (`GOAnnotate_db.sqlite`, always built; pass it with `--db`):

| Table | Content | Source |
|---|---|---|
| `go_mapping` | UniProt accession -> GO terms, each with the evidence codes behind it (`GO:0005515:IDA,IPI;GO:0005634:IEA`; NOT annotations excluded) | GOA `goa_uniprot_all.gaf.gz` |
| `ncbi_idmapping` | UniProt accession -> NCBI GeneID | UniProt `idmapping_selected.tab.gz` |
| `ncbi_geneinfo` | GeneID -> official symbol and NCBI gene type (protein-coding, pseudo, ...) | NCBI `gene_info.gz` |
| `ncbi_gene_symbols` | (taxon, name) -> GeneID, the name being the official symbol, the locus tag or a synonym (column `kind`) | NCBI `gene_info.gz` |
| `ncbi_orthologs` | GeneID -> ortholog GeneID | NCBI `gene_orthologs.gz` |
| `taxonomy_nodes`, `taxonomy_merged` | taxon -> parent, rank, name; merged taxon IDs (for the lineage of the annotated organism) | NCBI `taxdump.tar.gz` |
| `go_taxon_constraints`, `go_taxon_unions` | GO term -> `in_taxon` / `never_in_taxon` statements, and the member taxa of GO's taxon union classes | GO `go-computed-taxon-constraints.obo`, `go-taxon-groupings.obo` |
| `sources` | URL, size, MD5, verification, release and dates of every downloaded file | |
| `metadata` | build date, script version, UniProt release, row counts, protein databases | |

Table groups can be built or replaced on their own: `--tables go,taxonomy` builds only those groups (`go`, `idmapping`, `gene_info`, `orthologs`, `taxonomy`), and with `--update` they replace the same tables in an existing database (each table in one transaction; the rest of the database is untouched). A database built before the evidence codes were stored still works; `--go-evidence` then has nothing to act on, and the GO taxon constraints need the `taxonomy` tables or the two OBO files.

**Protein databases** (with `--protein-dbs LINEAGE`): UniProtKB Swiss-Prot and TrEMBL restricted to a lineage, given as an NCBI taxon name (`metazoa`) or taxon ID (`33208`), or `all`. The script writes `<lineage>_SwissProt.dmnd` and `<lineage>_TrEMBL.dmnd` (with taxonomy, so searches can be narrowed with `--taxonlist`), their FASTA files and taxon maps, and with `--combined-db` also `<lineage>_UniProtKB.dmnd`. Pass the two databases to GOAnnotate with `--diamond-swissprot` and `--diamond-trembl`.

**Verification.** Each file is downloaded to `<name>.part` (resumed after an interruption) and renamed only when complete. It is checked against the checksum the provider publishes: MD5 and size from UniProt's `RELEASE.metalink`, the `.md5` file beside each GOA and NCBI taxonomy file. NCBI `gene/DATA` publishes no checksums, so `gene_info.gz` and `gene_orthologs.gz` are checked against the server's size and by decompressing them completely. A verified file receives `<name>.verified.json`; a file without one is verified before use, and a file from an earlier release is downloaded again. The SQLite database is written as `GOAnnotate_db.sqlite.tmp` and renamed only when every table is complete; an existing database is never overwritten. `build_info.tsv` and the `sources` table record what was used.

Options:

```
-o, --output-dir      Directory for all files (required)
--protein-dbs LINEAGE Also build lineage-restricted Swiss-Prot and TrEMBL Diamond databases
--combined-db         With --protein-dbs, also build the combined Swiss-Prot + TrEMBL database
--threads N           Threads for Diamond (default: all CPUs)
--keep-downloads      Keep the downloaded archives after a successful build
--flat-files          Also write Uniprot_All_GO_mapping.tsv, idmapping_selected.tab and gene_info.tsv
                      (for --go-mapping, --ncbi-idmapping, --ncbi-geneinfo without --db)
--skip-sqlite         Build the protein databases only
--tables GROUPS       Table groups to build: go, idmapping, gene_info, orthologs, taxonomy, or all
--update              Replace the --tables groups in an existing database
--yes                 Do not ask for confirmation
```

Disk: about 80 GB of downloads (GAF 29 GB, TrEMBL FASTA 41 GB, ID mapping 7 GB, gene_info 1.6 GB), the SQLite database (about 30 GB) and, for Metazoa, about 45 GB of FASTA and Diamond files. The GAF annotations are sorted by accession on disk (GNU sort, about 60 GB of temporary files) and written one accession at a time.

## Input Modes

### Genome mode: `--gff3` and `--genome`

The GFF3 is parsed into genes, transcripts, exons and CDS. Every transcript type is accepted (mRNA, transcript, lnc_RNA, ncRNA, miRNA, tRNA, V_gene_segment, ...), as are pseudogenes, transcripts without a gene parent, and genes with CDS children but no transcript. The CDS of each coding transcript is assembled in translation order, trimmed by the GFF3 phase of its first segment (5' partial models), and translated. The genetic code is table 1 unless a `transl_table` attribute says otherwise or the sequence is named with `--mito-seqid` (table 2); `transl_except` attributes (selenocysteine, for example) are applied. Three FASTA files are written:

| File | Content |
|------|---------|
| `intermediate/{prefix}_AA.faa` | Translated CDS, one record per coding transcript; the Diamond query set |
| `intermediate/{prefix}_CDS.fna` | CDS nucleotide sequence per coding transcript |
| `intermediate/{prefix}_transcripts.fna` | Spliced exon sequence per transcript, coding or not |

Headers follow the format of GRAFT's output files, with fields that have no value left out:

```
>rna-MOCH1_014632-R1 gene=gene-MOCH1_014632 symbol=MYADM product="myeloid-associated differentiation marker isoform 1-T1" type=mRNA locus=chr1:70314112-70315068(+)
>rna-MOCH1_012492-R1 gene=gene-MOCH1_012492 product="olfactory receptor 4S2" type=mRNA locus=chr1:23246-24202(-)
```

`symbol` is the gene's `gene` attribute (or its `Name` when that is not just the ID or locus tag); `product` is the CDS `product`, else the transcript's; `type` is the transcript feature type (`mRNA` for any coding transcript that is not `mRNA`, `transcript` or a `*_gene_segment`); `locus` is 1-based and inclusive, the CDS span in the protein and CDS files and the exon span in the transcript file. On the GRAFT graft_09 output, all 43,357 protein headers are identical to GRAFT's. `reports/{prefix}_extraction_report.tsv` lists every transcript with its coding status, lengths, codon table and flags (`phase_trimmed`, `length_not_multiple_of_3`, `start_codon:XXX`, `no_stop_codon`, `internal_stop:n`, `partial`, `pseudo`, `transl_except`, `exception`). Nothing is dropped silently.

The gene universe is every gene in the GFF3, coding or not. Non-coding genes appear in all mapping files with empty annotations, so the files can serve directly as the background for enrichment analyses of gene-level counts. By default gene and transcript identifiers are the GFF3 `ID` values; `--gene-id-attr` and `--transcript-id-attr` select another attribute (`Name`, `locus_tag`, ...), which must be unique.

Proteins are searched with `diamond blastp`.

### GTF input: `--gtf`

`--gtf annotation.gtf` takes the place of `--gff3` in every input mode and selects the GTF parser. Features are grouped by their `gene_id` and `transcript_id` attributes, which are the gene and transcript identifiers of the outputs. `gene` and `transcript` lines are used when the file has them; a gene or a transcript without a line of its own is built from the span of its exon and CDS lines, with the `gene*` and `transcript*` attributes of its first line. Column 9 holding a bare identifier on gene and transcript lines (AUGUSTUS) is accepted. `stop_codon` features are joined to the CDS, so that the CDS includes the stop codon as in GFF3 whether or not the GTF excludes it; a GTF with CDS features and no `stop_codon` features is used as given, and a warning is logged. The frame column is read as the phase. The transcript type is `transcript_biotype` or `transcript_type` when present, and the gene type `gene_biotype` or `gene_type`. The existing gene symbol is read from `gene` (NCBI) or `gene_name` (Ensembl, GENCODE), products from `product` and `description`. A key that occurs several times on a line (`transl_except`, `db_xref`, `tag`) is read as one comma-separated value. The annotated output is a GTF (see Output Files).

### Pre-extracted proteins: `--gff3` and `--proteins`

A protein FASTA made from the same GFF3 (for example GRAFT's `<name>_graft_proteins.fasta`) is searched instead of an extraction. It is checked against the GFF3 before use:

- the run stops if an ID occurs twice, if an ID is not a transcript identifier in the GFF3 (after `--transcript-id-attr`), or if a `gene=` header field names a different gene than the transcript's parent in the GFF3;
- coding transcripts without a protein, proteins for transcripts without CDS, and proteins whose length differs from the CDS length by more than one residue are reported (proteins written with internal stops removed are shorter by their number).

With `--genome` as well, the transcript and CDS files are still extracted and every supplied protein is compared with GOAnnotate's own translation. Problems are listed in `reports/{prefix}_protein_input_check.tsv` and counted in the summary.

### Transcriptome mode: `--transcripts`

A nucleotide FASTA (CDS or mRNA) is searched with `diamond blastx`, and each transcript with a hit is translated in the frame of its best hit: the protein is the open reading frame that contains the alignment, extended upstream to the first methionine after the preceding in-frame stop (or to the transcript start) and downstream to the first in-frame stop (or to the transcript end). The proteins are then searched with `blastp`, and the rest of the run (reciprocal best hits, reference names and GO terms, the modules) proceeds as for a genome. Transcripts without a blastx hit get no protein and no annotation. `reports/{prefix}_orf_translation.tsv` lists each protein with its strand, frame, ORF coordinates, completeness (complete, 5_partial, 3_partial or partial) and the hit that defined it; transcripts whose top hits lie in more than one frame or in separate regions (possible chimeras or frameshifts) are flagged there and translated from the best hit. A translated transcript whose protein finds no blastp hit (short proteins with hits near the e-value threshold) is annotated from its blastx hits. Fragments, proteins lacking one or both ends, keep a voted reference symbol when the alignment to that symbol's reference protein covers `--symbol-min-fragment-coverage` (0.6) of the fragment, instead of the 0.8 required of a gene model: a fragment's protein runs to the transcript ends, which rarely align. The transcript table gives each protein's `completeness` (complete, 5_partial, 3_partial or partial) and the evidence for each end (`completeness_evidence`): the 3' end is present when the reading frame ends in a stop codon; the 5' end when it begins with a methionine and either an in-frame stop precedes it in the transcript or the alignment to the reference protein (that of the gene's symbol, else the best reference hit) begins within 10 residues of the reference start. The gene table gives the most complete of the gene's transcripts, and the FASTA headers carry `completeness=`. A fragment's confidence score is lowered by 0.05 for each missing end. `--orf-translation off` annotates from the blastx hits alone, as earlier versions did, without the reference checks or the modules. A protein FASTA is detected automatically and searched with `blastp`.

Isoforms are grouped into genes through a `gene=` field in the headers, through `--tx2gene FILE` (two columns: transcript_id, gene_id), or, without either, through Trinity identifiers (`DN1000_c0_g1_i1` and `DN1000_c0_g1_i2` are isoforms of `DN1000_c0_g1`; the `TRINITY_` prefix, Trinity 2.0's `TR1000|c0_g1_i1` and the older `comp1000_c0_seq1` are read as well). `--isoform-grouping off` keeps every sequence as its own gene.

### Precomputed results: `--transcripts` and `--blast-results`

Skips the search and annotates existing Diamond or BLAST results (format 6 with the 13 fields `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle`). Adding `--gff3` produces the annotated GFF3 output as well. For nucleotide transcripts with the translation on, the file is taken as the blastx results that define the proteins (it then needs the `qlen` and `slen` columns), and the proteins are searched with the Diamond databases.

## Running GOAnnotate

### Genome mode

```bash
goannotate \
    --gff3 Annotations.gff3 \
    --genome Genome.fasta \
    --diamond-swissprot SwissProt.dmnd \
    --diamond-trembl TrEMBL.dmnd \
    --db GOAnnotate_db.sqlite \
    --mito-seqid chrM \
    --threads 48 \
    --diamond-args "-b 8 -c 1" \
    --output vole_annotation
```

### Transcriptome mode

```bash
goannotate \
    --transcripts Trinity.fasta \
    --diamond-swissprot SwissProt.dmnd \
    --diamond-trembl TrEMBL.dmnd \
    --db GOAnnotate_db.sqlite \
    --taxid 35019 \
    --modules signalp,orthology,interproscan \
    --output trinity_annotation
```

### Precomputed results

```bash
goannotate \
    --transcripts CDS.fasta \
    --blast-results UniProt_results.tsv \
    --db GOAnnotate_db.sqlite \
    --gff3 Annotations.gff3 \
    --output my_annotations
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--gff3` or `--gtf` + `--genome`, `--gff3` or `--gtf` + `--proteins`, or `--transcripts` | Input (see Input Modes) |
| `--db` or `--go-mapping` | SQLite database or flat-file GO mapping |
| `--diamond-swissprot` and/or `--diamond-trembl` | Diamond databases (not needed with `--blast-results`) |
| `--output` | Output directory |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--blast-results` | none | Precomputed search results; skips Diamond |
| `--tx2gene` | none | Transcript-to-gene file for FASTA input without `gene=` fields |
| `--gene-id-attr`, `--transcript-id-attr` | `ID` | GFF3 attributes used as identifiers |
| `--codon-table` | 1 | NCBI genetic code; `transl_table` attributes override it |
| `--mito-seqid` | none | Sequence name(s) translated with table 2 |
| `--sensitivity` | fast | Diamond sensitivity mode |
| `--max-target-seqs` | 100 | Hits per query and per database |
| `--diamond-args` | none | Extra Diamond arguments, quoted |
| `--diamond-bin` | `diamond` on PATH | Diamond executable |
| `--force-search` | off | Rerun Diamond even when results exist |
| `--go-obo` | auto-download | GO hierarchy file (downloaded if missing or stale) |
| `--ncbi-idmapping`, `--ncbi-geneinfo` | none | NCBI cross-references (not needed with `--db`) |
| `--evalue` | 1e-5 | E-value threshold (search and filtering) |
| `--top-n` | 100 | Top hits to keep per transcript after filtering |
| `--hit-coverage` | permissive | Coverage a hit needs for naming and voting: `strict` (query and subject), `permissive` (either, or a long alignment), `off` |
| `--min-query-cov`, `--min-subject-cov` | 0.5, 0.5 | Fractions of the query and of the database protein the alignment must span |
| `--min-identity` | 40 | Percent identity a hit needs under the coverage filter |
| `--min-align-length` | 100 | Alignment length that passes the permissive filter |
| `--no-coverage-fallback` | off | Leave transcripts without hits when none passes the coverage filter (default: use their unfiltered hits and mark the gene) |
| `--word-weighting` | idf | Token weights in the name similarity: `idf`, `stoplist`, `off` |
| `--stop-words` | bundled file | Stop word list for the name similarity |
| `--bad-patterns` | bundled file | Name, symbol and product pattern file |
| `--bad-patterns-extra` | none | Additional pattern file layered on the bundled one (repeatable) |
| `--write-bad-patterns [FILE]` | `bad_patterns.txt` | Write a copy of the bundled pattern file and exit; an existing file is not overwritten |
| `--consensus-threshold` | 0.5 | Similarity threshold for name clustering |
| `--namespace` | BP MF CC | GO namespaces to include |
| `--go-evidence` | all | GO annotations used by evidence code: `all` (except ND), `non-electronic` (no IEA), `experimental` |
| `--go-iea-min-hits` | 1 | Hits that must carry a term whose only evidence is IEA |
| `--taxid` | none | NCBI taxon ID of the annotated organism; enables the GO taxon constraints |
| `--taxon-constraints` | on | Remove GO terms the ontology states do not apply to the `--taxid` organism |
| `--taxon-constraints-obo`, `--taxon-groupings-obo` | from `--db` | The two GO constraint files, instead of the database tables |
| `--taxdump` | from `--db` | Directory with `nodes.dmp` for the `--taxid` lineage, instead of the database tables |
| `--reference-taxid` | auto | Species whose nomenclature symbols follow. auto chooses from the `--taxid` lineage: vertebrates human (9606), arthropods Drosophila melanogaster (7227), fungi Saccharomyces cerevisiae (559292), otherwise human |
| `--ortholog-anchor-taxid` | auto | Anchor species of the NCBI gene_orthologs rows (the reference species when it is an NCBI anchor, else 9606) |
| `--secondary-reference-taxid` | auto | Second nomenclature for the `secondary_symbol` column (human when the reference species is not human; `none` disables it) |
| `--secondary-reference-proteome` | from SwissProt DB | Protein FASTA of the secondary reference species |
| `--no-resolve-orthologs` | off | Keep each hit's own GN= symbol |
| `--symbol-case` | upper | Output casing of symbols: upper, title, lower, source |
| `--reference-proteome` | from SwissProt DB | Reference-species protein FASTA for reciprocal best hits |
| `--no-reciprocal` | off | Skip the checks of voted symbols against the reference proteome |
| `--symbol-support` | 0.9 | Minimum score of the voted symbol's reference hit relative to the gene's best reference hit |
| `--symbol-min-coverage` | 0.5 | Minimum fraction of the reference protein covered by the gene's alignment |
| `--symbol-policy` | ortholog | Voted symbols outside the reference nomenclature: `ortholog` replaces them by a reference symbol only with evidence that it is the same gene (unique NCBI synonym, unambiguous reciprocal best hit) and keeps the rest; `reference` uses the closest covered reference protein or removes them; `any` keeps them |
| `--rbh-margin` | 0.9 | A reciprocal best hit is ortholog evidence only when the second-best reference protein scores below this fraction of the best; 0 accepts all |
| `--drop-unmatched-symbols` | off | With `ortholog`, remove non-reference symbols without ortholog evidence |
| `--symbol-min-query-coverage` | 0.8 | A voted reference symbol is also supported when the alignment covers this fraction of the gene (partial gene models) |
| `--symbol-min-fragment-coverage` | 0.6 | The same for translated transcripts whose protein lacks one or both ends |
| `--keep-paralog-designators` | off | Keep other genomes' paralog designators (`O3-B`, `beta a`) in voted product names |
| `--annotation-mode` | fill | How results enter the annotated GFF3 (see Output Files) |
| `--keep-existing-identifiers` | off | In fill and add modes, keep existing LOC-type symbols and placeholder products instead of filling them |
| `--hypothetical` | off | Write `product=hypothetical protein` for unnamed coding models |
| `--no-product-cleanup` | off | Do not apply the `product_regex:` naming rules to products |
| `--no-compare` | off | Skip the comparison with existing annotation |
| `--existing-symbol-attr` | gene | Attribute(s) holding an existing gene symbol |
| `--existing-product-attr` | product description | Attribute(s) holding an existing product name |
| `--strata-gene`, `--strata-transcript` | graft_status, graft_protein | Attributes the comparison is stratified by |
| `--modules` | none | Optional modules to run (`signalp`, `orthology`, `interproscan`; `auto` runs every module whose tool is found); see Optional Modules |
| `--check-tools` | off | Report which external tools are available, then exit |
| `--force-modules` | off | Rerun the module tools even when their results exist |
| `--signalp-bin`, `--signalp-mode`, `--signalp-features` | `signalp6` on PATH, fast, off | SignalP module settings |
| `--signalp-min-probability` | 0 (SignalP's call) | Minimum SignalP probability for a signal peptide call |
| `--signalp-processes` | one per physical core | Parallel `signalp6` processes, one thread each |
| `--signalp-alt-starts` | off | Flag proteins with a signal peptide from an internal start codon |
| `--interproscan-workflow`, `--interproscan-revision` | none | InterProScan 6 workflow: a local clone at a release tag, or the GitHub name with a release |
| `--interproscan-datadir`, `--interproscan-interpro` | none, the single release found | InterPro data directory and release |
| `--interproscan-nextflow`, `--interproscan-profile` | `nextflow` on PATH, auto | Nextflow executable and container runtime |
| `--interproscan-apps` | all except deep-learning analyses | InterProScan analyses to run |
| `--interproscan-local` | off | Run every analysis locally instead of looking up precalculated matches at EBI |
| `--interproscan-workers` | half of `--threads` | CPUs for InterProScan's parallel tasks |
| `--interproscan-keep-work` | off | Keep Nextflow's work directory |
| `--interpro-go` | add | InterPro GO terms: `add` to every gene, `fill` genes without GO, `off` |
| `--interpro-products` | fill | Name coding genes left without a product from their InterPro domains |
| `--nxf-home` | `~/.nextflow` | Directory for Nextflow's own files (NXF_HOME) |
| `--relatives` | 3 | Close relatives chosen automatically for the orthology module |
| `--relative` | none | Relative to use instead (NCBI assembly accession or taxon ID; repeatable) |
| `--relative-min-named` | 0.2 | Lowest informative fraction of a relative's protein-coding genes (symbol, or a link to a named reference gene) |
| `--relative-named-ratio` | 0.75 | A relative must reach this fraction of the best candidate at its level or closer |
| `--relative-named-margin` | 0.10 | Per lineage level, how much better named a farther relative must be to rank above a closer one |
| `--same-species` | on | Use the RefSeq annotation of the query species itself, when one exists, as one of the `--relatives` |
| `--relatives-dir` | `intermediate/{prefix}_orthology/relatives` | Where relatives' annotations are downloaded |
| `--orthofinder-bin`, `--orthofinder-threads`, `--datasets-bin` | on PATH, `--threads`, on PATH | OrthoFinder and NCBI datasets settings |
| `--ortholog-symbol` | weighted | How the relatives' ortholog symbols are used: `weighted` (close relatives decide, distant ones vote and fill), `fill` (fill only), `off` |
| `--ortholog-close-rank` | family | A relative is close when its common ancestor with the query species is at or below this rank |
| `--reference-go` | add | GO terms of the reference gene named by a gene's supported symbol: `add`, `fill` (genes without GO from hits), `off` |
| `--threads` | all CPUs | Threads for Diamond, workers for annotation and threads for the module tools |
| `--prefix` | output dir name | Prefix for output filenames |

## How It Works

### Workflow overview

The stages run in this order (modules only when selected with `--modules`):

1. **Input** (Phase 0): gene universe and query proteins from the GFF3 and genome, or from FASTA files.
2. **Diamond search** (Phase 1): query proteins against SwissProt and TrEMBL, or precomputed results with `--blast-results`.
3. **Orthology module, first part**: relatives are chosen from the query's lineage, their RefSeq annotations downloaded, and OrthoFinder run on the query and relatives. Each gene's one-to-one ortholog symbol, and whether close relatives support it, is passed to the vote.
4. **Per-gene vote** (Phases 2 to 8): hit cleaning, isoform aggregation, clustering, cluster and symbol selection, GO terms from the hits, symbol propagation. A symbol supported by close relatives decides the vote; a distant relative's symbol counts as one more vote.
5. **Reference checks** (Phase 8b): reciprocal best hits and sequence support against the reference proteome, the secondary nomenclature, symbol aliases.
6. **Orthology module, second part**: the close relatives' symbols are applied (`confirmed`, `taken`, `conflict_kept`, `filled`, `kept`), with orthology evidence added to the confidence.
7. **Reference-gene GO**: GO terms of the reference gene named by each supported symbol (Phase 7, `--reference-go`).
8. **SignalP and InterProScan modules**: signal peptides; InterPro entries, InterPro GO terms, and domain-based products for genes left without one.
9. **Concordance and confidence** of the final annotation (Phase 9).
10. **Comparison** with an existing annotation in the GFF3 (Phase 10), and the output files.

### Phase 0: Input Parsing

In genome mode, the GFF3 defines the gene universe and the sequences are extracted and translated as described under Input Modes. In the FASTA modes, GOAnnotate reads the headers: the transcript ID is the first whitespace-delimited field and a `gene=` tag groups transcripts into genes. Every gene appears in the output files, whether or not it has hits.

### Phase 1: Diamond Search

Proteins (or transcripts) are searched against SwissProt and TrEMBL in two separate runs and the results are concatenated. Running the databases separately matters because `--max-target-seqs` limits hits per database: in a combined database a query would often receive 100 TrEMBL hits and no SwissProt hits, losing the curated SwissProt annotations that carry the best GO terms. Hits are parsed and filtered by e-value, keyed by transcript ID. Each hit's subject ID is checked for `sp|` or `tr|` prefixes to identify SwissProt vs. TrEMBL origin. The protein name and gene symbol (from the `GN=` field in the stitle) are extracted from each hit.

### Phase 2: Hit Cleaning and Top-N Selection

Each hit's protein name is cleaned and tested with the patterns in `bad_patterns.txt`. The purpose is to stop designations that belong to another genome or annotation pipeline from being carried into this one, and to keep names that carry no information from being used as products:

1. **Qualifiers** (`qualifier_regex:`) are stripped: database prefixes, `PREDICTED:`, `LOW QUALITY PROTEIN:`, `(Fragment)`, isoform and variant designations, `putative`/`probable`, and a trailing `-like` (the NCBI pipeline's paralog marker for one genome).
2. **Identifiers** (`identifier_regex:`) are removed as tokens: LOC numbers, RIKEN clone IDs, C#orf# names in short and long form, KIAA/FLJ/MGC/DKFZ clone IDs, UPF and DUF numbers, Gm numbers, clone accessions. "Bombesin receptor-activated protein C6orf89" becomes "Bombesin receptor-activated protein".
3. The remaining name is **discarded** when it is a placeholder (`placeholder:`, e.g. "uncharacterized protein", "protein") or matches an `exact:`, `contains:` or `regex:` pattern. A hit whose name is discarded contributes nothing.

Gene symbols (`GN=`) are carried as found (a FlyBase transcript or protein suffix, `kay-RA`, `CG1647-PB`, which cDNA-clone entries carry in the GN= field, is dropped) until ortholog resolution has been attempted, so that a LOC number or species-numbered symbol (`Olr1468`, `Zfp1006`) can still reach the human symbol of its gene. Only symbols that fail to resolve and match a `symbol_regex:` pattern are blanked.

Family-level names (`generic_regex:` patterns, plus a structural rule for one- or two-word names ending in a protein-class word without a member designator) are kept, but within the winning symbol group the best-scoring specific name is used as the product: "Olfactory receptor 4K36" in preference to "Olfactory receptor".

**Hit coverage filter.** A hit then has to cover enough of both sequences to stand for the whole protein rather than one shared domain (after PANNZER2's hit filter). Under `--hit-coverage permissive` (default) the alignment must span at least 50% of the query or 50% of the database protein (`--min-query-cov`, `--min-subject-cov`), or be at least 100 residues long (`--min-align-length`), at 40% identity or more (`--min-identity`), so fragments and multidomain proteins pass; `strict` requires both coverages; `off` disables the test. The search records query and subject lengths for this (`qlen`, `slen`; results written by earlier versions lack them and are searched again). A transcript none of whose hits pass keeps its unfiltered hits by default, and the gene is marked `fallback` in the evidence table with its confidence lowered by 0.1; `--no-coverage-fallback` leaves such transcripts without hits. The best query and subject coverage reached within the winning symbol group are reported for every gene.

Top-N selection happens *after* bad name filtering and the coverage filter. This is important: if filtering happened after truncation, uninformative hits would consume slots and push out real annotations. With the current order, you get up to N informative hits per transcript.

### Phase 3: Isoform Aggregation

Hits from all transcripts of the same gene are merged into a single pool, each tagged with its source transcript.

### Phase 4: Agglomerative Clustering

The merged hits for each gene are clustered by protein name similarity using average-linkage agglomerative clustering. This groups hits that describe the same (or very similar) proteins together, even when the exact wording differs between species or databases.

Similarity between two hits is the larger of the Jaccard token similarity and the `SequenceMatcher` ratio of the normalized protein names, with a modification for gene symbols. The Jaccard similarity is weighted by word (`--word-weighting idf`, default): stop words (`goannotate/data/stop_words.txt`: protein, putative, homolog, family, fragment, isoform and the like) count 0, and every other word counts log((N + 1) / (df + 1)), where N is the number of distinct names in the run and df the number of names containing the word, so a designator that appears in few names ("4K36", "FOXO3") outweighs a class word that appears in many ("receptor", "kinase"). `stoplist` drops the stop words only; `off` restores the unweighted similarity.

- If both hits have valid gene symbols and the symbols match, similarity gets a +0.2 bonus (capped at 1.0).
- If both hits have valid gene symbols and the symbols are dissimilar (< 0.8 similarity), the overall similarity is **capped at 0.3**, regardless of how similar the protein names are.
- If either hit lacks a gene symbol, name similarity is used directly.

The symbol-aware cap is what prevents paralogs from merging. For example, coagulation factors F7 and F9 have very similar protein names ("Coagulation factor VII" vs. "Coagulation factor IX") but different gene symbols. Without the cap, they would cluster together and one would dominate. With it, they stay in separate clusters.

Clustering proceeds by iteratively merging the two most similar clusters until no pair exceeds the similarity threshold (default 0.5).

### Phase 5: Winning Cluster Selection

Each cluster is scored and the highest-scoring cluster wins. There is no fallback: the winning cluster always determines the annotation. The scoring formula is:

```
score = sum(bitscores) * sqrt(cluster_size) * (1 + 0.1 * (n_isoforms - 1))
```

This rewards clusters with strong evidence (high bitscores), broad support (many hits), and corroboration across isoforms. The square root on cluster size prevents a large but low-scoring cluster from beating a small cluster with very strong hits.

### Phase 6: Name and Symbol Selection

Before voting, the symbol of every hit is resolved to the nomenclature of one reference species (`--reference-taxid`, default human). A hit from another species is mapped to its NCBI GeneID, either through the UniProt ID mapping or through its (species, GN= symbol) pair in NCBI gene_info, then to the reference-species ortholog in NCBI gene_orthologs, and the ortholog's symbol replaces the hit's own. Hits from the reference species keep their symbol, and hits that cannot be resolved keep their GN= symbol. This removes species-specific variants (Zfp428, Olfr577, Ick) in favor of the current reference symbols (ZNF428, OR51G2, CILK1) and makes symbols comparable across genes. The two tables this needs (`ncbi_orthologs`, `ncbi_gene_symbols`) are built by `build_databases.py`; an existing database is extended with `build_databases.py -o DIR --add-orthologs`. `--no-resolve-orthologs` disables the step.

Within the winning cluster, hits are grouped by gene symbol. The symbol group with the highest score is selected, and both the gene symbol and protein name come from that group. The score is the mean bitscore of the group times the square root of the group size, so that agreement among many species counts once symbols are in one nomenclature.

A **single-hit guard** prevents annotation errors in large clusters: if a symbol group has only one hit and there are other groups with two or more hits, the single-hit group cannot win. This stops a lone misannotated hit from overriding the majority. An additional guard rejects a single-hit symbol that represents less than 10% of a cluster with 10+ hits.

If no valid gene symbol exists in the winning cluster, GOAnnotate falls back to NCBI cross-referencing: it looks up the best hit's UniProt accession in the NCBI ID mapping to find a GeneID, then looks up that GeneID in NCBI's gene_info to find a symbol. This step requires either the `--db` SQLite database or the `--ncbi-idmapping` and `--ncbi-geneinfo` flat files.

The protein name comes from the highest-bitscore hit in the winning symbol group. Output symbols are written in one casing convention (`--symbol-case`, default upper).

### Phase 7: GO Term Collection

GO terms are collected in two steps:

1. **Primary collection:** The top hits (by bitscore) from the winning symbol group, 10 by default (`--go-source-hits`), are looked up in the GO mapping. GO annotations belong to proteins: each hit contributes the annotations of its UniProt entry, and the symbol group makes sure these are proteins of the same gene.

2. **SwissProt supplement:** The pipeline then scans all remaining hits in the winning cluster for SwissProt entries that share the winning gene symbol or protein name. GO terms from matching SwissProt hits are added to the pool. This is important because TrEMBL hits typically dominate the bitscore rankings (especially for non-model organisms with closely related species in TrEMBL), but SwissProt entries carry richer, manually curated GO annotations.

**Evidence.** Databases built by the current `build_databases.py` store the GO evidence codes behind every accession-term pair. `--go-evidence` chooses which annotations are used: `all` (default; every code except ND, whose terms are the ontology roots), `non-electronic` (no IEA annotations), or `experimental` (EXP, IDA, IPI, IMP, IGI, IEP and the high-throughput codes). **Support.** Once every source of GO terms has been added (the hits, the reference gene, InterProScan), a term that comes from the hits alone must be supported by `--go-min-hits` (default 5) of the source hits, or by a majority of the gene's source hits with GO annotations when that is fewer (counting at most 10, and never fewer than 2: a gene with two such hits needs both, one with four needs three), and a term whose only evidence is electronic (IEA) by at least `--go-iea-min-hits` (default 2), never by fewer hits than other terms of the gene. A gene with a single source hit keeps that hit's non-electronic terms when the alignment covers at least 70% of both proteins (`--go-single-hit-coverage`) at 50% identity or more (`--go-single-hit-identity`); such a hit stands for the same protein, the inference a curator records as ISS. A hit supports a term when its protein is annotated with the term or with a more specific one. Exempt are terms with an experimental annotation (EXP, IDA, IPI, IMP, IGI, IEP and the high-throughput codes), terms that a SwissProt entry of the gene itself (GN= symbol equal to the gene's symbol) carries with non-electronic evidence, and terms confirmed by the reference gene or by InterPro. A removed term is replaced by its most specific ancestors that pass the requirement, so the supported part of the claim stays ("kinase activity" where only one hit says "MAP kinase kinase activity"). Terms from the reference gene are kept, except one that none of the gene's hits carries, whose only evidence for the reference gene is electronic (IEA), and that no close relative's one-to-one ortholog carries either (with the orthology module, each close relative whose ortholog's UniProt entries carry a term, or a more specific one, also counts as one supporting source for terms from the hits). Removed terms are listed in `reports/{prefix}_go_evidence.tsv` with status `removed_support` and their source (`hits` or `reference`); `--go-min-hits 1 --go-iea-min-hits 1` keeps every term. The codes behind each assigned term, its strongest evidence class and the number of supporting hits are written to `reports/{prefix}_go_evidence.tsv`; the evidence table summarizes them per gene (`go_evidence_classes`).

**From the reference gene.** Once the symbols are final (after Phase 8b and the orthology module), a gene whose symbol is supported or unverified (`rbh`, `rbh_override`, `rbh_added`, `vote_supported`, `vote_unverified`, `reference_best_hit`, `foreign_ortholog`, `ortholog`, `ortholog_fill`) also receives the GO terms of the reference gene that symbol names: the terms of the reference proteome's proteins carrying the symbol, under the same `--go-evidence` policy, taxon constraints and namespaces. A reference gene whose only name is an identifier (FlyBase `BcDNA:SD01241`, human `C1orf122`) gives no symbol; a gene left without a symbol that is its reciprocal best hit still receives its terms. With `--reference-go add` (default) they join the gene's terms without redundancy (a term already present, or an ancestor of one, adds nothing); `fill` gives them only to genes without GO terms from hits; `off` leaves GO to the hits and InterProScan. For a human reference these are the curated GOA annotations of the human gene; `reports/{prefix}_go_evidence.tsv` marks their source (`reference`) and lists the proteins behind each term.

**Taxon constraints.** The Gene Ontology states for many terms the taxa they apply to (`in_taxon`, `never_in_taxon`; constraints pass to subtypes and parts). With `--taxid` (the NCBI taxon of the organism being annotated) GOAnnotate computes the organism's lineage from the taxonomy tables of the database (or `--taxdump`), expands GO's taxon union classes, and removes every GO term that a constraint on it or on one of its ancestors excludes for that lineage: photosynthesis terms from a mammal, lactation from a fish. The constraints come from the database tables or from `--taxon-constraints-obo` and `--taxon-groupings-obo`; removed terms are listed in `reports/{prefix}_go_evidence.tsv` with the constraint that removed them. `--taxon-constraints off` disables the step; without `--taxid` it is skipped.

After collection, GO terms are filtered by namespace (BP, MF, CC) and then reduced to the most specific terms by removing redundant ancestor terms using the GO hierarchy from the OBO file. For example, if a gene has both GO:0005634 (nucleus) and GO:0005737 (cytoplasm), both are kept because neither is an ancestor of the other. But if it has both GO:0005634 (nucleus) and GO:0005622 (intracellular anatomical structure, a parent of nucleus), the parent is removed.

### Phase 8: Symbol Propagation

After all genes are annotated independently, a post-annotation pass looks for genes that share the same normalized protein name. If some of those genes have a gene symbol and others do not, the most common symbol that also occurs among a symbol-less sibling's own hits is propagated to it. This fills in gaps where the vote for a particular gene left it without a symbol although its hits carry one. Identical product names alone are not used: TrEMBL names such as "Disks large homolog 5-like" are shared by unrelated genes that resemble one domain of the named protein.

### Phase 8b: Symbol Checks Against the Reference Proteome

The vote decides identity from many database hits, but its naming can be wrong in two ways: a gene that shares a domain with a well-studied protein can be named after it, and for paralog pairs the better-studied sibling tends to win the vote. A second pass checks every voted symbol against the complete reference-species proteome (the SwissProt entries of the reference taxon, which for human cover essentially every protein-coding gene; extracted from the SwissProt Diamond database with `diamond getseq`, or given with `--reference-proteome`). The reference proteome is searched against the query proteins and the query proteins against it; both searches take seconds. A gene is the reciprocal best hit (RBH) of reference symbol X when X's protein has that gene as its best query hit and X is among the gene's hits.

- A voted symbol equal to the RBH symbol is confirmed (`rbh`), and the reference protein's name becomes the product.
- Product of a confirmed symbol (`rbh`, `rbh_added`, `vote_supported`): the reference protein's name replaces the voted name, which often comes from an automatic family name of TrEMBL ("RNA helicase" for eIF4A3). The exception is a reference species' own gene name: when SwissProt hits of other species than the reference carry the voted name and none carries the reference name ("minidiscs" and "polychaetoid" with a Drosophila reference, where other species say "Y+L amino acid transporter 2" and "Tight junction protein ZO-2"), the hits' name stays. Curated entries of another reference gene (a paralog: UCP2 entries for a UCP3 gene) do not count for the hits' name. When the check changes the symbol to another gene, that gene's reference name is used.
- An RBH symbol that differs from the vote but occurs among the resolved symbols of the winning cluster replaces it (`rbh_override`). A gene the vote left unnamed takes its RBH symbol (`rbh_added`).
- Any other voted reference symbol must be supported by the gene's own sequence: the gene's hit to that symbol's reference protein must score at least 0.9 of its best reference hit (`--symbol-support`) and cover at least half of the reference protein (`--symbol-min-coverage`). Supported votes are kept (`vote_supported`). Unsupported votes are replaced by the closest supported reference protein (`reference_best_hit`; the original is kept in the `voted_symbol` column). When no reference protein is covered, as for a shared domain, a fragment or a lineage-specific gene, the gene's reference alignments decide: when the voted symbol's protein is still the gene's best or near-best reference hit, the reference agrees on the identity and cannot confirm the coverage, so the voted symbol is kept (`vote_unverified`) with 0.1 less confidence; when other reference proteins score higher, or the gene does not hit the voted symbol's protein at all, the symbol is removed (`unsupported_removed`). The product name is kept either way.
- Genes without any reference hit keep their vote (`vote`).
- Voted symbols outside the reference nomenclature (a zebrafish `slc6a4b`, a rodent `Vmn1r237`, an outdated `AARS`) follow `--symbol-policy`. Under `ortholog` (default) a symbol is replaced by a reference symbol only with evidence that it names the same gene, and the reference protein must be covered (`--symbol-min-coverage` of it, or `--symbol-min-query-coverage` of the gene): it is a unique NCBI synonym of a reference symbol (eif2c4 -> AGO4, AARS -> AARS1), or the gene is the reciprocal best hit of the reference protein and the second-best reference protein scores below `--rbh-margin` (0.9) of it (Lmtk1 -> AATK). The margin keeps families expanded in one lineage (olfactory and vomeronasal receptors, MHC, cytochromes P450) from being named after a human paralog. Replaced symbols are `foreign_ortholog`; symbols without evidence are kept (`foreign_kept`), so a gene without a human ortholog still receives the best-supported symbol of another species (fgf24, Vmn1r5). Human pseudogenes (NCBI gene type `pseudo`: CMAHP, VN1R17P) are never assigned by these checks or by ortholog resolution; a functional gene whose human counterpart is a pseudogene keeps its own symbol (Cmah).
- A voted reference symbol is also supported when the gene's alignment covers at least `--symbol-min-query-coverage` (0.8) of the gene itself, so partial gene models keep the symbol of the protein they are part of (RCCD1, ABCA2), while a gene that shares only one domain covers little of either protein.
- In the vote, symbol groups in the reference nomenclature are tried first; when they yield no symbol (a single resolved hit among many), all groups are tried, so a symbol carried by many hits of other species (slc6a4b on 45 of 104 hits) is not lost. All symbol comparisons are case-insensitive.
- **Aliases.** Every symbol gets a `symbol_aliases` list (tables, evidence file, GFF3 `goannotate_symbol_aliases`): the GN= symbols of the hits behind the vote, the replaced voted symbol, and the reference species' NCBI synonyms that read as gene symbols, belong to one gene, and are not disease-locus designations, clone identifiers or pseudogene names (FOXO3: FKHRL1, FOXO2, FOXO3A, foxo3b).

**Reference species by lineage.** NCBI gene_orthologs converts symbols only through its anchor species: human for vertebrates, Drosophila melanogaster for insects, Saccharomyces cerevisiae for fungi, and a few others. Human anchor rows contain no insect genes, so an insect annotated against human keeps the symbols of its hits, mostly genome-project locus tags. With `--reference-taxid auto` (default) the reference is chosen from the lineage of `--taxid` (taxonomy tables of the database, or `--taxdump`). For a reference species with fewer than 10,000 SwissProt proteins (Drosophila: about 3,900) the reference proteome also takes the species' TrEMBL entries. In a test on the paper wasp Polistes fuscatus, Drosophila as reference raised the share of coding genes with a symbol from 66.8% (human) to 73.9%, with 28% of hits converted to the reference against 1%.

**Second nomenclature.** When the reference species is not human, every gene also receives, in the `secondary_symbol` column (tables, evidence file, GFF3 `goannotate_secondary_symbol`), the symbol of its human ortholog when there is an unambiguous one: the gene is the reciprocal best hit of a human protein, that protein is its best human hit, the second-best scores below `--rbh-margin` of it, and the alignment covers the protein or the gene. An insect gene can then carry Dad (FlyBase) and SMAD6 (human). `--secondary-reference-taxid` chooses another species or disables it. The secondary reference also names products: a gene whose product is its primary reference protein's name, used by none of its hits (a reference species' own gene name such as "larval translucida" in Drosophila), takes the name of its secondary reference protein when at least one of its hits uses that name (`product_note` `secondary_reference`).

**Paralog designators.** Other genomes' nomenclatures add designators to the names of their own duplicates: teleost `a`/`b` ("Forkhead box protein O3-B", "protein kinase C beta a"), ", tandem duplicate 1". These mean nothing in the genome being annotated. ", tandem duplicate N" is stripped from hit names like any qualifier. A trailing ` a`/` b`, a letter glued to a trailing number ("O3b") or `-A`/`-B` after a number ("O3-B") is removed from a voted product name only when the name without it is a reference protein name or a name carried by the gene's own hits, so that "Cytochrome b" and "Annexin A" are left alone (`product_note` = `designator_stripped` in the evidence table; `--keep-paralog-designators` disables this).

Several genes may carry the same symbol. Paralogs, recent duplicates and retrocopies are real, and each copy is named for the gene it is homologous to; GOAnnotate does not reserve a symbol for one gene or rename copies with a "-like" suffix. Nucleotide transcripts annotated from their blastx hits alone (`--orf-translation off`) skip this pass, since the reverse search needs protein queries; `--no-reciprocal` disables it.

**Pseudogenes, readthroughs and weakly evidenced paralogs.** A readthrough locus of the reference species (an official symbol made of two official symbols: IQCJ-SCHIP1, KLRC4-KLRK1) never names a gene. A reference pseudogene symbol (NCBI type `pseudo`) does not win the vote when only reference-species hits carry it (FAM90A27P); when hits of other species carry it, it names their functional gene and stays (rodent Uox, Gucy1b2, Fer1l4, whose human orthologs are pseudogenes). Neither is ever a reciprocal best hit or a replacement. After the checks, a symbol whose reference proteins have only weak UniProt protein evidence (PE 2 to 5) is replaced by the established reference gene (PE=1) it is named after (its symbol begins with that gene's symbol) when the gene covers it and hits it at least 0.95 as well: a retrocopy or a recent human duplicate scoring almost like its parent gives way to the parent (NANOGP8 to NANOG, CD8B2 to CD8B). The replaced symbol is kept in `voted_symbol` and among the aliases.

**Ortholog symbols.** With the orthology module, OrthoFinder runs before the vote, the relatives' ortholog symbols enter the vote in Phase 6 (a close relative's symbol decides it, a distant one is one more vote), and they are applied after these checks: a close relatives' reference symbol for another gene is taken over a supported core symbol when the gene's own reference alignments allow it (`ortholog`, decision `taken`, listed for review) unless a single close relative alone names it against a reciprocal best hit (`conflict_kept`); a distant one fills what is left (`ortholog_fill`). See [orthology](#orthology-symbols-from-relatives-weighted-by-their-closeness).

### Phase 9: Confidence

Each gene receives a confidence score in [0, 1] and a label, computed once the annotation is final (after the reference checks and the modules) from the concordance of the final symbol and product with the gene's hits:

```
score = 0.45 * name_concordance
      + 0.25 * min(1, winning_cluster_size / 10)
      + 0.15 * (symbol assigned)
      + 0.15 * symbol_concordance
```

The score is lowered by 0.1 when no hit of the gene passed the hit coverage filter (coverage status `fallback`). With the orthology module, orthology evidence enters the score: a symbol confirmed by close relatives' one-to-one orthologs adds 0.15 (distant relatives 0.05); a symbol taken from close relatives against the hits is placed between 0.55 and 0.69 (medium, for review); a filled symbol scores at least 0.55 (close relatives) or 0.4 (distant), also for genes without hits; a kept conflict (`conflict_kept`) subtracts 0.15. Labels: `high` (>= 0.7), `medium` (>= 0.4), `low` (< 0.4), `none` (no hits). The label is written to the evidence table and to the annotated GFF3 (`goannotate_confidence`).

### Phase 10: Comparison with Existing Annotation

When the GFF3 already carries gene symbols (`gene=` by default) or product names (`product=` or `description=`), GOAnnotate compares them with its own result for every gene. Symbols are compared case-insensitively. Products are compared after normalization (case, punctuation, and qualifiers such as "transcript variant X2", "isoform", "-like", "precursor" are removed) with the same similarity measure used for clustering: `match` (>= 0.8), `partial` (>= 0.5), `mismatch`. Genes are also classified as `goannotate_only`, `existing_only` or `neither`.

The comparison is stratified by any gene- or transcript-level attributes named with `--strata-gene` and `--strata-transcript` (defaults: GRAFT's `graft_status` and `graft_protein`); attributes absent from the file are skipped. Transcript-level values are collapsed per gene, with several distinct values reported as `mixed(a|b)`.

Resolution follows from the annotation mode. In the default `fill` mode the existing annotation is kept wherever it exists, genes without an existing symbol or product (LOC-style genes, novel genes) receive the GOAnnotate values, and every gene where the two symbols both exist and disagree is listed in `{prefix}_review.tsv` with both annotations, the stratification attributes and the confidence, while the GFF3 keeps the existing symbol and records the GOAnnotate one as `goannotate_symbol`. `--annotation-mode replace` takes the GOAnnotate values everywhere instead.

### Hit Support (Concordance) and the Evaluation Report

For each annotated gene with hits, GOAnnotate computes two hit-support (concordance) scores against the gene's full merged hit pool (not just the winning cluster), once the symbol and product are final (after the reference checks, the orthology module and InterProScan), so they describe the annotation that is written:

- **Name concordance:** the fraction of the gene's hits whose cleaned protein name matches the final product (at or above the similarity threshold).
- **Symbol concordance:** the fraction of the gene's hits that carry a symbol and name the gene of the final symbol: the symbol itself, the replaced voted symbol, or one of the gene's aliases (hit symbols of the same gene, reference-species synonyms). Hits without a symbol (many TrEMBL entries) are not counted.

A name concordance of 1.0 means every hit agrees with the product. Low concordance (e.g., < 0.2) indicates the gene's hits are split across divergent annotations, which can happen with multidomain proteins, gene fusions, or genuinely ambiguous homology, or that the product came from another source (a close relative's ortholog, an InterPro domain). Genes without hits (symbols filled from close relatives) get 0 and are left out of the summary statistics. The confidence score is computed from these final values.

**Evaluation report.** `reports/{prefix}_evaluation.tsv` scores every gene's final symbol and product against each source of evidence independently, without changing anything: the gene's own hits (`hits_symbol`, `hits_name`: agree at half the hits or more, minority, disagree), the reference proteome (`reference`: the symbol's reference protein is the gene's best reference hit, near_best within 90% of it, below, absent when the gene hits other reference proteins but not this one, not_reference for a lineage symbol, none), the relatives' one-to-one orthologs and the same-species annotation (`relatives`, `same_species`, with the symbols), the InterPro entries (`domains`: the product shares an informative word with an entry description; entry descriptions name domains and families, so a shared word supports the product and `no_shared_word` says nothing against it), and the product against the symbol (`product_symbol`: the product names the gene of the symbol by the reference protein's name or the curated hits carrying the symbol). `product_generic` marks family-level products, which name a family rather than the symbol's gene. `decided_by_relatives` marks genes whose symbol the relatives decided: the orthology prior picked a cluster or symbol group that the hits had scored lower, or the orthology module took or filled the symbol; the relatives' agreement is then read separately for the genes they did not decide. `conflict` lists the sources that contradict a gene's symbol; the summary counts each source's agreement and the conflicts by kind. This report is the check that the decisions made along the way are supported by the evidence, and the list of genes to review.

## Output Files

The annotation files are written to the top level of the `--output` directory; reports and intermediate files go to subdirectories. `{prefix}` is `--prefix`, or the name of the output directory.

```
{output}/
    {prefix}_annotated.gff3                  input GFF3 with GOAnnotate results
    {prefix}_annotated.gtf                   input GTF with GOAnnotate results (with --gtf)
    {prefix}_annotated_AA.faa                proteins with annotated headers
    {prefix}_annotated_CDS.fna               coding sequences with annotated headers (genome extraction)
    {prefix}_annotated_transcripts.fna       spliced transcripts with annotated headers (genome extraction,
                                             or nucleotide --transcripts input, with ORF coordinates)
    {prefix}_gene_annotations.tsv            one row per gene: symbol, product, locus, ..., GO terms
    {prefix}_transcript_annotations.tsv      one row per transcript
    {prefix}_gene2symbol.tsv, _gene2name.tsv, _gene2go.tsv
    {prefix}_summary.txt
    reports/
        {prefix}_comparison.tsv, _comparison_summary.txt, _review.tsv   comparison with the input annotation
        {prefix}_annotation_evidence.tsv     per-gene evidence for every decision
        {prefix}_go_evidence.tsv             per gene and GO term: evidence codes, class, supporting hits,
                                             source (hits, reference gene, InterPro); terms removed by the
                                             taxon constraints
        {prefix}_name_audit.txt              what the name and symbol filters removed and kept
        {prefix}_evaluation.tsv              each gene's final symbol and product against every source of evidence
        {prefix}_extraction_report.tsv       per-transcript extraction flags (genome mode)
        {prefix}_orf_translation.tsv         strand, frame, ORF and completeness of each translated transcript
        {prefix}_protein_input_check.tsv     problems found in a --proteins file
    intermediate/
        {prefix}_AA.faa, _CDS.fna, _transcripts.fna   extracted sequences (genome mode); _AA.faa holds the
                                             translated proteins of nucleotide transcripts
        {prefix}_blastx_diamond_results.tsv  blastx results of nucleotide transcripts (before translation)
        {prefix}_diamond_SwissProt.tsv, _diamond_TrEMBL.tsv, _diamond_results.tsv
        {prefix}_reference_proteome.faa, _vs_reference.tsv, _reference_vs_query.tsv, *.dmnd
```

The files in `intermediate/` are reused by later runs into the same output directory (`--force-search` repeats the searches). An output directory written by an earlier version, with everything at the top level, is reorganized on the next run: search and reference files move to `intermediate/`, outdated sequence and report files to `intermediate/old_layout/`.

### Annotation Tables

Two tab-delimited tables for everyday lookup (grep, a spreadsheet, R): one row per gene and one row per transcript, identifiers first and the specific GO terms (as in gene2go, semicolon-delimited) in the last column. Every gene and transcript of the input is listed, annotated or not.

| File | Columns |
|------|---------|
| `{prefix}_gene_annotations.tsv` | gene_id, symbol, symbol_aliases, secondary_symbol, product, gene_type, coding, locus, n_transcripts, transcripts, confidence, symbol_source, product_source, uniprot_accession, go_terms, then one column per module that ran |
| `{prefix}_transcript_annotations.tsv` | gene_id, transcript_id, symbol, symbol_aliases, secondary_symbol, product, transcript_type, coding, locus, cds_locus, protein_length, confidence, symbol_source, product_source, uniprot_accession, go_terms, then the module columns |

`symbol` and `product` follow `--annotation-mode`, as in the annotated GFF3 and FASTA: with `fill` (default) and `add` the input annotation is kept where it exists and GOAnnotate's result fills the rest; with `replace` GOAnnotate's result is used wherever there is one. `symbol_source` and `product_source` say which (`input` or `goannotate`). The symbol and the product name the same gene: where GOAnnotate gives the symbol, its product is used unless the input product agrees with it (at least half of the informative words of the shorter name) or GOAnnotate's product is family-level (a `generic_regex:` name such as "ig-like domain-containing protein"), in which cases the input product is kept; where the input symbol is kept, GOAnnotate's product fills a missing product only when GOAnnotate names the same symbol or none. Transcript rows keep the input transcript's own product, including its isoform label, under `fill` and `add` when the gene keeps the input product; the gene row uses the gene description or the first coding transcript's product without the isoform label. Loci are 1-based and inclusive (`chr1:23246-24202(-)`); `protein_length` is the length of the searched protein. `uniprot_accession` is the best hit of the winning symbol group. Without a GFF3 (transcript FASTA input) the rows come from the FASTA, with input values read from `symbol=`, `product=`, `type=` and `locus=` header fields.

### Simple Mapping Files

These are two-column tab-delimited files. Every gene in the universe is included. Genes without annotations have an empty second column.

| File | Columns | Description |
|------|---------|-------------|
| `{prefix}_gene2symbol.tsv` | gene_id, gene_symbol | Gene ID to gene symbol mapping |
| `{prefix}_gene2name.tsv` | gene_id, protein_name | Gene ID to protein name mapping |
| `{prefix}_gene2go.tsv` | gene_id, go_terms | Gene ID to GO terms (semicolon-delimited) |

Product names in all outputs are cleaned for submission: database prefixes, bracketed qualifiers and evidence codes are removed, unbalanced parentheses and non-ASCII characters (Greek letters are spelled out) are fixed, GFF3 separators are replaced, trailing punctuation is removed, and the first letter is lowercased unless the first word is an acronym (NCBI convention: "olfactory receptor 4S2", "DNA ligase 1"). The `product_regex:` rules of the pattern file then apply the INSDC/NCBI naming conventions that UniProt names do not always follow: EC and TC numbers and parenthesized molecular weights are removed, species names are removed ("(Homo sapiens)", "Drosophila homolog" -> "homolog"), "<SYMBOL> protein" becomes "protein <SYMBOL>", doubled words and a trailing "precursor" are removed. The rules that fired are counted in the summary; `--no-product-cleanup` disables them.

### Annotation Evidence File

`reports/{prefix}_annotation_evidence.tsv` contains detailed per-gene evidence for every annotation decision. This is the file to look at if you want to understand why a particular gene received a particular annotation.

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| gene_symbol | Assigned gene symbol (empty if none) |
| protein_name | Assigned protein name (empty if none) |
| n_isoforms | Number of transcripts that contributed hits |
| n_total_hits | Total merged hits for this gene |
| n_cluster_hits | Number of hits in the winning cluster |
| winning_cluster_names | Protein names from winning cluster hits (pipe-delimited) |
| winning_cluster_symbols | Symbols from winning cluster hits (pipe-delimited) |
| top_hit_sseqid | Subject ID of the best hit in the winning symbol group |
| top_hit_bitscore | Bitscore of that hit |
| top_hit_source | `sp` (SwissProt) or `tr` (TrEMBL) |
| mean_cluster_bitscore | Mean bitscore of hits in the winning cluster |
| name_concordance | Fraction of the gene's hits whose name matches the final product |
| symbol_concordance | Fraction of the gene's hits with a symbol that name the final symbol's gene |
| go_source | How GO terms were collected (`cluster_top5`, `cluster_top5+swissprot_supplement`, or empty) |
| n_go_terms | Number of specific GO terms assigned |
| confidence | `high`, `medium`, `low` or `none` |
| confidence_score | Numeric confidence in [0, 1] |
| gene_type | Gene biotype or feature type from the GFF3 (genome mode) |
| coding | `yes` or `no` (genome mode) |
| symbol_status | Outcome of the reference checks: `rbh`, `rbh_override`, `rbh_added`, `vote_supported`, `reference_best_hit`, `vote_unverified`, `unsupported_removed`, `foreign_ortholog`, `foreign_kept`, `foreign_replaced`, `foreign_removed`, `vote`, `vote_rbh_conflict`, `propagated`, `none` |
| voted_symbol | The voted symbol, when the reference check replaced or removed it |
| coverage_status | `pass` (the annotation rests on hits that passed the coverage filter), `fallback` (no hit passed; unfiltered hits used), empty when the filter is off |
| best_query_coverage, best_subject_coverage | Best fraction of the query and of the database protein spanned by a hit of the winning symbol group |
| go_evidence_classes | Specific GO terms by strongest evidence class (`experimental:12,electronic:30`) |
| n_go_taxon_removed | GO terms removed by the taxon constraints |
| product_note | `designator_stripped` when a paralog designator was removed from the product |

### Comparison Files

`reports/{prefix}_comparison.tsv` has one row per gene with the existing and GOAnnotate symbol and product, their status, the stratification attributes, the confidence and the number of GO terms. `{prefix}_comparison_summary.txt` tabulates agreement overall and by stratum. `{prefix}_review.tsv` holds the rows where the existing and GOAnnotate symbols both exist and differ. All three are written only when the GFF3 carries existing annotation.

### Annotated Sequence Files

| File | Description |
|------|-------------|
| `{prefix}_annotated_AA.faa` | The searched proteins (extracted, `--proteins`, protein `--transcripts`, or the proteins translated from nucleotide transcripts, whose headers carry `orf=START-END strand= frame= orf_status=`) |
| `{prefix}_annotated_CDS.fna` | Extracted coding sequences (genome mode) |
| `{prefix}_annotated_transcripts.fna` | Extracted spliced transcripts, all types (genome mode), or the nucleotide `--transcripts` input, with the ORF fields on translated transcripts |

The FASTA files carry GOAnnotate's symbol and product in the header format described under Input Modes; `--annotation-mode` applies as for the GFF3 (`add`: `goannotate_symbol=` and `goannotate_product=` appended; `fill`: `symbol=` and `product=` set where absent; `replace`: overwritten).

| File | Description |
|------|-------------|
| `{prefix}_annotated.gff3` | The input GFF3 with GOAnnotate results added (genome mode, or `--gff3` with `--transcripts`) |

The annotated GFF3 is the input file passed through line by line: every feature, attribute, comment and `###` separator is kept, including non-coding and pseudogene models. Gene features receive `goannotate_symbol`, `goannotate_product`, `goannotate_confidence`, `Ontology_term` (the specific GO terms) and a `Dbxref=UniProtKB:ACC` for the top hit. `--annotation-mode` controls the standard attributes:

| Mode | Behaviour |
|------|-----------|
| `add` | Only the attributes above are added |
| `fill` (default) | `gene`, `Name` and `description` on genes and `product` on transcripts and CDS are written where the input has none, and, where GOAnnotate gives the symbol, an existing product that names another protein is replaced by GOAnnotate's so that symbol and product name the same gene; with `--hypothetical`, coding models still without a product receive `product=hypothetical protein` (off by default, since table2asn assigns it itself). Existing symbols that are identifiers (LOC numbers of an NCBI annotation, locus tags) and placeholder products ("uncharacterized protein LOC122519524") count as missing and are filled; `--keep-existing-identifiers` keeps them |
| `replace` | The same attributes are overwritten with GOAnnotate values wherever GOAnnotate has one |

With `--gtf` the file is `{prefix}_annotated.gtf`, the input GTF passed through in the same way; lines that receive nothing are written unchanged. A gene's results (`goannotate_symbol`, `goannotate_product`, `goannotate_confidence`, one `Ontology_term "GO:..."` per specific GO term, `db_xref "UniProtKB:ACC"`, module attributes) are written on its gene line. A gene without a gene line has them on its transcript lines, and a transcript without a line of its own on its exon and CDS lines. The modes act as for GFF3 with the GTF keys: the symbol is written on every line of the gene, under `gene` when the input uses that key and under `gene_name` otherwise; `description` is written on the gene line and `product` on transcript and CDS lines. List values of modules (`Dbxref`) are written as one `db_xref` per value, and `sig_peptide` features as GTF lines with `gene_id` and `transcript_id`.

### Summary

`{prefix}_summary.txt` contains a text summary of the run: command, inputs, parameters, GFF3 and extraction statistics, Diamond statistics, annotation coverage, confidence, GO term statistics, concordance statistics, and the comparison summary.

## Optional Modules

Modules wrap external tools that add information to the annotation. The orthology module weighs the symbols of close relatives' orthologs against the hit-based vote (see below); SignalP adds signal peptide calls and no symbols, products or GO terms; the InterProScan module adds GO terms and names genes left without a product. They are off by default and selected with `--modules NAME[,NAME]`, or `--modules auto` to run every module whose tool is found. The tools are installed separately by the user; `goannotate --check-tools` reports which are available. The orthology module runs OrthoFinder before the vote, so that the relatives' symbols can steer it, and applies its symbols after the symbol checks; SignalP and InterProScan run after the symbol checks and the reference-gene GO terms, so that their reports list the final symbol and product of every gene, and a gene named from its orthologs is not given a domain-based name. Their raw outputs are kept in `intermediate/{prefix}_{module}/` and reused by later runs into the same directory (`--force-modules` reruns them). Modules need protein queries and are skipped for nucleotide transcripts annotated from their blastx hits alone (`--orf-translation off`).

Each module appends its columns to the end of the gene and transcript tables, adds `goannotate_*` attributes (or standard INSDC features) to the annotated GFF3, writes a report to `reports/`, and is summarized in a Modules section of the summary file.

### signalp: signal peptides and the secretome

Runs [SignalP 6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (Teufel et al. 2022, Nature Biotechnology 40:1023) on the query proteins in eukarya mode and records, for every protein, whether a Sec/SPI signal peptide is predicted, its probability and the cleavage site.

SignalP 6 is distributed under an academic licence: register at the DTU page, download the package (the `fast` model at least; the `slow` models are a separate download), install it with `pip install signalp-6-package/` and copy its model weights into the installed package as the SignalP installation notes describe. The package supports Python up to 3.10 and PyTorch below 2.0, so it will usually live in its own conda environment; point `--signalp-bin` at that environment's `signalp6` executable (for example `~/miniconda3/envs/signalp6/bin/signalp6`). The module is skipped with a note when the executable is not found.

SignalP reads only the first 70 residues of a protein. GOAnnotate therefore sends it only those 70 residues (stop codons and non-standard letters replaced), predicts each distinct 70-residue N-terminus once, and copies the result to every transcript that shares it; isoforms often do (about a third fewer predictions for the prairie vole). The predictions run as parallel `signalp6` processes on chunks of the input, each with one PyTorch thread: one process per physical core by default, limited by `--threads` and by available memory (about 2.5 GB per process in `fast` mode), or the number set with `--signalp-processes`. Progress is written to `signalp.log` in each chunk directory and summarized in the GOAnnotate log every minute.

`--signalp-mode` chooses the model: `fast` (default), `slow` (the full model, needs more than 14 GB RAM per process) or `slow-sequential` (the full model within the memory of `fast`, about six times slower). A protein is called as having a signal peptide when SignalP calls one; `--signalp-min-probability P` additionally requires a SignalP probability of at least P.

Runtime guide (`fast` mode, Threadripper 3970X, 32 cores): about 15 distinct N-termini per second with 32 processes of one thread each, against 9.6 per second for one process with 32 threads; about 30 minutes for the prairie vole (43,357 proteins, 27,685 distinct N-termini). Predictions are identical for every process and thread arrangement.

**Alternative starts (`--signalp-alt-starts`).** A protein without a signal peptide at its annotated start may have one from a later start codon when the gene model begins too far upstream. With this option, every such protein whose first 60 residues contain an internal methionine followed by a signal-peptide-like hydrophobic core (6 of 8 consecutive hydrophobic residues within the next 30; 98% of the signal peptides SignalP calls in the prairie vole pass this test from their own start) is predicted again from each such methionine, and the first one with a signal peptide is reported. The gate depends on the sequence only, so it also finds secreted proteins without known homologs. It adds about a third to the SignalP runtime (prairie vole: 7,941 additional distinct starts, about 9 minutes on 32 cores). A signal peptide found this way is reported as a flag in `sp_alternative_start` and in a separate report; it does not change the main call or the gene model. Where the transcript's best hit that aligns from the database protein's first residues starts its alignment at or after that methionine, the flag notes this homolog (the residues before the methionine have no counterpart in it).

Outputs:

| Where | Content |
|-------|---------|
| gene table | `signal_peptide`: `yes`, `no`, `yes (1/3 isoforms)` when only some coding isoforms carry one, or `not_assessable` when no isoform has its N-terminus (a 5' partial transcript, a gene model without a start codon: a negative call on such a protein says nothing); empty for genes without a scored protein |
| transcript table | `signal_peptide` (`yes`/`no`/`not_assessable`), `sp_probability`, `sp_cleavage` (`22-23`: the signal peptidase cleaves between residues 22 and 23); with `--signalp-alt-starts`, `sp_alternative_start` (`M34: SP 0.99, cleavage 55-56; homolog P99999 aligns from residue 34`, positions in the annotated protein) |
| annotated GFF3 | `goannotate_signal_peptide=yes` on genes and transcripts with a signal peptide, `goannotate_sp_cleavage` on the transcript, `goannotate_sp_alternative_start=M34` on transcripts flagged by `--signalp-alt-starts`; with `--signalp-features`, a `sig_peptide` feature (INSDC type, accepted by table2asn) under each such mRNA, with genomic coordinates mapped through the CDS segments (one line per segment when the signal peptide spans an intron) |
| `reports/{prefix}_secretome.tsv` | gene_id, transcript_id, symbol, product, sp_probability, cleavage_site, cs_probability for every protein with a signal peptide |
| `reports/{prefix}_signalp_alternative_starts.tsv` | with `--signalp-alt-starts`: gene_id, transcript_id, symbol, product, sp_probability_annotated_start, alternative_start, sp_probability_alternative_start, cleavage_site, homology_support |
| `intermediate/{prefix}_signalp/` | `annotated_starts/` (and `alternative_starts/`): the merged `prediction_results.txt` and one directory per chunk with SignalP's own outputs and its log |

No GO terms are added: a signal peptide is not a GO assertion by itself.

### orthology: symbols from relatives, weighted by their closeness

Finds the orthologs of every gene in the RefSeq annotations of relatives with [OrthoFinder](https://github.com/davidemms/OrthoFinder) (Emms & Kelly 2019, Genome Biology 20:238) and weighs their symbols against the hit-based vote according to how close the relatives are. An ortholog in a close, well-named relative is the most direct evidence of a gene's identity: the gene tree separates paralogs that a best hit confuses, and symbols of expanded or lineage-specific families are often present there and missing from UniProt hits. An ortholog in a distant relative is weaker evidence, because one-to-one pairs between distant genomes are often lineage-specific paralogs, so it counts as one more vote and fills gaps.

**Requirements.** OrthoFinder (for example `conda create -n orthofinder -c bioconda orthofinder`; give its executable with `--orthofinder-bin`, and GOAnnotate puts that environment first on PATH while it runs) and the NCBI `datasets` command-line tool on PATH (or `--datasets-bin`). Relatives' annotations are downloaded from NCBI.

**Choosing relatives.** By default they are chosen from the lineage of `--taxid` (which needs the taxonomy tables of the SQLite database), closest first. The lineage is walked outward from the genus; at each level up to 10 NCBI reference genomes with a RefSeq annotation are listed, newest annotation first, the query species itself is skipped (see `--same-species`), and each candidate's GFF3 is downloaded to measure its informative fraction: the protein-coding genes that have a symbol rather than an identifier (LOC numbers and the other symbol patterns), plus those whose NCBI GeneID links through gene_orthologs to a named reference gene. The second part matters for older annotations without symbols: the paper wasp *Polistes canadensis* (RefSeq 2015) has no symbols but is 39.8% informative through its Drosophila orthologs. A candidate is rejected only below a low floor, `--relative-min-named` (default 0.2), or below `--relative-named-ratio` (default 0.75) of the best candidate at its level or closer. The remaining candidates are ranked by informative fraction less `--relative-named-margin` (default 0.10) per lineage level step, so a farther relative is preferred only when it is clearly better named, and the best `--relatives` (default 3) are used. When no candidate is usable, the module is skipped with a warning and the annotation rests on the Diamond vote alone. `--relative` names relatives instead (an assembly accession such as `GCF_000001635.27`, or a taxon ID, whose newest annotated reference genome is used). Every candidate examined, with its level, annotation, symbol and informative fractions and whether it was used (`chosen`, `not_needed`, `rejected_named`, `rejected_ratio`), is listed in `reports/{prefix}_orthology_relatives.tsv`. Downloads are kept in `--relatives-dir` and reused. When NCBI cannot be reached (a `datasets summary` query fails), the relatives chosen in the previous run in the same output directory are used again, provided their annotations are still in `--relatives-dir`; the summary says so.

**Orthologs.** One protein per coding gene of the query (the isoform with the best Diamond hit, the longest among equals or without hits: a mis-extended model aligns worse than a shorter isoform) and the longest protein of each gene of each relative goes to OrthoFinder (Diamond searches, default settings: DendroBLAST gene trees, no RAxML). `--orthofinder-threads` limits its threads (default `--threads`); its progress lines in the log include the CPU temperature when `sensors` is available.

**Closeness (`--ortholog-close-rank`, default `family`).** A relative is close when its lowest common ancestor with the query species, in the NCBI taxonomy of the database, has a rank at or below the given one (an unranked ancestor counts as the first ranked one above it): for the vole, *M. pennsylvanicus* (genus), *Myodes* and *Alexandromys* (subfamily Arvicolinae) are close, mouse and rat (order Rodentia) are distant. The query species' own annotation (`--same-species`) is close. Each relative's common ancestor and closeness are listed in the summary and in `reports/{prefix}_orthology_relatives.tsv`.

**The ortholog symbol.** For every gene the module derives one ortholog symbol: the symbol on which its one-to-one orthologs agree among the relatives that name them (case-insensitive; identifiers such as LOC numbers do not count). Genes whose orthologs are one-to-many, or whose relatives disagree, get none. With transcript input, a gene of a relative that several transcripts share (many-to-one) is used like a one-to-one ortholog when the transcript's own evidence carries the symbol (its vote, one of its hits, or its alignment to the symbol's reference protein scoring at least 90% of its best reference hit): fragments of one gene are grouped this way, and so are transcripts of repeats inserted in a relative's gene, which the transcript's own evidence tells apart (decision `shared_unsupported` when it does not carry the symbol). With gene models, many-to-one orthologs are paralogs and are not used. The summary counts the ortholog symbols that rest on a shared ortholog and those not used. The symbol is supported by close relatives when at least one close relative names it. Official reference symbols of uncharacterized genes (C19orf81, KIAA1549L), which the bad-name rules reject when they come from hits, are accepted from close relatives: the one-to-one ortholog mapped to the reference gene identifies the gene, and the symbol is the official one (products of such genes still follow the product rules). Such symbols only fill genes without a supported symbol; they never replace one.

**Weighting (`--ortholog-symbol weighted`, default).** OrthoFinder runs before the Diamond vote, and the ortholog symbols enter it:

- A symbol supported by close relatives decides the vote: when the gene's hits carry it, the cluster and symbol group that carry it win (the log counts these genes as decided by a close relative's ortholog, and how many of them over a group the hits scored higher; the report marks them `confirmed`). When the hits do not carry it, or the reference check (Phase 8b) chose another symbol, the ortholog symbol may still be taken after Phase 8b, but only when the gene's own sequence allows it: its alignment to the ortholog symbol's reference protein must score at least 90% of its alignment to the core symbol's reference protein. The ortholog symbol of a gene without any reference hit is taken. A taken symbol gives `symbol_status` `ortholog`, the core's symbol kept in `voted_symbol` (not among the aliases, which list other names of the same gene; the vote's aliases are removed as well, since they name the core's gene), the reference protein's name as the product, and decision `taken` in `reports/{prefix}_orthology.tsv`. When the sequence scores the ortholog symbol's reference protein below 90% of the core symbol's, or does not hit it at all, the core symbol stays (decision `conflict_kept`): the relatives' annotation or the orthogroup is then the likelier error, and the gene is listed for review. A second guard covers a single close relative alone naming a symbol that no hit carries against a reciprocal best hit (also `conflict_kept`). The summary counts both.
- The gene's own evidence names the gene; the relatives confirm it or fill gaps. When the relatives name the same reference gene in a form of their own (a lineage name, or a symbol with a copy suffix such as `kdm2aa` for `kdm2a`), the gene keeps the reference symbol and the relatives' form is shown in the orthology report. A symbol of the relatives that has no reference protein cannot be judged by the gene's sequence and does not replace a symbol the gene's reference alignments support (decision `kept`); it still names a gene left without a supported symbol. A taken symbol (another gene) receives the GO terms of its own reference gene, never those of the replaced symbol.
- A symbol supported only by distant relatives is one more vote: in the Diamond vote its symbol group wins when it scores at least 75% of the best group; afterwards it fills genes left without a symbol, or with a symbol outside the reference nomenclature kept for lack of a reference ortholog (`foreign_kept`), as `ortholog_fill`. It never replaces a supported core symbol; where the two differ the gene is listed with decision `kept`.
- Every supported core symbol is compared with the ortholog symbol, and the agreement rate is reported as a check of both (former names of the same reference gene count as agreement).

The product follows the symbol: a gene that takes an ortholog symbol, or is filled with one while the vote left it unnamed, takes the reference protein's name of that symbol when it is a reference symbol, else the relatives' most frequent product (`product_note` `ortholog`); a gene filled with the reference form of its own lineage symbol keeps the hits' product. A symbol given to a gene without a supported one is checked against the gene's own evidence: when none of the gene's hits carries the ortholog symbol, the symbol is a reference symbol, and the gene's alignment to its reference protein scores below 90% of the gene's best reference hit, or the gene does not hit it at all, the gene is left as it was (decision `fill_refused`). A symbol that the hits carry, a lineage symbol, or a gene without any reference hit cannot be judged this way and is filled.

`--ortholog-symbol fill` is the earlier behaviour: no influence on the vote, fills only. `off` records the orthologs and the checks without changing a symbol.

**Reference nomenclature.** Before they are compared, the relatives' symbols are put in the reference nomenclature in the same way as hit symbols: each relative gene's NCBI GeneID (from its GFF3) is mapped through NCBI gene_orthologs to its ortholog in the reference species, whose symbol is used (with a human reference, mouse `Zfp558` and vole `Znf558` both become `ZNF558`). A gene without a reference ortholog keeps its own symbol. The per-gene report shows both forms (`XP_006509297.1:Zfp558>ZNF558`). `--no-resolve-orthologs` turns this off as it does for hits.

**Previous annotation of the same species (`--same-species`).** When the query species already has a RefSeq annotation (of the same assembly or an earlier one), its newest one is by default the first of the `--relatives`, so that its symbols carry over to the new annotation (`--same-species off` leaves it out; without one the run proceeds without it). Its one-to-one ortholog's symbol is used when the other relatives agree with it, name no one-to-one ortholog, or disagree among themselves; when they agree on a different symbol, nothing is filled and the gene is counted among the disagreements. Supported core symbols are also compared with the previous annotation, and the agreement is reported as a continuity check. The other check (supported symbols against the other relatives) does not use the same-species annotation.

Outputs:

| Where | Content |
|-------|---------|
| gene table | `orthogroup`, `relative_orthologs` (`Microtus_pennsylvanicus:Avp(1:1);...`), `ortholog_symbol` (the agreed one-to-one symbol, used or not), `ortholog_support` (`close` or `distant`) |
| annotated GFF3 | `goannotate_orthogroup`, `goannotate_ortholog_symbol` on the gene |
| `reports/{prefix}_orthology.tsv` | per gene: symbol and status, orthogroup, ortholog symbol with its agreement and support (close or distant), the decision (`confirmed`, `taken`, `conflict_kept`, `filled`, `fill_refused`, `shared_unsupported`, `kept`, `recorded`, `none`), the close relatives' own form, the core's previous symbol when one was taken or filled, and per relative the ortholog type and orthologs with their symbols |
| `reports/{prefix}_orthology_relatives.tsv` | every candidate relative examined with its common ancestor with the query species and closeness, and, for those used, genes used and query genes with one-to-one and any orthologs |
| `intermediate/{prefix}_orthology/` | OrthoFinder input and results, its log, and the downloaded annotations |

### interproscan: domains, InterPro GO terms and domain-based names

Runs [InterProScan 6](https://github.com/ebi-pf-team/interproscan6) (EMBL-EBI) on the query proteins and adds, for every gene and transcript, the InterPro entries and Pfam families found, the GO terms of those InterPro entries (InterPro2GO), and, for coding genes that the hit-based annotation left without a product, a name built from their InterPro domains.

**Installation.** InterProScan 6 is a Nextflow workflow whose analyses run in containers. It needs Nextflow 25.10.4 or later (which needs Java 17 or later), one container runtime (Docker, Singularity, Apptainer or Podman), the workflow itself, and the InterPro data for one release. InterProScan 5 ended with release 5.78-109.0 (June 2026) and is not supported.

```bash
# Workflow, pinned at a release
git clone https://github.com/ebi-pf-team/interproscan6.git /path/to/InterProScan
git -C /path/to/InterProScan checkout 6.0.3.0

# Data for one InterPro release, downloaded beforehand (about 34 GB unpacked for InterPro 110.0)
D=/path/to/InterProScan_data
URL=https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0
mkdir -p "$D/.download" && cd "$D/.download"
get() {
    name=$1; ver=$2
    curl -sS -OJ "$URL/$name/$name-$ver.tar.gz" && curl -sS -OJ "$URL/$name/$name-$ver.tar.gz.md5" &&
    md5sum -c "$name-$ver.tar.gz.md5" && tar -C "$D" -zxf "$name-$ver.tar.gz" &&
    rm "$name-$ver.tar.gz" "$name-$ver.tar.gz.md5" && echo "done: $name $ver"
}
export -f get; export D URL
printf '%s\n' "interpro 110.0" "antifam 8.0" "cath 4.3.0" "cdd 3.21" "hamap 2026_01" "ncbifam 20.0" \
    "panther 19.0" "pfam 38.2" "pirsf 3.10" "pirsr 2025_05" "prints 42.0" "prosite 2026_01" "sfld 4" \
    "smart 9.0" "superfamily 1.75" | xargs -P 4 -L 1 bash -c 'get "$0" "$1"'
cd "$D" && rmdir .download
```

The member-database versions for an InterPro release are listed in `interpro/<release>/databases.json`. The data must be in place before GOAnnotate runs: InterProScan 6.0.3.0 can download missing data itself, but its download of databases that share a directory (PROSITE patterns and profiles, CATH-Gene3D and CATH-FunFam) can start a search on an empty directory. GOAnnotate therefore checks that every analysis has its data and does not run InterProScan otherwise.

**Running.** `--modules interproscan --interproscan-workflow /path/to/InterProScan --interproscan-datadir /path/to/InterProScan_data`. GOAnnotate writes the cleaned proteins to `intermediate/{prefix}_interproscan/` and runs Nextflow there (its `work/` directory, logs and run cache stay with the GOAnnotate output; `work/` is deleted after a successful run unless `--interproscan-keep-work` is given), with `--goterms`, the pinned InterPro release, `-resume` (an interrupted run reuses its finished tasks), `--max-workers` set to half of `--threads` (`--interproscan-workers` sets it; analyses computed locally, such as CATH's, use more CPU than they request), and the container profile found (`--interproscan-profile`). `--max-workers` is the CPU budget of Nextflow's local executor. Most InterProScan tasks request one CPU; GOAnnotate passes a small configuration (`goannotate_resources.config` in the same directory) that gives the multi-threaded MobiDB-lite step 4 CPUs and, with Docker or Podman, limits each task's container to the CPUs it requested, so the load stays within `--threads`. Results are reused by later runs into the same output directory. By default InterProScan first looks up precalculated matches at EBI for proteins already in UniParc, sending only sequence checksums; proteins not found are analysed locally. `--interproscan-local` runs every analysis locally (on the 103 test proteins of InterProScan: 7 minutes against 2). `--nxf-home` moves Nextflow's own files out of `~/.nextflow`.

**GO terms (`--interpro-go`).** With `add` (default) the GO terms of a gene's InterPro entries join the terms from its hits, and the most specific terms are recomputed: a term the gene already has, or an ancestor of one of its terms, adds nothing, and a more specific InterPro term replaces its ancestor. Terms outside `--namespace`, root terms and terms excluded by the taxon constraints are not added. With `fill` only genes without GO terms from their hits receive them. InterPro terms carry evidence code IEA; `reports/{prefix}_go_evidence.tsv` gives each term's source (`hits`, `reference`, `interpro`, or combinations such as `hits+interpro`) and the InterPro entries behind it.

**Products (`--interpro-products`).** A coding gene without a product after the hit-based annotation is named from the InterPro entries of type Domain (then Repeat) of its longest protein, in order of position, at most two: `uncharacterized EGF-like domain-containing protein`, `uncharacterized EGF-like and protein kinase domain-containing protein`, `uncharacterized WD40 repeat-containing protein`. Genes with a product are not changed; family and superfamily entries do not give names. Such genes are marked `interpro_domain` in the evidence table's `product_note`.

Outputs:

| Where | Content |
|-------|---------|
| gene table | `interpro` (`IPR000742:EGF-like domain;...`, in order of position), `pfam` (Pfam accessions), `go_from_interpro` (specific GO terms the gene received from InterPro only) |
| transcript table | `interpro`, `pfam` |
| annotated GFF3 | `Dbxref=InterPro:...,Pfam:...` on the gene, merged with existing Dbxref values; InterPro GO terms in `Ontology_term` |
| `reports/{prefix}_interproscan.tsv` | every match: gene_id, transcript_id, analysis, signature, signature_description, start, end, score, interpro, interpro_description, interpro_type, go_terms |
| `intermediate/{prefix}_interproscan/` | the protein input, InterProScan's TSV, the Nextflow console output and `.nextflow.log` |

## The bad_patterns.txt File

`bad_patterns.txt` is the list of rules that decide which text from database hits is unusable in another genome: which words are removed from protein names, which names and gene symbols are rejected, which names count as family-level, and how final product names are tidied. Its purpose is to prevent carryover of designations that belong to another genome or annotation pipeline (LOC numbers, locus tags, clone IDs, chromosome ORF names, the NCBI `-like` suffix) and to keep uninformative names ("uncharacterized protein") from being applied as products.

Each rule is one line of the file. A rule is switched off by putting `#` in front of its line or deleting it, and a rule is added by adding a line. A user who wants zebrafish clone symbols kept, for example, comments out the `symbol_regex:^SI:` line.

| To do this | Use |
|---|---|
| Run with the bundled rules | nothing; `goannotate/data/bad_patterns.txt` is read by default |
| Change or switch off bundled rules | `goannotate --write-bad-patterns my_patterns.txt` writes a copy of the bundled file (an existing file is not overwritten); edit the copy and pass it with `--bad-patterns my_patterns.txt` |
| Add rules and keep the bundled ones | put the new lines in a file of their own and pass it with `--bad-patterns-extra FILE` (repeatable) |
| See what the rules did in a run | `reports/{prefix}_name_audit.txt`: discards by reason, identifiers and qualifiers removed by pattern, the most frequent surviving names (family-level names marked `G`), and surviving names and symbols that look like identifiers |

The pattern file used is named in the log and in the summary of every run. `tests/audit_bad_names.py` writes the audit report for a Diamond results file without running the pipeline.

### Text classes and pattern types

| Prefix | Class | Effect |
|---|---|---|
| `qualifier_regex:` (alias `strip_regex:`) | Genome-specific qualifier | Removed from the name; the rest is used |
| `identifier_regex:` | Genome-specific identifier | Removed from the name as a token; name discarded only if nothing informative remains |
| `placeholder:` | No information | Whole-name match after cleaning (case and punctuation ignored); name discarded |
| `exact:` / `contains:` / `regex:` | No information | Name discarded |
| `generic_regex:` | Family-level descriptor | Kept; used as the product only when no specific name is available |
| `normalize_regex:` | Comparison noise | Ignored when names are compared (e.g. `homolog`); never removed from output |
| `symbol_regex:` | Genome-specific symbol | Symbol blanked, after ortholog resolution was attempted |
| `product_regex:` | Naming convention | `PATTERN => REPLACEMENT` applied to the final product names only (without `=>` the match is removed) |
| `symbol_strip_regex:` | Database-specific symbol suffix | Removed from a hit's gene symbol (FlyBase `kay-RA` becomes `kay`) |
| `alias_reject_regex:` | Synonym that is not a gene name | A reference gene's synonym matching it is not listed as an alias (human disease-locus names such as `CMT2N`) |
| `designator_regex:` | Paralog designator of another genome | Group 1, the name without the designator, replaces the product when it is a reference protein name or a name of the gene's own hits ("claudin 15 a") |
| `ortholog_symbol_regex:` | Exception to `symbol_regex:` | A rejected symbol is accepted when a close relative's one-to-one ortholog maps to a reference gene with that official symbol (`C19orf81`) |
| `interpro_identifier_regex:` | InterPro entry named after an uncharacterized protein | The entry gives no product name (InterProScan module) |

All patterns are case-insensitive; `(?-i:...)` makes part of a pattern case-sensitive. Lines without a prefix are read as in 1.x (6 characters or fewer as `exact:`, longer as `contains:`), with a note in the log.

```
qualifier_regex:[\s-]like\s*$
identifier_regex:\bLOC\d{4,}\b
identifier_regex:\bC(?:\d{1,2}|[XY])(?:H\d{1,2})?orf\d+[A-Z]?\b(?: homolog)?
placeholder:uncharacterized protein
generic_regex:^olfactory receptor$
normalize_regex:\bhomolog\b
symbol_regex:^LOC\d+$
product_regex:^((?-i:[A-Z][A-Z0-9]{2,})(?:-\d+)?) protein$ => protein \1
```

Comments (lines starting with `#`) and blank lines are ignored.

## Parameters

### `--top-n` (default: 100)

The number of informative hits kept per transcript after the name and coverage filters, from the up to `--max-target-seqs` hits per database (100 each from SwissProt and TrEMBL). In tests on the prairie vole and Poecilia parae (2026-10-31), 100 instead of 50 raised the genes with a human-nomenclature symbol from 83.4% to 84.9% (vole) and from 69.8% to 73.4% (P. parae), and the genes with GO terms from 91.7% to 92.2% and from 88.4% to 90.1%, for about 4 minutes more annotation time per genome. `--top-n` cannot exceed the hits the search returned.

### `--consensus-threshold` (default: 0.5)

The similarity threshold for clustering protein names. Two names must have at least this similarity score (0 to 1) to be placed in the same cluster. Lower values produce larger, more inclusive clusters. Higher values produce smaller, more specific clusters. The default of 0.5 works well for most use cases.

### `--evalue` (default: 1e-5)

E-value cutoff passed to Diamond and applied again when results are parsed.

### `--sensitivity` (default: fast)

Diamond's sensitivity mode. `fast` is Diamond's default and is adequate for complete proteins searched against a database containing close relatives; `sensitive` or `very-sensitive` recover more distant homologs at several times the runtime.

### `--namespace` (default: BP MF CC)

Which GO namespaces to include: biological_process (BP), molecular_function (MF), cellular_component (CC). You can restrict to a subset if needed.

## Utility Scripts

### clean_gff3.py

A utility for stripping all existing attributes from a GFF3 file and re-assigning clean, sequentially numbered IDs. Useful for preparing messy GFF3 files before annotation. Parent references are updated to match the new IDs.

```bash
python3 clean_gff3.py input.gff3 output.gff3
```

### tests/evaluate_swissprot_preference.py

Evaluates, on a finished run that was compared with existing annotation, whether taking symbols from SwissProt hits of preferred species would raise agreement. Reports agreement before and after and the genes that would change.

### tests/validate_extraction.py

Compares GOAnnotate's extracted proteins with an independent protein FASTA for the same GFF3 (for example the protein set written by the annotation pipeline), matched by transcript ID, and reports identical, prefix/suffix and differing sequences.

```bash
python3 tests/validate_extraction.py --gff3 Annotations.gff3 --genome Genome.fasta \
    --reference proteins.fasta -o validation
```

## Files

| File | Description |
|------|-------------|
| `GOAnnotate.py` | Wrapper script; runs the package from a checkout |
| `goannotate/` | The package: `cli.py`, `pipeline.py`, `gff3.py`, `sequences.py`, `search.py`, `blast.py`, `filters.py`, `annotate.py`, `assign.py`, `refgo.py` (GO terms from the reference gene), `reference.py`, `go.py`, `taxon.py`, `mappings.py`, `compare.py`, `output.py`, `models.py` |
| `goannotate/modules/` | Optional modules: `base.py` (framework), `signalp.py`, `orthology.py`, `interproscan.py` |
| `goannotate/data/bad_patterns.txt` | Default name, symbol and product patterns |
| `goannotate/data/stop_words.txt` | Default stop words for the name similarity |
| `build_databases.py` | Database download and build script |
| `clean_gff3.py` | GFF3 ID cleanup utility |
| `tests/test_annotation_rules.py` | Checks of the coverage filter, word weighting, GO evidence, taxon constraints, ortholog symbols in the vote, reference-gene GO terms, product rules and symbol policy on small inputs |
| `tests/test_modules.py` | Checks of the module framework and the SignalP module with a stand-in `signalp6` script |
| `tests/test_interproscan.py` | Checks of the InterProScan module with a stand-in `nextflow` script and a miniature InterPro data directory |
| `tests/test_orthology.py` | Checks of the orthology module with stand-in `datasets` and `orthofinder` scripts and a miniature taxonomy |
| `tests/validate_extraction.py` | Extraction validation against an external protein set |

## License

GPL-3.0-or-later. See LICENSE.
