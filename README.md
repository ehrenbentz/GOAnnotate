# GOAnnotate

Functional annotation of gene models using Diamond searches against UniProt.

GOAnnotate assigns protein names, gene symbols, and Gene Ontology (GO) terms to genes from Diamond hits against SwissProt and TrEMBL. Starting from a GFF3 and a genome, it extracts and translates the coding sequences itself, runs the Diamond searches, aggregates hits across isoforms, clusters them by protein name similarity, and uses a voting/scoring system to pick the best annotation for each gene. Genomes that already carry functional annotation (for example from EGAPx) are compared gene by gene with the GOAnnotate result, and the existing annotation is kept in the output.

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
- [The bad_names.txt File](#the-bad_namestxt-file)
- [Parameters](#parameters)
- [Utility Scripts](#utility-scripts)

## Requirements

- Python 3.9+
- [Diamond](https://github.com/bbuchfink/diamond) 2.x (not needed when precomputed results are supplied with `--blast-results`)
- Optional, for the modules: SignalP 6 (see Optional Modules)

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

### Pre-extracted proteins: `--gff3` and `--proteins`

A protein FASTA made from the same GFF3 (for example GRAFT's `<name>_graft_proteins.fasta`) is searched instead of an extraction. It is checked against the GFF3 before use:

- the run stops if an ID occurs twice, if an ID is not a transcript identifier in the GFF3 (after `--transcript-id-attr`), or if a `gene=` header field names a different gene than the transcript's parent in the GFF3;
- coding transcripts without a protein, proteins for transcripts without CDS, and proteins whose length differs from the CDS length by more than one residue are reported (proteins written with internal stops removed are shorter by their number).

With `--genome` as well, the transcript and CDS files are still extracted and every supplied protein is compared with GOAnnotate's own translation. Problems are listed in `reports/{prefix}_protein_input_check.tsv` and counted in the summary.

### Transcriptome mode: `--transcripts`

A nucleotide FASTA (CDS or mRNA) is searched with `diamond blastx`; a protein FASTA is detected automatically and searched with `blastp`. Isoforms are grouped into genes through a `gene=` field in the headers, or through `--tx2gene FILE` (two columns: transcript_id, gene_id). Without either, each sequence is its own gene.

### Precomputed results: `--transcripts` and `--blast-results`

Skips the search and annotates existing Diamond or BLAST results (format 6 with the 13 fields `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle`). Adding `--gff3` produces the annotated GFF3 output as well. This is the mode of earlier GOAnnotate versions and is unchanged.

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
    --tx2gene Trinity.gene_trans_map.tsv \
    --diamond-swissprot SwissProt.dmnd \
    --diamond-trembl TrEMBL.dmnd \
    --db GOAnnotate_db.sqlite \
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
| `--gff3` + `--genome`, `--gff3` + `--proteins`, or `--transcripts` | Input (see Input Modes) |
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
| `--bad-names` | bundled file | Name, symbol and product pattern file |
| `--bad-names-extra` | none | Additional pattern file layered on the bundled one (repeatable) |
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
| `--symbol-score` | mean_sqrt | Symbol group score within the winning cluster: mean bitscore, or mean times the square root of the group size |
| `--reference-proteome` | from SwissProt DB | Reference-species protein FASTA for reciprocal best hits |
| `--no-reciprocal` | off | Skip the checks of voted symbols against the reference proteome |
| `--symbol-support` | 0.9 | Minimum score of the voted symbol's reference hit relative to the gene's best reference hit |
| `--symbol-min-coverage` | 0.5 | Minimum fraction of the reference protein covered by the gene's alignment |
| `--symbol-policy` | ortholog | Voted symbols outside the reference nomenclature: `ortholog` replaces them by a reference symbol only with evidence that it is the same gene (suffix-free form, unique NCBI synonym, unambiguous reciprocal best hit) and keeps the rest; `reference` uses the closest covered reference protein or removes them; `any` keeps them |
| `--rbh-margin` | 0.9 | A reciprocal best hit is ortholog evidence only when the second-best reference protein scores below this fraction of the best; 0 accepts all |
| `--drop-unmatched-symbols` | off | With `ortholog`, remove non-reference symbols without ortholog evidence |
| `--symbol-min-query-coverage` | 0.8 | A voted reference symbol is also supported when the alignment covers this fraction of the gene (partial gene models) |
| `--keep-paralog-designators` | off | Keep other genomes' paralog designators (`O3-B`, `beta a`) in voted product names |
| `--annotation-mode` | fill | How results enter the annotated GFF3 (see Output Files) |
| `--keep-existing-identifiers` | off | In fill and add modes, keep existing LOC-type symbols and placeholder products instead of filling them |
| `--hypothetical` | off | Write `product=hypothetical protein` for unnamed coding models |
| `--no-product-cleanup` | off | Do not apply the `product_regex:` naming rules to products |
| `--no-compare` | off | Skip the comparison with existing annotation |
| `--existing-symbol-attr` | gene | Attribute(s) holding an existing gene symbol |
| `--existing-product-attr` | product description | Attribute(s) holding an existing product name |
| `--strata-gene`, `--strata-transcript` | graft_status, graft_protein | Attributes the comparison is stratified by |
| `--modules` | none | Optional modules to run (`signalp`; `auto` runs every module whose tool is found); see Optional Modules |
| `--check-tools` | off | Report which external tools are available, then exit |
| `--force-modules` | off | Rerun the module tools even when their results exist |
| `--signalp-bin`, `--signalp-mode`, `--signalp-features` | `signalp6` on PATH, fast, off | SignalP module settings |
| `--threads` | all CPUs | Threads for Diamond, workers for annotation and threads for the module tools |
| `--prefix` | output dir name | Prefix for output filenames |

## How It Works

### Phase 0: Input Parsing

In genome mode, the GFF3 defines the gene universe and the sequences are extracted and translated as described under Input Modes. In the FASTA modes, GOAnnotate reads the headers: the transcript ID is the first whitespace-delimited field and a `gene=` tag groups transcripts into genes. Every gene appears in the output files, whether or not it has hits.

### Phase 1: Diamond Search

Proteins (or transcripts) are searched against SwissProt and TrEMBL in two separate runs and the results are concatenated. Running the databases separately matters because `--max-target-seqs` limits hits per database: in a combined database a query would often receive 100 TrEMBL hits and no SwissProt hits, losing the curated SwissProt annotations that carry the best GO terms. Hits are parsed and filtered by e-value, keyed by transcript ID. Each hit's subject ID is checked for `sp|` or `tr|` prefixes to identify SwissProt vs. TrEMBL origin. The protein name and gene symbol (from the `GN=` field in the stitle) are extracted from each hit.

### Phase 2: Hit Cleaning and Top-N Selection

Each hit's protein name is cleaned and tested with the patterns in `bad_names.txt`. The purpose is to stop designations that belong to another genome or annotation pipeline from being carried into this one, and to keep names that carry no information from being used as products:

1. **Qualifiers** (`qualifier_regex:`) are stripped: database prefixes, `PREDICTED:`, `LOW QUALITY PROTEIN:`, `(Fragment)`, isoform and variant designations, `putative`/`probable`, and a trailing `-like` (the NCBI pipeline's paralog marker for one genome).
2. **Identifiers** (`identifier_regex:`) are removed as tokens: LOC numbers, RIKEN clone IDs, C#orf# names in short and long form, KIAA/FLJ/MGC/DKFZ clone IDs, UPF and DUF numbers, Gm numbers, clone accessions. "Bombesin receptor-activated protein C6orf89" becomes "Bombesin receptor-activated protein".
3. The remaining name is **discarded** when it is a placeholder (`placeholder:`, e.g. "uncharacterized protein", "protein") or matches an `exact:`, `contains:` or `regex:` pattern. A hit whose name is discarded contributes nothing.

Gene symbols (`GN=`) are carried as found until ortholog resolution has been attempted, so that a LOC number or species-numbered symbol (`Olr1468`, `Zfp1006`) can still reach the human symbol of its gene. Only symbols that fail to resolve and match a `symbol_regex:` pattern are blanked.

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

Within the winning cluster, hits are grouped by gene symbol. The symbol group with the highest score is selected, and both the gene symbol and protein name come from that group. The default score is the mean bitscore of the group times the square root of the group size (`--symbol-score mean_sqrt`), so that agreement among many species counts once symbols are in one nomenclature; `--symbol-score mean` uses the mean alone, which favors small groups of the closest species.

A **single-hit guard** prevents annotation errors in large clusters: if a symbol group has only one hit and there are other groups with two or more hits, the single-hit group cannot win. This stops a lone misannotated hit from overriding the majority. An additional guard rejects a single-hit symbol that represents less than 10% of a cluster with 10+ hits.

If no valid gene symbol exists in the winning cluster, GOAnnotate falls back to NCBI cross-referencing: it looks up the best hit's UniProt accession in the NCBI ID mapping to find a GeneID, then looks up that GeneID in NCBI's gene_info to find a symbol. This step requires either the `--db` SQLite database or the `--ncbi-idmapping` and `--ncbi-geneinfo` flat files.

The protein name comes from the highest-bitscore hit in the winning symbol group. Output symbols are written in one casing convention (`--symbol-case`, default upper).

### Phase 7: GO Term Collection

GO terms are collected in two steps:

1. **Primary collection:** The top 5 hits (by bitscore) from the winning symbol group are looked up in the GO mapping. Their GO terms are pooled together.

2. **SwissProt supplement:** The pipeline then scans all remaining hits in the winning cluster for SwissProt entries that share the winning gene symbol or protein name. GO terms from matching SwissProt hits are added to the pool. This is important because TrEMBL hits typically dominate the bitscore rankings (especially for non-model organisms with closely related species in TrEMBL), but SwissProt entries carry richer, manually curated GO annotations.

**Evidence.** Databases built by the current `build_databases.py` store the GO evidence codes behind every accession-term pair. `--go-evidence` chooses which annotations are used: `all` (default; every code except ND, whose terms are the ontology roots), `non-electronic` (no IEA annotations), or `experimental` (EXP, IDA, IPI, IMP, IGI, IEP and the high-throughput codes). `--go-iea-min-hits N` additionally requires a term whose only evidence is electronic to be carried by at least N of the hits it is collected from. The codes behind each assigned term, its strongest evidence class and the number of supporting hits are written to `reports/{prefix}_go_evidence.tsv`; the evidence table summarizes them per gene (`go_evidence_classes`).

**Taxon constraints.** The Gene Ontology states for many terms the taxa they apply to (`in_taxon`, `never_in_taxon`; constraints pass to subtypes and parts). With `--taxid` (the NCBI taxon of the organism being annotated) GOAnnotate computes the organism's lineage from the taxonomy tables of the database (or `--taxdump`), expands GO's taxon union classes, and removes every GO term that a constraint on it or on one of its ancestors excludes for that lineage: photosynthesis terms from a mammal, lactation from a fish. The constraints come from the database tables or from `--taxon-constraints-obo` and `--taxon-groupings-obo`; removed terms are listed in `reports/{prefix}_go_evidence.tsv` with the constraint that removed them. `--taxon-constraints off` disables the step; without `--taxid` it is skipped.

After collection, GO terms are filtered by namespace (BP, MF, CC) and then reduced to the most specific terms by removing redundant ancestor terms using the GO hierarchy from the OBO file. For example, if a gene has both GO:0005634 (nucleus) and GO:0005737 (cytoplasm), both are kept because neither is an ancestor of the other. But if it has both GO:0005634 (nucleus) and GO:0005622 (intracellular anatomical structure, a parent of nucleus), the parent is removed.

### Phase 8: Symbol Propagation

After all genes are annotated independently, a post-annotation pass looks for genes that share the same normalized protein name. If some of those genes have a gene symbol and others do not, the most common symbol is propagated to a symbol-less sibling when that symbol also occurs among the sibling's own hits. This fills in gaps where the vote for a particular gene left it without a symbol although its hits carry one. Identical product names alone are not used: TrEMBL names such as "Disks large homolog 5-like" are shared by unrelated genes that resemble one domain of the named protein.

### Phase 8b: Symbol Checks Against the Reference Proteome

The vote decides identity from many database hits, but its naming can be wrong in two ways: a gene that shares a domain with a well-studied protein can be named after it, and for paralog pairs the better-studied sibling tends to win the vote. A second pass checks every voted symbol against the complete reference-species proteome (the SwissProt entries of the reference taxon, which for human cover essentially every protein-coding gene; extracted from the SwissProt Diamond database with `diamond getseq`, or given with `--reference-proteome`). The reference proteome is searched against the query proteins and the query proteins against it; both searches take seconds. A gene is the reciprocal best hit (RBH) of reference symbol X when X's protein has that gene as its best query hit and X is among the gene's hits.

- A voted symbol equal to the RBH symbol is confirmed (`rbh`), and the reference protein's name becomes the product.
- An RBH symbol that differs from the vote but occurs among the resolved symbols of the winning cluster replaces it (`rbh_override`). A gene the vote left unnamed takes its RBH symbol (`rbh_added`).
- Any other voted reference symbol must be supported by the gene's own sequence: the gene's hit to that symbol's reference protein must score at least 0.9 of its best reference hit (`--symbol-support`) and cover at least half of the reference protein (`--symbol-min-coverage`). Supported votes are kept (`vote_supported`). Unsupported votes are replaced by the closest supported reference protein (`reference_best_hit`; the original is kept in the `voted_symbol` column), or removed when no reference protein is covered, as for a shared domain or a fragment (`unsupported_removed`); the product name is kept either way.
- Genes without any reference hit keep their vote (`vote`).
- Voted symbols outside the reference nomenclature (a zebrafish `slc6a4b`, a rodent `Vmn1r237`, an outdated `AARS`) follow `--symbol-policy`. Under `ortholog` (default) a symbol is replaced by a reference symbol only with evidence that it names the same gene, and the reference protein must be covered (`--symbol-min-coverage` of it, or `--symbol-min-query-coverage` of the gene): the symbol without a genome suffix (`.2`, `a`, `b`, `l`) is a reference symbol (slc6a4b -> SLC6A4), it is a unique NCBI synonym of one (eif2c4 -> AGO4, AARS -> AARS1), or the gene is the reciprocal best hit of the reference protein and the second-best reference protein scores below `--rbh-margin` (0.9) of it (Lmtk1 -> AATK). The margin keeps families expanded in one lineage (olfactory and vomeronasal receptors, MHC, cytochromes P450) from being named after a human paralog. Replaced symbols are `foreign_ortholog`; symbols without evidence are kept (`foreign_kept`), so a gene without a human ortholog still receives the best-supported symbol of another species (fgf24, Vmn1r5). Human pseudogenes (NCBI gene type `pseudo`: CMAHP, VN1R17P) are never assigned by these checks or by ortholog resolution; a functional gene whose human counterpart is a pseudogene keeps its own symbol (Cmah).
- A voted reference symbol is also supported when the gene's alignment covers at least `--symbol-min-query-coverage` (0.8) of the gene itself, so partial gene models keep the symbol of the protein they are part of (RCCD1, ABCA2), while a gene that shares only one domain covers little of either protein.
- In the vote, symbol groups in the reference nomenclature are tried first; when they yield no symbol (a single resolved hit among many), all groups are tried, so a symbol carried by many hits of other species (slc6a4b on 45 of 104 hits) is not lost. All symbol comparisons are case-insensitive.
- **Aliases.** Every symbol gets a `symbol_aliases` list (tables, evidence file, GFF3 `goannotate_symbol_aliases`): the GN= symbols of the hits behind the vote, the replaced voted symbol, and the reference species' NCBI synonyms that read as gene symbols, belong to one gene, and are not disease-locus designations, clone identifiers or pseudogene names (FOXO3: FKHRL1, FOXO2, FOXO3A, foxo3b).

**Reference species by lineage.** NCBI gene_orthologs converts symbols only through its anchor species: human for vertebrates, Drosophila melanogaster for insects, Saccharomyces cerevisiae for fungi, and a few others. Human anchor rows contain no insect genes, so an insect annotated against human keeps the symbols of its hits, mostly genome-project locus tags. With `--reference-taxid auto` (default) the reference is chosen from the lineage of `--taxid` (taxonomy tables of the database, or `--taxdump`). For a reference species with fewer than 10,000 SwissProt proteins (Drosophila: about 3,900) the reference proteome also takes the species' TrEMBL entries. In a test on the paper wasp Polistes fuscatus, Drosophila as reference raised the share of coding genes with a symbol from 66.8% (human) to 73.9%, with 28% of hits converted to the reference against 1%.

**Second nomenclature.** When the reference species is not human, every gene also receives, in the `secondary_symbol` column (tables, evidence file, GFF3 `goannotate_secondary_symbol`), the symbol of its human ortholog when there is an unambiguous one: the gene is the reciprocal best hit of a human protein, that protein is its best human hit, the second-best scores below `--rbh-margin` of it, and the alignment covers the protein or the gene. An insect gene can then carry Dad (FlyBase) and SMAD6 (human). `--secondary-reference-taxid` chooses another species or disables it.

**Paralog designators.** Other genomes' nomenclatures add designators to the names of their own duplicates: teleost `a`/`b` ("Forkhead box protein O3-B", "protein kinase C beta a"), ", tandem duplicate 1". These mean nothing in the genome being annotated. ", tandem duplicate N" is stripped from hit names like any qualifier. A trailing ` a`/` b`, a letter glued to a trailing number ("O3b") or `-A`/`-B` after a number ("O3-B") is removed from a voted product name only when the name without it is a reference protein name or a name carried by the gene's own hits, so that "Cytochrome b" and "Annexin A" are left alone (`product_note` = `designator_stripped` in the evidence table; `--keep-paralog-designators` disables this).

Several genes may carry the same symbol. Paralogs, recent duplicates and retrocopies are real, and each copy is named for the gene it is homologous to; GOAnnotate does not reserve a symbol for one gene or rename copies with a "-like" suffix. Nucleotide (blastx) input skips this pass, since the reverse search needs protein queries; `--no-reciprocal` disables it.

### Phase 9: Confidence

Each gene receives a confidence score in [0, 1] and a label:

```
score = 0.45 * name_concordance
      + 0.25 * min(1, winning_cluster_size / 10)
      + 0.15 * (symbol assigned)
      + 0.15 * symbol_concordance
```

The score is lowered by 0.1 when no hit of the gene passed the hit coverage filter (coverage status `fallback`). Labels: `high` (>= 0.7), `medium` (>= 0.4), `low` (< 0.4), `none` (no hits). The label is written to the evidence table and to the annotated GFF3 (`goannotate_confidence`).

### Phase 10: Comparison with Existing Annotation

When the GFF3 already carries gene symbols (`gene=` by default) or product names (`product=` or `description=`), GOAnnotate compares them with its own result for every gene. Symbols are compared case-insensitively. Products are compared after normalization (case, punctuation, and qualifiers such as "transcript variant X2", "isoform", "-like", "precursor" are removed) with the same similarity measure used for clustering: `match` (>= 0.8), `partial` (>= 0.5), `mismatch`. Genes are also classified as `goannotate_only`, `existing_only` or `neither`.

The comparison is stratified by any gene- or transcript-level attributes named with `--strata-gene` and `--strata-transcript` (defaults: GRAFT's `graft_status` and `graft_protein`); attributes absent from the file are skipped. Transcript-level values are collapsed per gene, with several distinct values reported as `mixed(a|b)`.

Resolution follows from the annotation mode. In the default `fill` mode the existing annotation is kept wherever it exists, genes without an existing symbol or product (LOC-style genes, novel genes) receive the GOAnnotate values, and every gene where the two symbols both exist and disagree is listed in `{prefix}_review.tsv` with both annotations, the stratification attributes and the confidence, while the GFF3 keeps the existing symbol and records the GOAnnotate one as `goannotate_symbol`. `--annotation-mode replace` takes the GOAnnotate values everywhere instead.

### Concordance Metrics

For each annotated gene, GOAnnotate computes two concordance scores against the full merged hit pool (not just the winning cluster):

- **Name concordance:** The fraction of all merged hits whose cleaned protein name matches the final assigned name (at or above the similarity threshold).
- **Symbol concordance:** The fraction of all merged hits whose gene symbol exactly matches the final assigned symbol (case-insensitive).

A name concordance of 1.0 means every hit for that gene agrees with the chosen annotation. Low concordance (e.g., < 0.2) indicates the gene's hits are split across divergent annotations, which can happen with multidomain proteins, gene fusions, or genuinely ambiguous homology. Symbol concordance is typically much lower than name concordance because orthologous genes across species often share protein names but carry different gene symbols.

## Output Files

The annotation files are written to the top level of the `--output` directory; reports and intermediate files go to subdirectories. `{prefix}` is `--prefix`, or the name of the output directory.

```
{output}/
    {prefix}_annotated.gff3                  input GFF3 with GOAnnotate results
    {prefix}_annotated_AA.faa                proteins with annotated headers
    {prefix}_annotated_CDS.fna               coding sequences with annotated headers (genome extraction)
    {prefix}_annotated_transcripts.fna       spliced transcripts with annotated headers (genome extraction,
                                             or nucleotide --transcripts input)
    {prefix}_gene_annotations.tsv            one row per gene: symbol, product, locus, ..., GO terms
    {prefix}_transcript_annotations.tsv      one row per transcript
    {prefix}_gene2symbol.tsv, _gene2name.tsv, _gene2go.tsv
    {prefix}_summary.txt
    reports/
        {prefix}_comparison.tsv, _comparison_summary.txt, _review.tsv   comparison with the input annotation
        {prefix}_annotation_evidence.tsv     per-gene evidence for every decision
        {prefix}_go_evidence.tsv             per gene and GO term: evidence codes, class, supporting hits;
                                             terms removed by the taxon constraints
        {prefix}_name_audit.txt              what the name and symbol filters removed and kept
        {prefix}_extraction_report.tsv       per-transcript extraction flags (genome mode)
        {prefix}_protein_input_check.tsv     problems found in a --proteins file
    intermediate/
        {prefix}_AA.faa, _CDS.fna, _transcripts.fna   extracted sequences (genome mode)
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

`symbol` and `product` follow `--annotation-mode`, as in the annotated GFF3 and FASTA: with `fill` (default) and `add` the input annotation is kept where it exists and GOAnnotate's result fills the rest; with `replace` GOAnnotate's result is used wherever there is one. `symbol_source` and `product_source` say which (`input` or `goannotate`). Transcript rows keep the input transcript's own product, including its isoform label, under `fill` and `add`; the gene row uses the gene description or the first coding transcript's product without the isoform label. Loci are 1-based and inclusive (`chr1:23246-24202(-)`); `protein_length` is the length of the searched protein. `uniprot_accession` is the best hit of the winning symbol group. Without a GFF3 (transcript FASTA input) the rows come from the FASTA, with input values read from `symbol=`, `product=`, `type=` and `locus=` header fields.

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
| name_concordance | Fraction of all hits matching the assigned name |
| symbol_concordance | Fraction of all hits matching the assigned symbol |
| go_source | How GO terms were collected (`cluster_top5`, `cluster_top5+swissprot_supplement`, or empty) |
| n_go_terms | Number of specific GO terms assigned |
| confidence | `high`, `medium`, `low` or `none` |
| confidence_score | Numeric confidence in [0, 1] |
| gene_type | Gene biotype or feature type from the GFF3 (genome mode) |
| coding | `yes` or `no` (genome mode) |
| symbol_status | Outcome of the reference checks: `rbh`, `rbh_override`, `rbh_added`, `vote_supported`, `reference_best_hit`, `unsupported_removed`, `foreign_ortholog`, `foreign_kept`, `foreign_replaced`, `foreign_removed`, `vote`, `vote_rbh_conflict`, `propagated`, `none` |
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
| `{prefix}_annotated_AA.faa` | The searched proteins (extracted, `--proteins`, or protein `--transcripts`) |
| `{prefix}_annotated_CDS.fna` | Extracted coding sequences (genome mode) |
| `{prefix}_annotated_transcripts.fna` | Extracted spliced transcripts, all types (genome mode), or the nucleotide `--transcripts` input |

The FASTA files carry GOAnnotate's symbol and product in the header format described under Input Modes; `--annotation-mode` applies as for the GFF3 (`add`: `goannotate_symbol=` and `goannotate_product=` appended; `fill`: `symbol=` and `product=` set where absent; `replace`: overwritten).

| File | Description |
|------|-------------|
| `{prefix}_annotated.gff3` | The input GFF3 with GOAnnotate results added (genome mode, or `--gff3` with `--transcripts`) |

The annotated GFF3 is the input file passed through line by line: every feature, attribute, comment and `###` separator is kept, including non-coding and pseudogene models. Gene features receive `goannotate_symbol`, `goannotate_product`, `goannotate_confidence`, `Ontology_term` (the specific GO terms) and a `Dbxref=UniProtKB:ACC` for the top hit. `--annotation-mode` controls the standard attributes:

| Mode | Behaviour |
|------|-----------|
| `add` | Only the attributes above are added |
| `fill` (default) | `gene`, `Name` and `description` on genes and `product` on transcripts and CDS are written where the input has none; with `--hypothetical`, coding models still without a product receive `product=hypothetical protein` (off by default, since table2asn assigns it itself). Existing symbols that are identifiers (LOC numbers of an NCBI annotation, locus tags) and placeholder products ("uncharacterized protein LOC122519524") count as missing and are filled; `--keep-existing-identifiers` keeps them |
| `replace` | The same attributes are overwritten with GOAnnotate values wherever GOAnnotate has one |

### Summary

`{prefix}_summary.txt` contains a text summary of the run: command, inputs, parameters, GFF3 and extraction statistics, Diamond statistics, annotation coverage, confidence, GO term statistics, concordance statistics, and the comparison summary.

## Optional Modules

Modules wrap external tools that add information to the annotation without changing the symbol, product or GO assignment made by the core pipeline. They are off by default and selected with `--modules NAME[,NAME]`, or `--modules auto` to run every module whose tool is found. The tools are installed separately by the user; `goannotate --check-tools` reports which are available. Modules run after the symbol checks and before the outputs, so that their reports can list the final symbol and product of every gene. Their raw outputs are kept in `intermediate/{prefix}_{module}/` and reused by later runs into the same directory (`--force-modules` reruns them). Modules need protein queries and are skipped for nucleotide (`blastx`) input.

Each module appends its columns to the end of the gene and transcript tables, adds `goannotate_*` attributes (or standard INSDC features) to the annotated GFF3, writes a report to `reports/`, and is summarized in a Modules section of the summary file.

### signalp: signal peptides and the secretome

Runs [SignalP 6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (Teufel et al. 2022, Nature Biotechnology 40:1023) on the query proteins in eukarya mode and records, for every protein, whether a Sec/SPI signal peptide is predicted, its probability and the cleavage site.

SignalP 6 is distributed under an academic licence: register at the DTU page, download the package (the `fast` model at least; the `slow` models are a separate download), install it with `pip install signalp-6-package/` and copy its model weights into the installed package as the SignalP installation notes describe. The package supports Python up to 3.10 and PyTorch below 2.0, so it will usually live in its own conda environment; point `--signalp-bin` at that environment's `signalp6` executable (for example `~/miniconda3/envs/signalp6/bin/signalp6`). The module is skipped with a note when the executable is not found.

Only the first 70 residues of a protein are used by SignalP, so the proteins are truncated to 200 residues before the run; stop codons and non-standard letters are replaced. `--signalp-mode` chooses the model: `fast` (default), `slow` (the full model, needs more than 14 GB RAM) or `slow-sequential` (the full model within the memory of `fast`, about six times slower). The run uses `--threads` for the model and up to 8 writer processes.

Outputs:

| Where | Content |
|-------|---------|
| gene table | `signal_peptide`: `yes`, `no`, or `yes (1/3 isoforms)` when only some coding isoforms carry one; empty for genes without a scored protein |
| transcript table | `signal_peptide` (`yes`/`no`), `sp_probability`, `sp_cleavage` (`22-23`: the signal peptidase cleaves between residues 22 and 23) |
| annotated GFF3 | `goannotate_signal_peptide=yes` on genes and transcripts with a signal peptide, `goannotate_sp_cleavage` on the transcript; with `--signalp-features`, a `sig_peptide` feature (INSDC type, accepted by table2asn) under each such mRNA, with genomic coordinates mapped through the CDS segments (one line per segment when the signal peptide spans an intron) |
| `reports/{prefix}_secretome.tsv` | gene_id, transcript_id, symbol, product, sp_probability, cleavage_site, cs_probability for every protein with a signal peptide |
| `intermediate/{prefix}_signalp/` | SignalP's own `prediction_results.txt`, `output.gff3`, `processed_entries.fasta` (mature proteins) and the run log |

No GO terms are added: a signal peptide is not a GO assertion by itself.

## The bad_names.txt File

The `bad_names.txt` file controls how protein names and gene symbols from database hits are cleaned and which are rejected. Its purpose is to prevent carryover of designations that are specific to another genome or annotation pipeline (LOC numbers, locus tags, RIKEN and other clone IDs, chromosome ORF names, the NCBI `-like` paralog suffix) and to keep uninformative names ("uncharacterized protein") from being applied as products.

The copy shipped inside the package (`goannotate/data/bad_names.txt`) is used by default. `--bad-names FILE` replaces it; `--bad-names-extra FILE` (repeatable) adds patterns on top of it, which is the recommended way to tune for a taxon, for example with the locus-tag style of the genomes that dominate your hits. Every run writes `reports/{prefix}_name_audit.txt`, which reports discards by reason, identifiers and qualifiers removed by pattern, the most frequent surviving names (family-level names marked `G`), surviving names and symbols that look like identifiers, and the family-level products assigned to genes. `tests/audit_bad_names.py` produces the same report for a Diamond results file without running the pipeline.

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
| `goannotate/` | The package: `cli.py`, `pipeline.py`, `gff3.py`, `sequences.py`, `search.py`, `blast.py`, `filters.py`, `annotate.py`, `assign.py`, `reference.py`, `go.py`, `taxon.py`, `mappings.py`, `compare.py`, `output.py`, `models.py` |
| `goannotate/modules/` | Optional modules: `base.py` (framework), `signalp.py` |
| `goannotate/data/bad_names.txt` | Default name, symbol and product patterns |
| `goannotate/data/stop_words.txt` | Default stop words for the name similarity |
| `build_databases.py` | Database download and build script |
| `clean_gff3.py` | GFF3 ID cleanup utility |
| `tests/test_annotation_rules.py` | Checks of the coverage filter, word weighting, GO evidence, taxon constraints, product rules and symbol policy on small inputs |
| `tests/test_modules.py` | Checks of the module framework and the SignalP module with a stand-in `signalp6` script |
| `tests/validate_extraction.py` | Extraction validation against an external protein set |

## License

GPL-3.0-or-later. See LICENSE.
