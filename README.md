# GOAnnotate

Functional annotation of predicted gene models using Diamond BLASTx against UniProt.

GOAnnotate assigns protein names, gene symbols, and Gene Ontology (GO) terms to genes by analyzing BLAST hits against SwissProt and TrEMBL.
Handles isoform aggregation, clusters hits by protein name similarity, and uses a voting/scoring system to pick the best annotation for each gene.
Designed for newly assembled and annotated genomes where gene models need functional descriptions.

**Author:** E. J. Bentz (2025)

## Contents

- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Database Setup](#database-setup)
- [Preparing Inputs](#preparing-inputs)
- [Running GOAnnotate](#running-goannotate)
- [How It Works](#how-it-works)
- [Output Files](#output-files)
- [The bad_names.txt File](#the-bad_namestxt-file)
- [Parameters](#parameters)
- [Utility Scripts](#utility-scripts)

## Requirements

- Python 3.8+
- [Diamond](https://github.com/bbuchfink/diamond) (for BLASTx searches)
- [AGAT](https://github.com/NBISweden/AGAT) (Suggested method for CDS extraction from gene models)

No Python packages beyond the standard library are required.

## Quick Start

```bash
# 1. Build databases (one-time setup, several hours)
./build_databases.py -o /path/to/databases --threads 8

# 2. Extract CDS from your genome annotation
# (Note: Any extraction method will work, but AGAT generally performs the best)
agat_sp_extract_sequences.pl -g Annotations.gff3 -f Genome.fasta -t cds --mrna -o CDS.fasta

# 3. Run Diamond BLASTx against SwissProt and TrEMBL separately
diamond blastx -d /path/to/databases/SwissProt_Diamond_DB.dmnd \
    -q CDS.fasta -o SwissProt_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 -e 1e-5 --threads 32

diamond blastx -d /path/to/databases/TrEMBL_Diamond_DB.dmnd \
    -q CDS.fasta -o TrEMBL_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 -e 1e-5 --threads 32

cat SwissProt_results.tsv TrEMBL_results.tsv > UniProt_results.tsv

# 4. Annotate
./GOAnnotate.py \
    --transcripts CDS.fasta \
    --blast-results UniProt_results.tsv \
    --db /path/to/databases/GOAnnotate_db.sqlite \
    --bad-names bad_names.txt \
    --gff Annotations.gff3 \
    --output my_annotations
```

## Database Setup

`build_databases.py` downloads and builds everything GOAnnotate needs from UniProt and NCBI. This is a one-time setup step. The databases can be reused across projects and should be updated periodically to stay current with UniProt releases.

```bash
./build_databases.py -o /path/to/databases --threads 32
```

**Be warned:** this downloads roughly 115 GB of compressed data and may temporarily need over 1 TB of disk space. Expect it to take up to 12-24 hours depending on your connection and hardware. The script will show a confirmation prompt before starting.

The script produces:

| File | Description |
|------|-------------|
| `SwissProt_Diamond_DB.dmnd` | Diamond database for SwissProt BLASTx |
| `TrEMBL_Diamond_DB.dmnd` | Diamond database for TrEMBL BLASTx |
| `GOAnnotate_db.sqlite` | SQLite database with GO mappings and NCBI cross-references |
| `Uniprot_All_GO_mapping.tsv` | Flat-file GO mapping (alternative to SQLite) |
| `idmapping_selected.tab` | UniProt ID mapping (alternative to SQLite) |
| `gene_info.tsv` | NCBI gene info (alternative to SQLite) |

The SQLite database (`--db`) is the recommended way to provide these to GOAnnotate. It consolidates the GO mapping, UniProt ID mapping, and NCBI gene info into a single file and loads significantly faster than the flat files. If you prefer flat files, you can use `--go-mapping`, `--ncbi-idmapping`, and `--ncbi-geneinfo` instead.

Options:

```
-o, --output-dir      Directory for database files (required)
--threads             Threads for diamond makedb (default: all CPUs)
--keep-intermediates  Keep compressed downloads after processing
--rebuild-index       Rebuild the SQLite database even if it exists
```

## Preparing Inputs

### Extract CDS Sequences

Use AGAT to extract coding sequences from your genome annotation. AGAT is strongly recommended for this step.

```bash
agat_sp_extract_sequences.pl -g Annotations.gff3 -f Genome.fasta -t cds --mrna -o CDS.fasta
```

**Note:** Your genome FASTA must be line-wrapped for AGAT compatibility.

GOAnnotate reads the FASTA headers to build a transcript-to-gene mapping. It looks for a `gene=` field anywhere in the header. If found, all transcripts sharing the same gene ID are grouped together and their BLAST hits are merged. If no `gene=` field is present, each sequence is treated as an independent unit. AGAT includes the `gene=` field automatically.

### Run Diamond BLASTx

Run Diamond BLASTx separately against SwissProt and TrEMBL, then concatenate the results. Running them separately is important because Diamond's `-k` flag limits results per database. If you searched a combined database, you would often get 50 TrEMBL hits and zero SwissProt hits for a given query, losing curated SwissProt annotations that are valuable for GO term assignment.

```bash
diamond blastx \
    -d SwissProt_Diamond_DB.dmnd \
    -q CDS.fasta \
    -o SwissProt_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 -e 1e-5 --threads 32

diamond blastx \
    -d TrEMBL_Diamond_DB.dmnd \
    -q CDS.fasta \
    -o TrEMBL_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 -e 1e-5 --threads 32

cat SwissProt_results.tsv TrEMBL_results.tsv > UniProt_results.tsv
```

The output format must be tab-separated format 6 with the 13 fields shown above. The `stitle` field (column 13) is required because GOAnnotate extracts protein names and gene symbols from it.

## Running GOAnnotate

### Minimal command (with SQLite database)

```bash
./GOAnnotate.py \
    --transcripts CDS.fasta \
    --blast-results UniProt_results.tsv \
    --db GOAnnotate_db.sqlite \
    --bad-names bad_names.txt \
    --output my_output
```

### Full command with all options

```bash
./GOAnnotate.py \
    --transcripts CDS.fasta \
    --blast-results UniProt_results.tsv \
    --db GOAnnotate_db.sqlite \
    --bad-names bad_names.txt \
    --gff Annotations.gff3 \
    --go-obo go.obo \
    --evalue 1e-5 \
    --top-n 50 \
    --consensus-threshold 0.5 \
    --namespace BP MF CC \
    --threads 48 \
    --output my_output \
    --prefix my_species
```

### Using flat files instead of SQLite

```bash
./GOAnnotate.py \
    --transcripts CDS.fasta \
    --blast-results UniProt_results.tsv \
    --go-mapping Uniprot_All_GO_mapping.tsv \
    --ncbi-idmapping idmapping_selected.tab \
    --ncbi-geneinfo gene_info.tsv \
    --bad-names bad_names.txt \
    --output my_output
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--transcripts` | CDS FASTA file |
| `--blast-results` | Combined Diamond BLASTx results (format 6 + stitle) |
| `--db` or `--go-mapping` | SQLite database or flat-file GO mapping (one required) |
| `--bad-names` | Pattern file for filtering uninformative names |
| `--output` | Output directory |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--gff` | none | GFF3 file (needed for annotated GFF3 output) |
| `--go-obo` | auto-download | GO hierarchy file (downloaded if missing or stale) |
| `--ncbi-idmapping` | none | UniProt-to-NCBI ID mapping (not needed with `--db`) |
| `--ncbi-geneinfo` | none | NCBI gene info (not needed with `--db`) |
| `--evalue` | 1e-5 | E-value threshold |
| `--top-n` | 30 | Top hits to keep per transcript after filtering |
| `--consensus-threshold` | 0.5 | Similarity threshold for name clustering |
| `--namespace` | BP MF CC | GO namespaces to include |
| `--threads` | all CPUs | Number of parallel workers |
| `--prefix` | output dir name | Prefix for output filenames |

## How It Works

GOAnnotate runs through a series of phases to go from raw BLAST hits to final annotations. Here is what happens at each step.

### Phase 0: Input Parsing

GOAnnotate reads the CDS FASTA file to build the complete gene universe. Every gene defined in the FASTA will appear in the output files, whether or not it has BLAST hits. This guarantees that downstream tools can use the output as a complete lookup table.

For each FASTA header, GOAnnotate extracts the transcript ID (first whitespace-delimited field) and scans remaining fields for a `gene=` tag. If found, that gene ID is used to group transcripts. If not found, the transcript ID is used as the gene ID (no grouping).

### Phase 1: BLAST Hit Parsing

BLAST results are parsed and filtered by e-value. Hits are keyed by transcript ID (matching the qseqid from Diamond). Each hit's subject ID is checked for `sp|` or `tr|` prefixes to identify SwissProt vs. TrEMBL origin. The protein name and gene symbol (from the `GN=` field in the stitle) are extracted from each hit.

### Phase 2: Hit Cleaning and Top-N Selection

Each hit's protein name goes through a two-stage cleaning process controlled by the `bad_names.txt` file:

1. **Artifact stripping** (`strip_regex:` patterns): database prefixes like `sp|Q12345|GENE_HUMAN`, prediction tags like `PREDICTED:`, isoform/variant designations, `-like` suffixes, and `homolog` are stripped from the name. The remaining text is the cleaned protein name.

2. **Bad name filtering** (plain text and `regex:` patterns): if the cleaned name matches any uninformative pattern (e.g., "hypothetical protein", "uncharacterized", LOC numbers, Ensembl IDs, gene prediction tool output), the hit is discarded.

Gene symbols go through a separate filter (`symbol_regex:` patterns) that removes uninformative symbols like LOC numbers, internal genome annotation IDs, and accession-like strings.

Top-N selection happens *after* bad name filtering. This is important: if filtering happened after truncation, uninformative hits would consume slots and push out real annotations. With the current order, you get up to N informative hits per transcript.

### Phase 3: Isoform Aggregation

Hits from all transcripts of the same gene are merged into a single pool. When the same subject (same UniProt accession) appears in multiple transcripts, the hit with the highest bitscore is kept.

### Phase 4: Agglomerative Clustering

The merged hits for each gene are clustered by protein name similarity using average-linkage agglomerative clustering. This groups hits that describe the same (or very similar) proteins together, even when the exact wording differs between species or databases.

Similarity between two hits is computed using `SequenceMatcher` on normalized protein names, with a critical modification for gene symbols:

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

This rewards clusters with strong BLAST evidence (high bitscores), broad support (many hits), and corroboration across isoforms. The square root on cluster size prevents a large but low-scoring cluster from beating a small cluster with very strong hits.

### Phase 6: Name and Symbol Selection

Within the winning cluster, hits are grouped by gene symbol. The symbol group with the highest mean bitscore is selected, and both the gene symbol and protein name come from that group.

A **single-hit guard** prevents annotation errors in large clusters: if a symbol group has only one hit and there are other groups with two or more hits, the single-hit group cannot win. This stops a lone misannotated hit from overriding the majority. An additional guard rejects a single-hit symbol that represents less than 10% of a cluster with 10+ hits.

If no valid gene symbol exists in the winning cluster, GOAnnotate falls back to NCBI cross-referencing: it looks up the best hit's UniProt accession in the NCBI ID mapping to find a GeneID, then looks up that GeneID in NCBI's gene_info to find a symbol. This step requires either the `--db` SQLite database or the `--ncbi-idmapping` and `--ncbi-geneinfo` flat files.

The protein name comes from the highest-bitscore hit in the winning symbol group.

### Phase 7: GO Term Collection

GO terms are collected in two steps:

1. **Primary collection:** The top 5 hits (by bitscore) from the winning symbol group are looked up in the GO mapping. Their GO terms are pooled together.

2. **SwissProt supplement:** The pipeline then scans all remaining hits in the winning cluster for SwissProt entries that share the winning gene symbol or protein name. GO terms from matching SwissProt hits are added to the pool. This is important because TrEMBL hits typically dominate the bitscore rankings (especially for non-model organisms with closely related species in TrEMBL), but SwissProt entries carry richer, manually curated GO annotations.

After collection, GO terms are filtered by namespace (BP, MF, CC) and then reduced to the most specific terms by removing redundant ancestor terms using the GO hierarchy from the OBO file. For example, if a gene has both GO:0005634 (nucleus) and GO:0005737 (cytoplasm), both are kept because neither is an ancestor of the other. But if it has both GO:0005634 (nucleus) and GO:0005622 (intracellular anatomical structure, a parent of nucleus), the parent is removed.

### Phase 8: Symbol Propagation

After all genes are annotated independently, a post-annotation pass looks for genes that share the same normalized protein name. If some of those genes have a gene symbol and others do not, the most common symbol is propagated to the symbol-less siblings. This fills in gaps where BLAST hits for a particular gene happened to lack the symbol even though sibling genes with identical annotations had it.

### Concordance Metrics

For each annotated gene, GOAnnotate computes two concordance scores against the full merged hit pool (not just the winning cluster):

- **Name concordance:** The fraction of all merged hits whose cleaned protein name matches the final assigned name (at or above the similarity threshold).
- **Symbol concordance:** The fraction of all merged hits whose gene symbol exactly matches the final assigned symbol (case-insensitive).

A name concordance of 1.0 means every BLAST hit for that gene agrees with the chosen annotation. Low concordance (e.g., < 0.2) indicates the gene's hits are split across divergent annotations, which can happen with multidomain proteins, gene fusions, or genuinely ambiguous homology. Symbol concordance is typically much lower than name concordance because orthologous genes across species often share protein names but carry different gene symbols.

## Output Files

All output files are written to the `--output` directory with the specified prefix.

### Simple Mapping Files

These are two-column tab-delimited files. Every gene in the input FASTA is included. Genes without annotations have an empty second column.

| File | Columns | Description |
|------|---------|-------------|
| `{prefix}_gene2symbol.tsv` | gene_id, gene_symbol | Gene ID to gene symbol mapping |
| `{prefix}_gene2name.tsv` | gene_id, protein_name | Gene ID to protein name mapping |
| `{prefix}_gene2go.tsv` | gene_id, go_terms | Gene ID to GO terms (semicolon-delimited) |

### Annotation Evidence File

`{prefix}_annotation_evidence.tsv` contains detailed per-gene evidence for every annotation decision. This is the file to look at if you want to understand why a particular gene received a particular annotation.

Columns:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| gene_symbol | Assigned gene symbol (empty if none) |
| protein_name | Assigned protein name (empty if none) |
| n_isoforms | Number of transcripts that contributed BLAST hits |
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

### Annotated Sequence Files

| File | Description |
|------|-------------|
| `{prefix}_annotated.fasta` | CDS FASTA with protein names and gene symbols in headers |
| `{prefix}_annotated.gff3` | GFF3 with `Name=` (gene symbol) and `product=` (protein name) attributes (only produced if `--gff` is provided) |

**Note on GFF3:** The annotated GFF3 replaces all existing column 9 attributes with GOAnnotate's annotations. If your input GFF3 already has functional annotations, they will be overwritten. GO terms are not included in the GFF3 output; use `gene2go.tsv` for those.

### Summary

`{prefix}_summary.txt` contains a text summary of the run: input files, parameters, BLAST statistics, annotation coverage, GO term statistics, and concordance statistics.

## The bad_names.txt File

The `bad_names.txt` file controls which protein names and gene symbols are considered uninformative and should be filtered out. A well-tuned bad names file is important for annotation quality, because uninformative hits (e.g., "hypothetical protein", "uncharacterized protein") can crowd out real annotations.

The included `bad_names.txt` covers common cases for vertebrate genomes, with extra patterns for fish, zebrafish clone IDs, and various genome annotation artifacts. You may need to add patterns for your specific organism or annotation source.

### Pattern types

**`strip_regex:` patterns** are applied first. They remove substrings from protein names without discarding the hit. Use these for database artifacts and qualifiers that obscure an otherwise informative name:

```
strip_regex:^PREDICTED:\s*
strip_regex:\(Fragment\)
strip_regex:\bisoform\s+\S+
strip_regex:-like\s*$
strip_regex:\bhomolog\b
```

**Plain text patterns** are matched as substrings (case-insensitive) against the cleaned protein name. Short patterns (6 characters or fewer) are matched exactly. If a hit's name matches, the hit is discarded:

```
ypothetical protein
ncharacterized
putative protein
```

(These partial strings intentionally match both capitalized and lowercase variants, e.g., "ypothetical" matches both "Hypothetical" and "hypothetical".)

**`regex:` patterns** are matched against cleaned protein names (case-insensitive):

```
regex:^LOC\d+
regex:^ENS[A-Z]+\d+$
regex:^maker-
```

**`symbol_regex:` patterns** are matched against gene symbols only:

```
symbol_regex:^LOC\d+$
symbol_regex:^[A-Z0-9]+_\d{5,}$
symbol_regex:^ENS[A-Z]+\d+$
```

Comments (lines starting with `#`) and blank lines are ignored.

## Parameters

### `--top-n` (default: 30)

The number of informative hits to keep per transcript after bad-name filtering. Higher values give the clustering step more evidence to work with but increase runtime. For well-annotated reference organisms, 30 is usually sufficient. For non-model organisms with sparser database representation, 50 may give better results.

### `--consensus-threshold` (default: 0.5)

The similarity threshold for clustering protein names. Two names must have at least this similarity score (from SequenceMatcher, 0 to 1) to be placed in the same cluster. Lower values produce larger, more inclusive clusters. Higher values produce smaller, more specific clusters. The default of 0.5 works well for most use cases.

### `--evalue` (default: 1e-5)

Standard e-value cutoff for BLAST hits. This is applied during parsing, before any other filtering.

### `--namespace` (default: BP MF CC)

Which GO namespaces to include: biological_process (BP), molecular_function (MF), cellular_component (CC). You can restrict to a subset if needed.

## Utility Scripts

### clean_gff3.py

A utility for stripping all existing attributes from a GFF3 file and re-assigning clean, sequentially numbered IDs. Useful for preparing messy GFF3 files before annotation. Parent references are updated to match the new IDs.

```bash
python3 clean_gff3.py input.gff3 output.gff3
```

This script is provided as a convenience utility and is not required by the main pipeline.

## Files

| File | Description |
|------|-------------|
| `GOAnnotate.py` | Main annotation pipeline |
| `build_databases.py` | Database download and build script |
| `bad_names.txt` | Default patterns for filtering uninformative names |
| `clean_gff3.py` | GFF3 ID cleanup utility |

## License

[To be added]
