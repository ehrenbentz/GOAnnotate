# GOAnnotate

Gene Ontology Annotation Pipeline using Diamond BLASTx and UniProt

## Description

GOAnnotate is a Python-based pipeline for annotating genes with Gene Ontology (GO) terms using Diamond BLAST searches against UniProt databases. The tool filters uninformative annotations based on customizable bad name patterns and provides comprehensive GO term assignment with consensus-based filtering.

## Author

E. J. Bentz (2025)

## Files

- `GOAnnotate.py` - Main annotation script
- `create_go_mapping.py` - GO mapping file creation from UniProt GAF files
- `dat_to_fasta.py` - Convert UniProt DAT format to FASTA
- `bad_names.txt` - Pattern file for filtering uninformative gene/protein names
- `GO_annotation_pipeline.txt` - Complete pipeline documentation and workflow

## Dependencies

- Python 3
- Diamond BLASTx
- AGAT

## Pipeline Overview

### 1. Prepare GO Mapping File

Download GO annotations from UniProt and create mapping file:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
./create_go_mapping.py goa-gaf --gaf goa_uniprot_all.gaf.gz --output Uniprot_All_GO_mapping.tsv
```

### 2. Prepare UniProtKB Database

Download and build Diamond databases for Swiss-Prot and TrEMBL:
```bash
# Download databases
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

# Build Diamond databases
diamond makedb --in uniprot_sprot.fasta.gz -d SwissProt_Diamond_DB
diamond makedb --in uniprot_trembl.fasta.gz -d TrEMBL_Diamond_DB
```

### 3. Extract CDS from Gene Models

Extract coding sequences from genome annotations using AGAT:
```bash
# Extract CDS sequences
agat_sp_extract_sequences.pl -g Annotations.gff3 -f Genome.fasta -t cds --mrna -o temp_CDS.fasta

# Reformat headers with gene ID as first field
awk '/^>/ {
    for(i=2; i<=NF; i++) {
        if($i ~ /^gene=/) {
            gene_id = substr($i, 6);
            break;
        }
    }
    printf ">%s %s", gene_id, substr($1,2);
    for(i=2; i<=NF; i++) printf " %s", $i;
    print "";
    next
} {print}' temp_CDS.fasta > CDS.fasta
```

### 4. Run Diamond BLASTx

Search against Swiss-Prot and TrEMBL separately, then combine results:
```bash
diamond blastx \
    -d SwissProt_Diamond_DB.dmnd \
    -q CDS.fasta \
    -o SwissProt_diamond_blastx_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 \
    -e 1e-5 \
    --threads 60

diamond blastx \
    -d TrEMBL_Diamond_DB.dmnd \
    -q CDS.fasta \
    -o TrEMBL_diamond_blastx_results.tsv \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    -k 50 \
    -e 1e-5 \
    --threads 60

# Combine results
cat SwissProt_diamond_blastx_results.tsv TrEMBL_diamond_blastx_results.tsv > UniProt_diamond_blastx_results.tsv
```

### 5. Download GO Ontology
```bash
wget https://purl.obolibrary.org/obo/go.obo -O go.obo
```

### 6. Run GOAnnotate
```bash
./GOAnnotate.py \
    --gff Annotations.gff3 \
    --transcripts CDS.fasta \
    --blast UniProt_diamond_blastx_results.tsv \
    --go-mapping Uniprot_All_GO_mapping.tsv \
    --ncbi-idmapping idmapping_selected.tab \
    --ncbi-geneinfo gene_info \
    --go-obo go.obo \
    --bad-names bad_names.txt \
    --evalue 1e-5 \
    --top-n 30 \
    --consensus-threshold 0.5 \
    --min-consensus-fraction 0.4 \
    --output Annotations
```

## Key Parameters

- `--evalue`: E-value threshold for BLAST hits (default: 1e-5)
- `--top-n`: Number of top BLAST hits to consider per query (default: 30)
- `--consensus-threshold`: Fraction of hits required to assign GO term (default: 0.5)
- `--min-consensus-fraction`: Minimum fraction for consensus filtering (default: 0.4)

## Notes

- Gene IDs must be the first field in FASTA headers for proper aggregation
- AGAT is strongly recommended over gffread for sequence extraction
- Genome FASTA should be line-wrapped to 80 or 60 characters for AGAT compatibility
- Swiss-Prot and TrEMBL searches are performed separately to prioritize curated annotations

## License

[To be added]
