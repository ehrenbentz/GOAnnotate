"""Core data structures shared across GOAnnotate modules."""

import re
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class BlastStats:
    """Statistics from BLAST results parsing."""
    total_hits_in_file: int = 0
    hits_filtered_by_evalue: int = 0
    hits_after_evalue_filter: int = 0
    queries_in_file: int = 0
    queries_after_evalue_filter: int = 0
    swissprot_hits: int = 0
    trembl_hits: int = 0
    other_hits: int = 0
    hits_with_lengths: int = 0


@dataclass
class BlastHit:
    """A single BLAST/Diamond hit (format 6 with stitle)."""
    query_id: str
    subject_id: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    stitle: str = ""
    qlen: int = 0      # query length (residues for blastp, nucleotides for blastx); 0 when not recorded
    slen: int = 0      # subject length (residues); 0 when not recorded

    @property
    def has_lengths(self) -> bool:
        """True when the search recorded query and subject lengths (qlen, slen)."""
        return self.qlen > 0 and self.slen > 0

    @property
    def query_coverage(self) -> float:
        """Fraction of the query spanned by the alignment, or 0.0 when qlen is unknown."""
        if self.qlen <= 0:
            return 0.0
        return min(1.0, (abs(self.qend - self.qstart) + 1) / self.qlen)

    @property
    def subject_coverage(self) -> float:
        """Fraction of the subject spanned by the alignment, or 0.0 when slen is unknown."""
        if self.slen <= 0:
            return 0.0
        return min(1.0, (abs(self.send - self.sstart) + 1) / self.slen)

    @property
    def gene_name(self) -> str:
        """Extract the protein description from the subject title."""
        if not self.stitle:
            return self.subject_id

        # Common patterns in BLAST titles:
        # >sp|P12345|GENE_SPECIES Description OS=Species ...
        # >tr|A0A123|A0A123_SPECIES Description OS=Species ...
        # >gi|123456|ref|NP_001234.1| Description [Species]
        # >ENSP00000123456 Description
        title = self.stitle

        # Remove OS= and everything after (UniProt format)
        if ' OS=' in title:
            title = title.split(' OS=')[0]

        # Remove [Species] suffix
        title = re.sub(r'\s*\[.*?\]\s*$', '', title)

        # Remove GN= gene name tag if present (we want the description)
        title = re.sub(r'\s*GN=\S+', '', title)

        # Remove PE= and SV= tags
        title = re.sub(r'\s*PE=\d+', '', title)
        title = re.sub(r'\s*SV=\d+', '', title)

        # Remove ECO evidence codes from TrEMBL entries {ECO:...}
        title = re.sub(r'\s*\{ECO:[^}]+\}', '', title)

        return title.strip()

    @property
    def gene_symbol(self) -> str:
        """Extract the gene symbol from the GN= field in a UniProt title."""
        if not self.stitle:
            return ""
        match = re.search(r'GN=(\S+)', self.stitle)
        if match:
            return match.group(1)
        return ""

    @property
    def taxid(self) -> str:
        """NCBI taxonomy ID from the OX= field in a UniProt title, or empty."""
        if not self.stitle:
            return ""
        match = re.search(r'OX=(\d+)', self.stitle)
        return match.group(1) if match else ""

    @property
    def accession(self) -> str:
        """UniProt accession from sp|ACC|ENTRY or tr|ACC|ENTRY, else the subject ID."""
        parts = self.subject_id.split('|')
        if len(parts) >= 2 and parts[0] in ('sp', 'tr'):
            return parts[1]
        return self.subject_id

    @property
    def is_swissprot(self) -> bool:
        """True for Swiss-Prot (reviewed) subjects."""
        return self.subject_id.startswith('sp|')

    @property
    def is_trembl(self) -> bool:
        """True for TrEMBL (unreviewed) subjects."""
        return self.subject_id.startswith('tr|')


@dataclass
class QueryAnnotation:
    """Annotation results for a single gene."""
    query_id: str
    hits: list = field(default_factory=list)
    filtered_hits: list = field(default_factory=list)
    consensus_name: str = ""
    consensus_symbol: str = ""
    go_terms: set = field(default_factory=set)
    specific_go_terms: set = field(default_factory=set)
    name_concordance: float = 0.0     # fraction of the gene's hits whose name matches the final product
    symbol_concordance: float = 0.0   # fraction of the gene's hits with a symbol that name the final symbol's gene
    hit_names: dict = field(default_factory=dict)    # normalized hit name -> number of hits (for concordance)
    hit_symbols: dict = field(default_factory=dict)  # upper-case hit symbol -> number of hits (for concordance)
    hit_curated_names: dict = field(default_factory=dict)  # normalized names of SwissProt hits of other species
                                                           # than the reference -> their upper-case symbols
    ortholog_reference_symbol: str = ""  # close relatives' symbol in the reference nomenclature (orthology module)
    ortholog_evidence: str = ""  # orthology module's evidence for the symbol (annotate.ORTHOLOG_CONFIDENCE key)
    n_isoforms: int = 0
    n_cluster_hits: int = 0
    winning_cluster_names: list = field(default_factory=list)
    winning_cluster_symbols: list = field(default_factory=list)
    mean_cluster_bitscore: float = 0.0
    go_source: str = ""
    confidence: str = ""
    confidence_score: float = 0.0
    symbol_status: str = ""     # rbh, rbh_override, rbh_added, vote_supported, reference_best_hit,
                                # vote (no reference hit), vote_rbh_conflict, propagated,
                                # foreign_ortholog, foreign_kept, foreign_replaced, foreign_removed, none;
                                # orthology module: ortholog (close relatives' symbol taken), ortholog_fill
    voted_symbol: str = ""      # the voted symbol, when the reference check replaced it
    coverage_status: str = ""   # "" (filter off or no hits), pass (hits passed the coverage filter),
                                # fallback (no hit passed; the unfiltered hits were used)
    best_query_coverage: float = 0.0    # best query coverage among the hits of the winning symbol group
    best_subject_coverage: float = 0.0  # best subject coverage among the hits of the winning symbol group
    go_evidence: dict = field(default_factory=dict)   # GO ID -> evidence codes (comma-separated) behind it
    go_support: dict = field(default_factory=dict)    # GO ID -> number of source hits carrying it
    go_taxon_removed: set = field(default_factory=set)  # terms removed by the taxon constraints
    go_n_source: int = 0        # hits the GO terms were collected from (top hits of the symbol group plus
                                # the SwissProt supplement)
    go_need: int = 0             # supporting hits a term from the hits alone needs in this gene (0: --go-min-hits)
    go_single_strong: bool = False  # the only source hit is a near full-length, high-identity alignment
    go_source_accs: set = field(default_factory=set)   # UniProt accessions of the GO source hits
    go_ortholog_support: dict = field(default_factory=dict)  # GO ID -> close relatives whose ortholog carries it
    go_curated: set = field(default_factory=set)   # terms of SwissProt entries of the gene with curated evidence
    go_fallback: dict = field(default_factory=dict)  # weak term -> {its most specific passing ancestors: support}
    go_weak_removed: dict = field(default_factory=dict)  # GO ID -> (support, codes, source) of removed terms
                                                         # alone removed for weak support
    product_note: str = ""      # designator_stripped when a paralog designator was removed from the product
    secondary_symbol: str = ""  # symbol in the secondary nomenclature (unambiguous reciprocal best hit)
    symbol_aliases: list = field(default_factory=list)  # other symbols of the same gene (hits' GN= symbols,
                                                        # the replaced voted symbol, reference-species synonyms)
    go_interpro: dict = field(default_factory=dict)     # GO ID -> InterPro entries that carry it (InterProScan module)
    go_reference: dict = field(default_factory=dict)    # GO ID -> reference proteins of the gene's symbol carrying it


@dataclass
class GeneUniverse:
    """
    The complete set of genes that appear in the output, with the mapping
    between transcripts and genes.

    Built either from FASTA headers (transcript ID first, gene= field) or
    from a parsed GFF3 file. Every gene appears in the output whether or
    not it has BLAST hits.
    """
    gene_ids: list                 # unique gene IDs in first-seen order
    transcript_to_gene: dict       # transcript_id -> gene_id
    gene_to_transcripts: dict      # gene_id -> list of transcript_ids
    total_sequences: int = 0
    gene_types: dict = field(default_factory=dict)   # gene_id -> biotype or feature type
    coding: dict = field(default_factory=dict)       # gene_id -> True if any transcript has a CDS

    def is_coding(self, gene_id: str) -> Optional[bool]:
        """Coding status of a gene, or None when it is not known (FASTA input)."""
        return self.coding.get(gene_id)
