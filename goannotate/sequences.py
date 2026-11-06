"""Genome access, sequence extraction, translation, and FASTA handling."""

import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from goannotate import GOAnnotateError
from goannotate.models import GeneUniverse

logger = logging.getLogger(__name__)


####################
# NUCLEOTIDES
####################

_COMPLEMENT = str.maketrans("ACGTUNRYKMSWBDHVacgtunrykmswbdhv",
                            "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def load_genome(path: str) -> dict:
    """
    Read a genome FASTA into memory as {sequence_name: uppercase sequence}.
    The sequence name is the first whitespace-delimited token of the header.
    """
    logger.info(f"Loading genome from {path}")
    genome = {}
    name = None
    chunks = []
    total = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(chunks).upper()
                    genome[name] = seq
                    total += len(seq)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seq = "".join(chunks).upper()
        genome[name] = seq
        total += len(seq)
    logger.info(f"  {len(genome):,} sequences, {total:,} bp")
    return genome


def spliced_sequence(seq: str, intervals: list, strand: str) -> str:
    """
    Concatenate 1-based inclusive genomic intervals (ascending) and
    reverse-complement for the minus strand. The result is in 5' to 3'
    transcript orientation.
    """
    parts = [seq[s - 1:e] for s, e in intervals]
    joined = "".join(parts)
    if strand == "-":
        return reverse_complement(joined)
    return joined


####################
# TRANSLATION
####################

_BASES = "TCAG"

# NCBI genetic code tables as amino acid strings in TCAG order
_TABLE_STRINGS = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}


def _build_table(aa_string: str) -> dict:
    table = {}
    i = 0
    for b1 in _BASES:
        for b2 in _BASES:
            for b3 in _BASES:
                table[b1 + b2 + b3] = aa_string[i]
                i += 1
    return table


_TABLES = {tid: _build_table(s) for tid, s in _TABLE_STRINGS.items()}


def codon_table(table_id: int) -> dict:
    if table_id not in _TABLES:
        raise GOAnnotateError(f"Genetic code table {table_id} is not available "
                              f"(available: {sorted(_TABLES)})")
    return _TABLES[table_id]


def stop_codons(table_id: int) -> set:
    return {c for c, aa in codon_table(table_id).items() if aa == "*"}


def translate(seq: str, table_id: int = 1) -> str:
    """Translate in frame 0; an incomplete trailing codon is ignored. Ambiguous codons give X."""
    table = codon_table(table_id)
    seq = seq.upper()
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(table.get(seq[i:i + 3], "X"))
    return "".join(out)


_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "TERM": "*", "OTHER": "X", "Xaa": "X",
}

_TRANSL_EXCEPT_RE = re.compile(r"\(pos:\s*(complement\()?(\d+)(?:\.\.(\d+))?\)?\s*,\s*aa:\s*([A-Za-z]+)\s*\)")


def parse_transl_except(value: str) -> list:
    """
    Parse an NCBI transl_except attribute, for example
    ``(pos:1001..1003,aa:Sec)`` or ``(pos:complement(500..502),aa:TERM)``.
    Returns a list of (start, end, one_letter_aa).
    """
    out = []
    for m in _TRANSL_EXCEPT_RE.finditer(value):
        start = int(m.group(2))
        end = int(m.group(3)) if m.group(3) else start
        aa = _AA3.get(m.group(4), _AA3.get(m.group(4).capitalize(), "X"))
        out.append((start, end, aa))
    return out


def apply_transl_except(protein: str, cds_segments: list, strand: str, phase: int,
                        exceptions: list) -> tuple:
    """
    Replace residues named by transl_except. Positions are genomic and
    are mapped through the CDS segments to codon indices.

    Returns (protein, n_applied).
    """
    positions = []
    ordered = sorted(cds_segments, key=lambda c: c[0], reverse=(strand == "-"))
    for s, e, _ in ordered:
        if strand == "-":
            positions.extend(range(e, s - 1, -1))
        else:
            positions.extend(range(s, e + 1))
    index_of = {pos: i for i, pos in enumerate(positions)}
    applied = 0
    protein_list = list(protein)
    for start, end, aa in exceptions:
        first = end if strand == "-" else start
        if first not in index_of:
            continue
        tx_index = index_of[first] - phase
        if tx_index < 0:
            continue
        codon_index = tx_index // 3
        if 0 <= codon_index < len(protein_list):
            protein_list[codon_index] = aa
            applied += 1
    return "".join(protein_list), applied


####################
# EXTRACTION
####################

@dataclass
class ExtractionResult:
    protein_fasta: Path
    cds_fasta: Path
    transcript_fasta: Path
    report: Path
    n_transcripts: int = 0
    n_coding: int = 0
    n_proteins_written: int = 0
    flag_counts: Counter = field(default_factory=Counter)
    missing_sequences: set = field(default_factory=set)


def _write_fasta_record(fh, header: str, seq: str, width: int = 80) -> None:
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


####################
# FASTA HEADERS
####################
# Same format as GRAFT's final output (graft/final_output.py):
#   >ID gene=GENE symbol=SYMBOL product="PRODUCT" type=TYPE locus=CHROM:START-END(STRAND)
# Fields without a value are left out. The locus is the transcript span in
# transcript files and the CDS span in CDS and protein files (1-based,
# inclusive).

HEADER_FIELD_ORDER = ("gene", "symbol", "product", "type", "locus")
_HEADER_FIELD = re.compile(r'(\w+)=(?:"([^"]*)"|(\S+))')


def gene_symbol_of(gene) -> str:
    """The gene's symbol: the gene attribute, else Name when it is not just the ID or locus tag."""
    sym = gene.attrs.get("gene", "")
    if not sym:
        name = gene.attrs.get("Name", "")
        if name and name not in (gene.attrs.get("locus_tag", ""), gene.id):
            sym = name
    return sym


def transcript_type_of(tx) -> str:
    """Transcript feature type as GRAFT reports it: coding transcripts of other types are mRNA."""
    t = tx.ftype or "mRNA"
    if tx.coding and t not in ("mRNA", "transcript") and not t.endswith("_gene_segment"):
        t = "mRNA"
    return t


def format_header(seq_id: str, fields: dict) -> str:
    """Header text from an ID and fields; known fields in GRAFT order, others after; empty values dropped."""
    parts = [seq_id]
    keys = [k for k in HEADER_FIELD_ORDER if k in fields] + [k for k in fields if k not in HEADER_FIELD_ORDER]
    for k in keys:
        v = fields[k]
        if v is None or v == "":
            continue
        if k == "product" or (" " in str(v)):
            parts.append(f'{k}="' + str(v).replace('"', "'") + '"')
        else:
            parts.append(f"{k}={v}")
    return " ".join(parts)


def parse_header(header: str) -> tuple:
    """(ID, ordered dict of key=value fields) from a FASTA header (with or without '>'); quoted values allowed."""
    text = header[1:] if header.startswith(">") else header
    text = text.strip()
    if not text:
        return "", {}
    parts = text.split(None, 1)
    seq_id = parts[0]
    fields = {}
    if len(parts) > 1:
        for m in _HEADER_FIELD.finditer(parts[1]):
            fields[m.group(1)] = m.group(2) if m.group(2) is not None else m.group(3)
    return seq_id, fields


def sequence_header(gene, tx, kind: str) -> str:
    """Header for a transcript ('transcript') or its CDS / protein ('cds')."""
    if kind == "cds" and tx.cds:
        product = tx.cds_attrs.get("product") or tx.attrs.get("product", "")
        start = min(s for s, _e, _p in tx.cds)
        end = max(e for _s, e, _p in tx.cds)
    else:
        product = tx.attrs.get("product", "")
        if tx.exons:
            start = min(s for s, _e in tx.exons)
            end = max(e for _s, e in tx.exons)
        else:
            start, end = tx.start, tx.end
    fields = {
        "gene": gene.key,
        "symbol": gene_symbol_of(gene),
        "product": product,
        "type": transcript_type_of(tx),
        "locus": f"{tx.seqid}:{start}-{end}({tx.strand})",
    }
    return format_header(tx.key, fields)


def _codon_table_for(gene, tx, default_table: int, mito_seqids: set) -> int:
    for attrs in (tx.cds_attrs, tx.attrs, gene.attrs):
        value = attrs.get("transl_table")
        if value and value.isdigit():
            return int(value)
    if tx.seqid in mito_seqids:
        return 2
    return default_table


def extract_sequences(genes: list, genome: dict, output_dir: Path, prefix: str,
                      default_table: int = 1, mito_seqids: list = None,
                      report_dir: Path = None) -> ExtractionResult:
    """
    Write three FASTA files from the parsed GFF3 and the genome:

      {prefix}_AA.faa           translated CDS, one record per coding transcript
      {prefix}_CDS.fna          CDS nucleotide sequence per coding transcript
      {prefix}_transcripts.fna  spliced exon sequence per transcript (all types)

    The extraction report goes to report_dir (default: output_dir).

    Headers follow GRAFT: ``>transcript gene=GENE symbol=SYMBOL
    product="PRODUCT" type=TYPE locus=CHROM:START-END(STRAND)``, fields
    without a value left out; the locus is the CDS span in the CDS and
    protein files. A report TSV lists every
    transcript with its coding status, lengths, codon table and flags.
    Nothing is dropped silently: problems are flagged and counted.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    mito = set(mito_seqids or [])

    result = ExtractionResult(
        protein_fasta=output_dir / f"{prefix}_AA.faa",
        cds_fasta=output_dir / f"{prefix}_CDS.fna",
        transcript_fasta=output_dir / f"{prefix}_transcripts.fna",
        report=Path(report_dir or output_dir) / f"{prefix}_extraction_report.tsv",
    )

    logger.info("Extracting transcript, CDS and protein sequences")

    with open(result.protein_fasta, "w") as f_prot, \
            open(result.cds_fasta, "w") as f_cds, \
            open(result.transcript_fasta, "w") as f_tx, \
            open(result.report, "w") as f_rep:

        f_rep.write("transcript_id\tgene_id\tgene_type\ttranscript_type\tcoding\tseqid\tstrand\t"
                    "n_exons\ttranscript_length\tcds_length\tprotein_length\tcodon_table\tflags\n")

        for gene in genes:
            for tx in gene.transcripts:
                result.n_transcripts += 1
                flags = []
                seq = genome.get(tx.seqid)
                if seq is None:
                    result.missing_sequences.add(tx.seqid)
                    flags.append("sequence_missing_from_genome")
                    result.flag_counts["sequence_missing_from_genome"] += 1
                    f_rep.write(f"{tx.key}\t{gene.key}\t{gene.biotype}\t{tx.ftype}\t"
                                f"{'yes' if tx.coding else 'no'}\t{tx.seqid}\t{tx.strand}\t"
                                f"{len(tx.exons)}\t{tx.transcript_length}\t{tx.cds_length}\t0\t\t"
                                f"{';'.join(flags)}\n")
                    continue

                tx_seq = spliced_sequence(seq, tx.exons, tx.strand)
                _write_fasta_record(f_tx, sequence_header(gene, tx, "transcript"), tx_seq)
                header = sequence_header(gene, tx, "cds") if tx.coding else ""

                protein_length = 0
                table_id = ""
                if tx.coding:
                    result.n_coding += 1
                    table_id = _codon_table_for(gene, tx, default_table, mito)
                    cds_ranges = [(s, e) for s, e, _ in tx.cds]
                    cds_seq = spliced_sequence(seq, cds_ranges, tx.strand)
                    _write_fasta_record(f_cds, header, cds_seq)

                    phase = tx.first_phase
                    if phase:
                        flags.append(f"phase_trimmed:{phase}")
                    coding_seq = cds_seq[phase:]
                    if len(coding_seq) % 3 != 0:
                        flags.append("length_not_multiple_of_3")

                    protein = translate(coding_seq, table_id)

                    exc_value = tx.cds_attrs.get("transl_except") or tx.attrs.get("transl_except")
                    if exc_value:
                        exceptions = parse_transl_except(exc_value)
                        protein, n_applied = apply_transl_except(
                            protein, tx.cds, tx.strand, phase, exceptions)
                        flags.append(f"transl_except:{n_applied}/{len(exceptions)}")

                    if protein.endswith("*"):
                        protein = protein[:-1]
                    else:
                        flags.append("no_stop_codon")
                    if coding_seq[:3] != "ATG":
                        flags.append(f"start_codon:{coding_seq[:3]}")
                    internal = protein.count("*")
                    if internal:
                        flags.append(f"internal_stop:{internal}")
                    if table_id != 1:
                        flags.append(f"codon_table:{table_id}")
                    for attrs in (tx.attrs, tx.cds_attrs, gene.attrs):
                        if any(k in attrs for k in ("partial", "start_range", "end_range")):
                            flags.append("partial")
                            break
                    for attrs in (tx.attrs, tx.cds_attrs, gene.attrs):
                        if "pseudo" in attrs:
                            flags.append("pseudo")
                            break
                    if "exception" in tx.attrs or "exception" in tx.cds_attrs:
                        flags.append("exception")

                    protein_length = len(protein)
                    if protein_length == 0:
                        flags.append("empty_protein")
                    else:
                        _write_fasta_record(f_prot, header, protein)
                        result.n_proteins_written += 1

                for flag in flags:
                    result.flag_counts[flag.split(":")[0]] += 1

                f_rep.write(f"{tx.key}\t{gene.key}\t{gene.biotype}\t{tx.ftype}\t"
                            f"{'yes' if tx.coding else 'no'}\t{tx.seqid}\t{tx.strand}\t"
                            f"{len(tx.exons)}\t{len(tx_seq)}\t{tx.cds_length}\t{protein_length}\t"
                            f"{table_id}\t{';'.join(flags)}\n")

    logger.info(f"  Transcripts: {result.n_transcripts:,}; coding: {result.n_coding:,}; "
                f"proteins written: {result.n_proteins_written:,}")
    if result.flag_counts:
        logger.info("  Flags: " + ", ".join(
            f"{k}={v:,}" for k, v in result.flag_counts.most_common()))
    if result.missing_sequences:
        logger.warning(f"  {len(result.missing_sequences)} sequence names in the GFF3 are absent "
                       f"from the genome: {', '.join(sorted(result.missing_sequences)[:10])}")
    logger.info(f"  Wrote {result.protein_fasta}, {result.cds_fasta}, {result.transcript_fasta}")
    return result


####################
# FASTA UNIVERSE
####################

# Trinity transcript identifiers: the isoform suffix follows the gene
# (TRINITY_DN1000_c0_g1_i1 and DN1000_c0_g1_i1 since Trinity 2.1,
# TR1000|c0_g1_i1 in Trinity 2.0, comp1000_c0_seq1 before 2.0)
_TRINITY_ISOFORM = (
    re.compile(r"^((?:TRINITY_)?DN\d+_c\d+_g\d+)_i\d+$"),
    re.compile(r"^(TR\d+\|c\d+_g\d+)_i\d+$"),
    re.compile(r"^(comp\d+_c\d+)_seq\d+$"),
)


def trinity_gene(transcript_id: str) -> str:
    """The Trinity gene of a transcript identifier, or '' when the identifier is not Trinity's."""
    for pattern in _TRINITY_ISOFORM:
        m = pattern.match(transcript_id)
        if m:
            return m.group(1)
    return ""


def universe_from_fasta(fasta_file: str, isoform_grouping: str = "trinity") -> GeneUniverse:
    """
    Build the gene universe from FASTA headers.

    The transcript ID is the first field of each header; the gene ID comes
    from a gene= field. Without a gene= field, a Trinity identifier gives
    its gene (isoform_grouping trinity; off keeps every transcript as its
    own gene); any other transcript ID is used as the gene ID (1:1 mapping).
    """
    logger.info(f"Parsing FASTA headers from {fasta_file}")

    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(list)
    gene_ids_ordered = []
    gene_ids_seen = set()
    total_sequences = 0
    n_trinity = 0

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                total_sequences += 1
                transcript_id, fields = parse_header(line.rstrip("\n"))
                gene_id = fields.get("gene")
                if not gene_id and isoform_grouping == "trinity":
                    gene_id = trinity_gene(transcript_id)
                    if gene_id:
                        n_trinity += 1
                gene_id = gene_id or transcript_id

                transcript_to_gene[transcript_id] = gene_id
                gene_to_transcripts[gene_id].append(transcript_id)
                if gene_id not in gene_ids_seen:
                    gene_ids_ordered.append(gene_id)
                    gene_ids_seen.add(gene_id)

    logger.info(f"Found {total_sequences:,} sequences, {len(gene_ids_ordered):,} unique genes, "
                f"{len(transcript_to_gene):,} transcripts")
    if n_trinity:
        logger.info(f"  {n_trinity:,} transcripts grouped into genes by their Trinity identifiers "
                    f"(--isoform-grouping off keeps each transcript as a gene)")
    n_multi = sum(1 for txs in gene_to_transcripts.values() if len(txs) > 1)
    if n_multi > 0:
        logger.info(f"  {n_multi:,} genes have multiple transcripts")

    return GeneUniverse(gene_ids=gene_ids_ordered, transcript_to_gene=transcript_to_gene,
                        gene_to_transcripts=dict(gene_to_transcripts),
                        total_sequences=total_sequences)


####################
# PRE-EXTRACTED PROTEINS
####################

def read_fasta(path: str) -> list:
    """List of (header line without '>', sequence) records."""
    records = []
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].rstrip("\n")
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


@dataclass
class ProteinInputCheck:
    path: Path
    report: Path
    n_records: int = 0
    n_matched: int = 0
    n_coding_without_protein: int = 0
    n_noncoding_with_protein: int = 0
    n_length_mismatch: int = 0
    n_compared: int = 0
    n_identical: int = 0


def check_protein_fasta(path: str, genes: list, output_dir: Path, prefix: str,
                        extracted_fasta: Path = None) -> ProteinInputCheck:
    """
    Check a pre-extracted protein FASTA against the parsed GFF3 before it is used.

    Fatal (GOAnnotateError): duplicate IDs; IDs that are not transcripts in
    the GFF3 (with --transcript-id-attr applied); a gene= header field that
    names a different gene than the transcript's parent in the GFF3.

    Reported: coding transcripts without a protein (their genes are
    annotated from their other isoforms, or not at all), proteins for
    transcripts without CDS, and proteins whose length differs from the CDS
    length by more than one residue. With extracted_fasta (GOAnnotate's own
    extraction from --genome), every protein is also compared with the
    extracted one.

    Writes {prefix}_protein_input_check.tsv with one row per problem.
    """
    output_dir = Path(output_dir)
    result = ProteinInputCheck(path=Path(path), report=output_dir / f"{prefix}_protein_input_check.tsv")
    logger.info(f"Checking protein FASTA {path} against the GFF3")

    tx_index = {}
    for gene in genes:
        for tx in gene.transcripts:
            tx_index[tx.key] = (gene, tx)

    records = read_fasta(path)
    result.n_records = len(records)
    if not records:
        raise GOAnnotateError(f"No sequences in {path}")

    problems = []
    seen = {}
    duplicates = []
    unknown = []
    gene_mismatch = []
    supplied = {}
    for header, seq in records:
        seq_id, fields = parse_header(header)
        if seq_id in seen:
            duplicates.append(seq_id)
            continue
        seen[seq_id] = True
        supplied[seq_id] = seq.rstrip("*")
        hit = tx_index.get(seq_id)
        if hit is None:
            unknown.append(seq_id)
            continue
        gene, tx = hit
        tagged_gene = fields.get("gene")
        if tagged_gene and tagged_gene != gene.key:
            gene_mismatch.append((seq_id, tagged_gene, gene.key))
            continue
        result.n_matched += 1
        if not tx.coding:
            result.n_noncoding_with_protein += 1
            problems.append((seq_id, gene.key, "protein_for_noncoding_transcript", ""))
            continue
        expected = (tx.cds_length - tx.first_phase) // 3
        length = len(supplied[seq_id])
        # The CDS may or may not include the stop codon. Proteins written
        # with internal stops removed (GRAFT) are shorter by their number.
        if not (expected - 1 <= length <= expected):
            result.n_length_mismatch += 1
            problems.append((seq_id, gene.key, "length_differs_from_cds",
                             f"protein {length} aa; CDS {tx.cds_length} nt ({expected} codons)"))

    hint = (" The FASTA IDs must be the transcript identifiers of the GFF3 "
            "(see --transcript-id-attr / --gene-id-attr).")

    def _fail(msg, items, suffix=hint):
        shown = ", ".join(str(x) for x in items[:10])
        raise GOAnnotateError(f"{msg}: {len(items):,} (first: {shown}).{suffix}")

    if duplicates:
        _fail(f"Duplicate sequence IDs in {path}", duplicates, "")
    if unknown:
        _fail(f"Sequences in {path} whose ID is not a transcript in the GFF3", unknown)
    if gene_mismatch:
        _fail(f"Sequences in {path} whose gene= field differs from the GFF3 parent gene",
              [f"{t} gene={a} (GFF3: {b})" for t, a, b in gene_mismatch])

    for gene in genes:
        for tx in gene.transcripts:
            if tx.coding and tx.key not in supplied:
                result.n_coding_without_protein += 1
                problems.append((tx.key, gene.key, "coding_transcript_without_protein", ""))

    if extracted_fasta is not None and Path(extracted_fasta).exists():
        for header, seq in read_fasta(str(extracted_fasta)):
            seq_id, _ = parse_header(header)
            if seq_id in supplied:
                result.n_compared += 1
                if supplied[seq_id] == seq.rstrip("*"):
                    result.n_identical += 1
                else:
                    gene_key = tx_index[seq_id][0].key if seq_id in tx_index else ""
                    problems.append((seq_id, gene_key, "differs_from_genome_extraction", ""))

    with open(result.report, "w") as fh:
        fh.write("transcript_id\tgene_id\tproblem\tdetail\n")
        for row in problems:
            fh.write("\t".join(row) + "\n")

    logger.info(f"  {result.n_records:,} proteins; all IDs match GFF3 transcripts and genes")
    if result.n_coding_without_protein:
        logger.warning(f"  {result.n_coding_without_protein:,} coding transcripts in the GFF3 have no protein "
                       f"in the FASTA")
    if result.n_noncoding_with_protein:
        logger.warning(f"  {result.n_noncoding_with_protein:,} proteins belong to transcripts without CDS "
                       f"in the GFF3")
    if result.n_length_mismatch:
        logger.warning(f"  {result.n_length_mismatch:,} proteins differ in length from their GFF3 CDS by more "
                       f"than one residue")
    if result.n_compared:
        logger.info(f"  Compared with the extraction from the genome: {result.n_identical:,} of "
                    f"{result.n_compared:,} identical")
    if problems:
        logger.info(f"  Details: {result.report}")
    return result


def fasta_is_protein(path: str, n_records: int = 200) -> bool:
    """
    Heuristic: True when the first records contain residues that are not
    IUPAC nucleotide codes (E, F, I, L, P, Q, X, Z, * or similar).
    """
    protein_only = set("EFILPQXZJO*efilpqxzjo")
    seen = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                seen += 1
                if seen > n_records:
                    break
                continue
            if any(ch in protein_only for ch in line.strip()):
                return True
    return False
