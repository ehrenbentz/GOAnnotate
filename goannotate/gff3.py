"""GFF3 parsing into a gene > transcript > exon/CDS hierarchy.

The parser accepts gene and pseudogene parents (and any type ending in
``_gene``), every transcript type (mRNA, transcript, lnc_RNA, ncRNA,
miRNA, tRNA, V_gene_segment, ...), transcripts without a gene parent,
genes with exon or CDS children and no transcript, and CDS-only models.
The GFF3 phase of the first CDS segment in translation order is kept so
that 5' partial models translate correctly. Features nested under a
transcript (a miRNA under a primary_transcript, for example) are counted
and otherwise ignored. Parsing stops at a ``##FASTA`` directive.

Coordinates are kept as in GFF3: 1-based, inclusive.
"""

import logging
import urllib.parse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional

from goannotate import GOAnnotateError
from goannotate.models import GeneUniverse

logger = logging.getLogger(__name__)

GENE_TYPES = {"gene", "pseudogene"}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "lncRNA", "primary_transcript",
                    "pseudogenic_transcript", "unconfirmed_transcript", "processed_transcript",
                    "V_gene_segment", "C_gene_segment", "D_gene_segment", "J_gene_segment"}
# Features derived from exons and CDS; never modelled
DERIVED_TYPES = {"five_prime_UTR", "three_prime_UTR", "UTR", "five_prime_utr", "three_prime_utr",
                 "start_codon", "stop_codon", "intron", "polyA_sequence", "polyA_site"}
# Sequence-level features that never belong to a gene model
IGNORED_TYPES = {"region", "chromosome", "contig", "scaffold", "supercontig", "match",
                 "cDNA_match", "match_part", "repeat_region", "sequence_feature",
                 "biological_region", "centromere", "telomere", "gap"}


def is_gene_type(ftype: str) -> bool:
    return ftype in GENE_TYPES or ftype.endswith("_gene")


def is_transcript_type(ftype: str) -> bool:
    return (ftype in TRANSCRIPT_TYPES or ftype.endswith("RNA")
            or ftype.endswith("_transcript") or ftype.endswith("_gene_segment"))


def parse_attributes(text: str) -> dict:
    """Parse a GFF3 column-9 string into an ordered dict of decoded values."""
    attrs = {}
    for token in text.strip().strip(";").split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            attrs[key.strip()] = urllib.parse.unquote(value.strip())
        else:
            attrs[token] = ""
    return attrs


def encode_attribute_value(value: str) -> str:
    """
    Percent-encode a GFF3 attribute value. Tab, newline, carriage return,
    %, ;, = and & are encoded; spaces, commas and other punctuation are
    kept readable, which matches the style of NCBI GFF3 files (commas
    appear unencoded in product names).
    """
    return urllib.parse.quote(str(value), safe=" ,:._-|/()[]'*+!@$^{}<>?#~`\"")


@dataclass
class Transcript:
    id: str
    ftype: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: dict = field(default_factory=dict)
    parent_id: Optional[str] = None
    exons: list = field(default_factory=list)      # [(start, end)] ascending
    cds: list = field(default_factory=list)        # [(start, end, phase)] ascending
    cds_attrs: dict = field(default_factory=dict)  # attributes of the first CDS line
    key: str = ""                                   # identifier used in outputs
    implicit: bool = False                          # created for a gene with direct exon/CDS children

    @property
    def coding(self) -> bool:
        return bool(self.cds)

    @property
    def cds_length(self) -> int:
        return sum(e - s + 1 for s, e, _ in self.cds)

    @property
    def transcript_length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)

    @property
    def first_phase(self) -> int:
        """Phase of the first CDS segment in translation order."""
        if not self.cds:
            return 0
        ordered = sorted(self.cds, key=lambda c: c[0], reverse=(self.strand == "-"))
        phase = ordered[0][2]
        return phase if phase in (0, 1, 2) else 0


@dataclass
class Gene:
    id: str
    ftype: str
    seqid: str
    strand: str
    start: int
    end: int
    attrs: dict = field(default_factory=dict)
    transcripts: list = field(default_factory=list)
    key: str = ""
    wrapper: bool = False      # created for a transcript without a gene parent

    @property
    def coding(self) -> bool:
        return any(t.coding for t in self.transcripts)

    @property
    def biotype(self) -> str:
        """gene_biotype, biotype or gene_type attribute, else the feature type."""
        return (self.attrs.get("gene_biotype") or self.attrs.get("biotype")
                or self.attrs.get("gene_type") or self.ftype)


@dataclass
class GFF3Stats:
    feature_counts: Counter = field(default_factory=Counter)
    n_genes: int = 0
    n_wrapper_genes: int = 0
    n_genes_without_transcripts: int = 0
    n_transcripts: int = 0
    n_coding_transcripts: int = 0
    n_implicit_transcripts: int = 0
    n_nested_features: int = 0
    n_orphan_children: int = 0
    n_ignored_features: int = 0
    n_coding_genes: int = 0
    format: str = "GFF3"
    n_exons_as_cds: int = 0


def parse_gff3(path: str) -> tuple:
    """
    Parse a GFF3 file.

    Returns (genes, stats): genes is a list of Gene objects in order of
    first appearance in the file.
    """
    logger.info(f"Parsing GFF3 from {path}")

    gene_records = {}
    gene_order = []
    candidates = {}            # possible transcript records keyed by ID
    candidate_order = []
    exon_records = defaultdict(list)
    cds_records = defaultdict(list)
    cds_attrs = {}
    stats = GFF3Stats()
    noid_counter = 0

    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr_text = fields
            stats.feature_counts[ftype] += 1
            attrs = parse_attributes(attr_text)
            fid = attrs.get("ID", "")
            parents = [p for p in attrs.get("Parent", "").split(",") if p]
            s = int(start)
            e = int(end)

            if ftype == "exon":
                if not parents:
                    stats.n_orphan_children += 1
                for p in parents:
                    exon_records[p].append((s, e))
            elif ftype == "CDS":
                ph = int(phase) if phase in ("0", "1", "2") else 0
                if not parents:
                    stats.n_orphan_children += 1
                for p in parents:
                    cds_records[p].append((s, e, ph))
                    if p not in cds_attrs:
                        cds_attrs[p] = {k: v for k, v in attrs.items() if k not in ("ID", "Parent")}
            elif ftype in DERIVED_TYPES:
                continue
            elif is_gene_type(ftype):
                if not fid:
                    noid_counter += 1
                    fid = f"gene_noid_{noid_counter}"
                gene_records[fid] = {
                    "seqid": seqid, "ftype": ftype, "start": s, "end": e,
                    "strand": strand, "attrs": attrs,
                }
                gene_order.append(fid)
            elif ftype in IGNORED_TYPES and not parents:
                stats.n_ignored_features += 1
            else:
                if not fid:
                    if not parents:
                        stats.n_ignored_features += 1
                        continue
                    noid_counter += 1
                    fid = f"{ftype}_noid_{noid_counter}"
                candidates[fid] = {
                    "seqid": seqid, "ftype": ftype, "start": s, "end": e,
                    "strand": strand, "attrs": attrs,
                    "parent": parents[0] if parents else None,
                }
                candidate_order.append(fid)

    # Classify candidates: transcript of a gene, transcript without a gene,
    # nested under another candidate, or unusable.
    genes = {}
    for gid in gene_order:
        g = gene_records[gid]
        genes[gid] = Gene(id=gid, ftype=g["ftype"], seqid=g["seqid"], strand=g["strand"],
                          start=g["start"], end=g["end"], attrs=g["attrs"])

    transcripts_by_gene = defaultdict(list)
    for cid in candidate_order:
        c = candidates[cid]
        parent = c["parent"]
        if parent is not None and parent in genes:
            gene_id = parent
        elif parent is not None and parent in candidates:
            stats.n_nested_features += 1
            continue
        else:
            has_children = cid in exon_records or cid in cds_records
            if not has_children and not is_transcript_type(c["ftype"]):
                stats.n_ignored_features += 1
                continue
            if parent is not None:
                stats.n_orphan_children += 1
                logger.warning(f"Transcript {cid} references unknown parent {parent}; "
                               f"treated as a gene of its own")
            gene_id = cid
            if gene_id not in genes:
                genes[gene_id] = Gene(id=gene_id, ftype="gene", seqid=c["seqid"], strand=c["strand"],
                                      start=c["start"], end=c["end"], attrs={}, wrapper=True)
                gene_order.append(gene_id)
                stats.n_wrapper_genes += 1

        tx = _build_transcript(cid, c, gene_id, exon_records, cds_records, cds_attrs)
        transcripts_by_gene[gene_id].append(tx)

    # Genes with exon or CDS children and no transcript: implicit transcript
    for gid in gene_order:
        if gid in exon_records or gid in cds_records:
            tid = f"{gid}.mRNA"
            g = gene_records.get(gid)
            if g is None:
                continue
            c = {"seqid": g["seqid"], "ftype": "mRNA", "start": g["start"], "end": g["end"],
                 "strand": g["strand"], "attrs": {}, "parent": gid}
            tx = _build_transcript(tid, c, gid, exon_records, cds_records, cds_attrs)
            tx.implicit = True
            transcripts_by_gene[gid].append(tx)
            stats.n_implicit_transcripts += 1

    out = []
    for gid in gene_order:
        gene = genes[gid]
        gene.transcripts = sorted(transcripts_by_gene.get(gid, []), key=lambda t: (t.start, t.end, t.id))
        if not gene.transcripts:
            stats.n_genes_without_transcripts += 1
        out.append(gene)

    stats.n_genes = len(out)
    stats.n_coding_genes = sum(1 for g in out if g.coding)
    stats.n_transcripts = sum(len(g.transcripts) for g in out)
    stats.n_coding_transcripts = sum(1 for g in out for t in g.transcripts if t.coding)

    logger.info(f"  Features by type: " + ", ".join(
        f"{k}={v:,}" for k, v in stats.feature_counts.most_common()))
    logger.info(f"  Genes: {stats.n_genes:,} ({stats.n_coding_genes:,} coding, "
                f"{stats.n_wrapper_genes:,} created for parentless transcripts, "
                f"{stats.n_genes_without_transcripts:,} without transcripts)")
    logger.info(f"  Transcripts: {stats.n_transcripts:,} ({stats.n_coding_transcripts:,} coding, "
                f"{stats.n_implicit_transcripts:,} implicit)")
    if stats.n_nested_features:
        logger.info(f"  Nested features ignored: {stats.n_nested_features:,}")
    if stats.n_orphan_children:
        logger.warning(f"  Children with missing parents: {stats.n_orphan_children:,}")
    if stats.n_ignored_features:
        logger.info(f"  Non-model features ignored: {stats.n_ignored_features:,}")

    return out, stats


def _build_transcript(tid, rec, gene_id, exon_records, cds_records, cds_attrs) -> Transcript:
    exons = sorted(set(exon_records.get(tid, [])))
    cds_segs = cds_records.get(tid, [])
    cds = []
    if cds_segs:
        seen = set()
        for s, e, ph in sorted(cds_segs):
            if (s, e) in seen:
                continue
            seen.add((s, e))
            cds.append((s, e, ph))
    if not exons:
        if cds:
            exons = [(s, e) for s, e, _ in cds]
        else:
            exons = [(rec["start"], rec["end"])]
    return Transcript(
        id=tid, ftype=rec["ftype"], seqid=rec["seqid"], strand=rec["strand"],
        start=rec["start"], end=rec["end"], attrs=rec["attrs"], parent_id=gene_id,
        exons=exons, cds=cds, cds_attrs=cds_attrs.get(tid, {}),
    )


def assign_keys(genes: list, gene_id_attr: str = "ID", transcript_id_attr: str = "ID") -> None:
    """
    Set the identifiers used in all outputs (FASTA headers, mapping files).

    By default the GFF3 ID is used. Another attribute (Name, locus_tag, ...)
    can be chosen; it must be unique across genes and across transcripts.
    """
    seen_genes = {}
    seen_tx = {}
    for gene in genes:
        key = gene.id if gene_id_attr == "ID" else gene.attrs.get(gene_id_attr, "")
        if not key:
            key = gene.id
        if key in seen_genes and seen_genes[key] != gene.id:
            raise GOAnnotateError(
                f"Gene identifier attribute '{gene_id_attr}' is not unique: '{key}' is used by "
                f"{seen_genes[key]} and {gene.id}. Use --gene-id-attr ID.")
        seen_genes[key] = gene.id
        gene.key = key
        for tx in gene.transcripts:
            tkey = tx.id if transcript_id_attr == "ID" else tx.attrs.get(transcript_id_attr, "")
            if not tkey:
                tkey = tx.id
            if tkey in seen_tx and seen_tx[tkey] != tx.id:
                raise GOAnnotateError(
                    f"Transcript identifier attribute '{transcript_id_attr}' is not unique: "
                    f"'{tkey}' is used by {seen_tx[tkey]} and {tx.id}.")
            seen_tx[tkey] = tx.id
            tx.key = tkey


def universe_from_genes(genes: list) -> GeneUniverse:
    """Build the gene universe (every gene, coding or not) from parsed genes."""
    gene_ids = []
    transcript_to_gene = {}
    gene_to_transcripts = {}
    gene_types = {}
    coding = {}
    total_sequences = 0
    for gene in genes:
        gene_ids.append(gene.key)
        gene_to_transcripts[gene.key] = [t.key for t in gene.transcripts]
        gene_types[gene.key] = gene.biotype
        coding[gene.key] = gene.coding
        for tx in gene.transcripts:
            transcript_to_gene[tx.key] = gene.key
            total_sequences += 1

    universe = GeneUniverse(gene_ids=gene_ids, transcript_to_gene=transcript_to_gene,
                            gene_to_transcripts=gene_to_transcripts,
                            total_sequences=total_sequences, gene_types=gene_types, coding=coding)
    n_multi = sum(1 for txs in gene_to_transcripts.values() if len(txs) > 1)
    logger.info(f"Gene universe from the annotation: {len(gene_ids):,} genes, {total_sequences:,} transcripts, "
                f"{n_multi:,} genes with multiple transcripts")
    return universe
