"""GTF parsing into the gene > transcript > exon/CDS hierarchy of goannotate.gff3.

Features are grouped by their ``gene_id`` and ``transcript_id`` attributes.
``gene`` and ``transcript`` lines are used when present; a gene or a
transcript without a line of its own is built from the span of its exon
and CDS lines. Column 9 holding a bare identifier on gene and transcript
lines (AUGUSTUS) is accepted. ``stop_codon`` features are joined to the
CDS, so that the CDS includes the stop codon as it does in GFF3. A
transcript without CDS lines whose exon lines all carry a frame (0, 1 or
2; GTF defines the frame only for coding features) has its exons read as
the coding sequence, as written by gene predictors that emit coding
segments only (GeMoMa). The frame
column is read as the GFF3 phase. A key that occurs several times on a
line (``tag``, ``db_xref``) is kept as a comma-separated value.

Coordinates are 1-based and inclusive, as in GFF3.
"""

import logging
from collections import defaultdict

from goannotate.gff3 import (DERIVED_TYPES, GFF3Stats, Gene, Transcript, is_gene_type,
                             is_transcript_type)

logger = logging.getLogger(__name__)

# Attributes of a CDS line that identify the line and are not kept as CDS attributes
_LINE_KEYS = {"gene_id", "transcript_id", "exon_number", "exon_id", "exon_version"}


def split_gtf_attributes(text: str) -> list:
    """Split column 9 into [key, raw_value] pairs; raw values keep their quotes."""
    tokens = []
    current = []
    in_quotes = False
    previous = ""
    for ch in text:
        if ch == '"' and previous != "\\":
            in_quotes = not in_quotes
        if ch == ";" and not in_quotes:
            tokens.append("".join(current))
            current = []
        else:
            current.append(ch)
        previous = ch
    tokens.append("".join(current))

    pairs = []
    for token in tokens:
        token = token.strip()
        if not token:
            continue
        parts = token.split(None, 1)
        pairs.append([parts[0], parts[1].strip() if len(parts) > 1 else None])
    return pairs


def unquote_gtf_value(raw) -> str:
    if raw is None:
        return ""
    raw = raw.strip()
    if len(raw) >= 2 and raw[0] == '"' and raw[-1] == '"':
        raw = raw[1:-1]
    return raw.replace('\\"', '"')


def bare_identifier(text: str) -> str:
    """Column 9 when it is a single identifier without a key (AUGUSTUS gene and transcript lines)."""
    text = text.strip().rstrip(";").strip()
    if text and not any(c in text for c in ' \t";'):
        return text
    return ""


def gtf_line_pairs(ftype: str, text: str) -> list:
    """Column 9 as [key, raw_value] pairs; a bare identifier becomes gene_id or transcript_id."""
    bare = bare_identifier(text)
    if bare:
        key = "gene_id" if is_gene_type(ftype) else "transcript_id"
        return [[key, f'"{bare}"']]
    return split_gtf_attributes(text)


def parse_gtf_attributes(ftype: str, text: str) -> dict:
    attrs = {}
    for key, raw in gtf_line_pairs(ftype, text):
        value = unquote_gtf_value(raw)
        if key in attrs and attrs[key] and value:
            attrs[key] = attrs[key] + "," + value
        elif key not in attrs or value:
            attrs[key] = value
    return attrs


def _join_stop_codons(cds: list, stops: list, strand: str) -> list:
    """
    Extend the CDS segment that a stop_codon segment overlaps or adjoins;
    a stop_codon segment apart from every CDS segment becomes a segment of
    its own. CDS segments are never joined with one another.
    """
    out = list(cds)
    for s, e, ph in sorted(set(stops)):
        for i, (cs, ce, cph) in enumerate(out):
            if s <= ce + 1 and e >= cs - 1:
                # The joined segment keeps the phase of its part translated first
                if strand == "-" and e > ce:
                    cph = ph
                elif strand != "-" and s < cs:
                    cph = ph
                out[i] = (min(cs, s), max(ce, e), cph)
                break
        else:
            out.append((s, e, ph))
    return sorted(out)


def parse_gtf(path: str) -> tuple:
    """
    Parse a GTF file.

    Returns (genes, stats) like gff3.parse_gff3: genes is a list of Gene
    objects in order of first appearance in the file.
    """
    logger.info(f"Parsing GTF from {path}")

    gene_lines = {}
    gene_order = []
    gene_seen = set()
    tx_lines = {}
    tx_order = []
    tx_seen = set()
    tx_gene = {}
    first_attrs = {}                 # transcript or gene ID -> attributes of its first child line
    spans = {}                       # transcript ID -> [seqid, strand, start, end] over its child lines
    exon_records = defaultdict(list)
    exon_frames = defaultdict(list)    # (start, end, frame or None) of exon lines
    cds_records = defaultdict(list)
    stop_records = defaultdict(list)
    cds_attrs = {}
    implicit_tx = set()
    stats = GFF3Stats()
    stats.format = "GTF"

    def see_gene(gid):
        if gid not in gene_seen:
            gene_seen.add(gid)
            gene_order.append(gid)

    def see_transcript(tid, gid):
        if tid not in tx_seen:
            tx_seen.add(tid)
            tx_order.append(tid)
        if gid and tid not in tx_gene:
            tx_gene[tid] = gid

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, frame, attr_text = fields
            stats.feature_counts[ftype] += 1
            attrs = parse_gtf_attributes(ftype, attr_text)
            gid = attrs.get("gene_id", "")
            tid = attrs.get("transcript_id", "")
            s = int(start)
            e = int(end)

            if is_gene_type(ftype):
                if not gid:
                    stats.n_ignored_features += 1
                    continue
                see_gene(gid)
                if gid not in gene_lines:
                    gene_lines[gid] = {"seqid": seqid, "ftype": ftype, "start": s, "end": e,
                                       "strand": strand, "attrs": attrs}
                continue

            if ftype in ("exon", "CDS", "stop_codon"):
                if not tid:
                    if not gid:
                        stats.n_orphan_children += 1
                        continue
                    tid = f"{gid}.mRNA"
                    implicit_tx.add(tid)
                if gid:
                    see_gene(gid)
                see_transcript(tid, gid)
                first_attrs.setdefault(tid, attrs)
                if gid:
                    first_attrs.setdefault(("gene", gid), attrs)
                span = spans.get(tid)
                if span is None:
                    spans[tid] = [seqid, strand, s, e]
                else:
                    span[2] = min(span[2], s)
                    span[3] = max(span[3], e)
                if ftype == "exon":
                    exon_records[tid].append((s, e))
                    exon_frames[tid].append((s, e, int(frame) if frame in ("0", "1", "2") else None))
                else:
                    ph = int(frame) if frame in ("0", "1", "2") else 0
                    if ftype == "CDS":
                        cds_records[tid].append((s, e, ph))
                        if tid not in cds_attrs:
                            cds_attrs[tid] = {k: v for k, v in attrs.items() if k not in _LINE_KEYS}
                    else:
                        stop_records[tid].append((s, e, ph))
            elif ftype in DERIVED_TYPES:
                continue
            elif tid and (ftype == "transcript" or is_transcript_type(ftype)):
                if gid:
                    see_gene(gid)
                see_transcript(tid, gid)
                if tid not in tx_lines:
                    tx_lines[tid] = {"seqid": seqid, "ftype": ftype, "start": s, "end": e,
                                     "strand": strand, "attrs": attrs}
            else:
                stats.n_ignored_features += 1

    # Transcripts
    transcripts_by_gene = defaultdict(list)
    for tid in tx_order:
        gid = tx_gene.get(tid)
        if not gid:
            # A transcript that names no gene is a gene of its own
            gid = tid
            tx_gene[tid] = gid
            see_gene(gid)
        rec = tx_lines.get(tid)
        inferred = rec is None
        if inferred:
            if tid not in spans:
                continue
            seqid, strand, s, e = spans[tid]
            child = first_attrs.get(tid, {})
            attrs = {k: v for k, v in child.items()
                     if k.startswith("transcript") or k == "product"}
            rec = {"seqid": seqid, "ftype": "transcript", "start": s, "end": e,
                   "strand": strand, "attrs": attrs}
            stats.n_implicit_transcripts += 1
        attrs = rec["attrs"]
        ftype = attrs.get("transcript_biotype") or attrs.get("transcript_type") or rec["ftype"]

        exons = sorted(set(exon_records.get(tid, [])))
        cds = []
        seen = set()
        for s, e, ph in sorted(cds_records.get(tid, [])):
            if (s, e) in seen:
                continue
            seen.add((s, e))
            cds.append((s, e, ph))
        if not cds and exon_frames.get(tid) and all(ph is not None for _, _, ph in exon_frames[tid]):
            # No CDS lines and a frame on every exon line: the exons are the
            # coding segments (GTF defines the frame only for coding features)
            cds = sorted(set(exon_frames[tid]))
            stats.n_exons_as_cds += 1
        if cds:
            cds = _join_stop_codons(cds, stop_records.get(tid, []), rec["strand"])
        if not exons:
            exons = [(s, e) for s, e, _ in cds] if cds else [(rec["start"], rec["end"])]
        tx = Transcript(id=tid, ftype=ftype, seqid=rec["seqid"], strand=rec["strand"],
                        start=rec["start"], end=rec["end"], attrs=attrs, parent_id=gid,
                        exons=exons, cds=cds, cds_attrs=cds_attrs.get(tid, {}))
        # implicit: the transcript has no line of its own in the file
        tx.implicit = inferred
        transcripts_by_gene[gid].append(tx)

    # Genes
    out = []
    for gid in gene_order:
        txs = sorted(transcripts_by_gene.get(gid, []), key=lambda t: (t.start, t.end, t.id))
        rec = gene_lines.get(gid)
        if rec is not None:
            gene = Gene(id=gid, ftype=rec["ftype"], seqid=rec["seqid"], strand=rec["strand"],
                        start=rec["start"], end=rec["end"], attrs=rec["attrs"])
        else:
            if not txs:
                continue
            child = first_attrs.get(("gene", gid), {})
            attrs = {k: v for k, v in child.items() if k.startswith("gene") or k == "description"}
            gene = Gene(id=gid, ftype="gene", seqid=txs[0].seqid, strand=txs[0].strand,
                        start=min(t.start for t in txs), end=max(t.end for t in txs),
                        attrs=attrs, wrapper=True)
            stats.n_wrapper_genes += 1
        gene.transcripts = txs
        if not txs:
            stats.n_genes_without_transcripts += 1
        out.append(gene)

    stats.n_genes = len(out)
    stats.n_coding_genes = sum(1 for g in out if g.coding)
    stats.n_transcripts = sum(len(g.transcripts) for g in out)
    stats.n_coding_transcripts = sum(1 for g in out for t in g.transcripts if t.coding)

    logger.info("  Features by type: " + ", ".join(
        f"{k}={v:,}" for k, v in stats.feature_counts.most_common()))
    logger.info(f"  Genes: {stats.n_genes:,} ({stats.n_coding_genes:,} coding, "
                f"{stats.n_wrapper_genes:,} without a gene line, "
                f"{stats.n_genes_without_transcripts:,} without transcripts)")
    logger.info(f"  Transcripts: {stats.n_transcripts:,} ({stats.n_coding_transcripts:,} coding, "
                f"{stats.n_implicit_transcripts:,} without a transcript line)")
    if stats.n_exons_as_cds:
        logger.warning(f"  Transcripts without CDS lines whose exon lines all carry a frame: "
                       f"{stats.n_exons_as_cds:,}; their exons are read as the coding sequence")
    if (stats.feature_counts.get("CDS") or stats.n_exons_as_cds) and not stats.feature_counts.get("stop_codon"):
        logger.warning("  The GTF has no stop_codon features; the coding sequence is used as given")
    if stats.n_orphan_children:
        logger.warning(f"  Exon or CDS lines without gene_id and transcript_id: {stats.n_orphan_children:,}")
    if stats.n_ignored_features:
        logger.info(f"  Non-model features ignored: {stats.n_ignored_features:,}")

    return out, stats
