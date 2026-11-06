"""Proteins for nucleotide transcripts, translated in the frame of their best Diamond blastx hit."""

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from goannotate import sequences

logger = logging.getLogger(__name__)

# Hits considered when checking whether a transcript's hits fall in more
# than one frame or in separate regions of the transcript
FRAME_CHECK_HITS = 10


@dataclass
class OrfRecord:
    transcript_id: str
    strand: str              # + or -
    frame: int               # 1-3 on the given strand
    orf_start: int           # 1-based, inclusive, on the transcript as given (start < end)
    orf_end: int
    protein_length: int
    completeness: str        # complete, 5_partial, 3_partial, partial
    subject_id: str          # hit that defined the frame
    bitscore: float
    n_hits: int
    note: str = ""           # multiple_frames, separate_regions, internal_stop, no_start_codon


@dataclass
class TranslationResult:
    protein_fasta: Path
    report: Path
    n_transcripts: int = 0
    n_with_hits: int = 0
    n_translated: int = 0
    completeness: Counter = field(default_factory=Counter)
    notes: Counter = field(default_factory=Counter)
    n_short: int = 0         # proteins shorter than SHORT_PROTEIN residues
    n_blastx_fallback: int = 0   # transcripts annotated from their blastx hits (no blastp hit)
    final_completeness: Counter = field(default_factory=Counter)
    confidence: Counter = field(default_factory=Counter)   # confidence of the translated genes, fragments lowered


SHORT_PROTEIN = 50


def _read_hits(blastx_results: str, evalue: float) -> dict:
    """query_id -> list of (bitscore, evalue, qstart, qend, qlen, subject_id) under the e-value threshold."""
    hits = defaultdict(list)
    with open(blastx_results) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15 or not (parts[12].isdigit() and parts[13].isdigit()):
                continue
            ev = float(parts[10])
            if ev > evalue:
                continue
            hits[parts[0]].append((float(parts[11]), ev, int(parts[6]), int(parts[7]), int(parts[12]), parts[1]))
    return hits


def _oriented(qstart: int, qend: int, qlen: int) -> tuple:
    """(strand, start, end) of an alignment on the strand it lies on, 1-based, start < end."""
    if qstart <= qend:
        return "+", qstart, qend
    return "-", qlen - qstart + 1, qlen - qend + 1


def orf_from_hit(seq: str, qstart: int, qend: int, table_id: int = 1) -> tuple:
    """
    The protein encoded around a blastx alignment: the open reading frame in
    the alignment's frame that contains it, extended upstream to the first
    methionine after the preceding in-frame stop (or to the transcript start)
    and downstream to the first in-frame stop (or to the transcript end).

    Returns (strand, frame, orf_start, orf_end, protein, completeness, note)
    with the coordinates on the given sequence, 1-based and inclusive, the
    stop codon included when present. completeness is complete, 5_partial
    (no methionine at the start), 3_partial (no stop) or partial (neither).
    note lists internal_stop and no_start_codon, and upstream_stop when an
    in-frame stop precedes the ORF (the methionine is then the first
    possible start, not one inside a longer reading frame cut by the
    transcript end).
    """
    qlen = len(seq)
    strand, a_nt, b_nt = _oriented(qstart, qend, qlen)
    oriented = seq if strand == "+" else sequences.reverse_complement(seq)
    offset = (a_nt - 1) % 3
    protein = sequences.translate(oriented[offset:], table_id)
    a = (a_nt - 1 - offset) // 3
    b = min((b_nt - 1 - offset) // 3, len(protein) - 1)
    note = []

    # a stop inside the aligned region ends the protein there
    inner_stop = protein.find("*", a, b + 1)
    if inner_stop != -1:
        note.append("internal_stop")
    # upstream: the first methionine after the preceding in-frame stop
    upstream_stop = protein.rfind("*", 0, a)
    region_start = upstream_stop + 1
    met = protein.find("M", region_start, a + 1)
    if met != -1:
        start = met
        if upstream_stop != -1:
            note.append("upstream_stop")
    else:
        start = region_start
        if upstream_stop != -1:
            note.append("no_start_codon")
    # downstream: the first in-frame stop
    end = protein.find("*", a) if inner_stop == -1 else inner_stop
    if end == -1:
        aa = protein[start:]
        last_nt = offset + 3 * len(protein)      # last complete codon
        has_stop = False
    else:
        aa = protein[start:end]
        last_nt = offset + 3 * (end + 1)         # stop codon included
        has_stop = True
    first_nt = offset + 3 * start + 1
    if aa.startswith("M") and has_stop:
        completeness = "complete"
    elif aa.startswith("M"):
        completeness = "3_partial"
    elif has_stop:
        completeness = "5_partial"
    else:
        completeness = "partial"
    if strand == "+":
        orf_start, orf_end = first_nt, last_nt
    else:
        orf_start, orf_end = qlen - last_nt + 1, qlen - first_nt + 1
    return strand, offset + 1, orf_start, orf_end, aa, completeness, ";".join(note)


def _frame_notes(hits: list, qlen: int) -> list:
    """Notes on the top hits: more than one frame, or alignments outside the best hit's region."""
    top = sorted(hits, key=lambda h: (-h[0], h[1]))[:FRAME_CHECK_HITS]
    best = _oriented(top[0][2], top[0][3], qlen)
    frames = set()
    separate = False
    for h in top:
        strand, s, e = _oriented(h[2], h[3], qlen)
        frames.add((strand, (s - 1) % 3))
        if strand != best[0] or e < best[1] or s > best[2]:
            separate = True
    notes = []
    if len(frames) > 1:
        notes.append("multiple_frames")
    if separate:
        notes.append("separate_regions")
    return notes


def translate_from_hits(transcripts: str, blastx_results: str, protein_fasta: Path, report: Path,
                        table_id: int = 1, evalue: float = 1e-5) -> TranslationResult:
    """
    Write one protein per transcript with a blastx hit: the open reading
    frame around the best hit (orf_from_hit). Transcripts without a hit get
    no protein. The report lists each translated transcript with its
    strand, frame, ORF coordinates, completeness and the defining hit.
    """
    logger.info(f"Translating transcripts in the frame of their best blastx hit: {transcripts}")
    hits = _read_hits(blastx_results, evalue)
    result = TranslationResult(protein_fasta=Path(protein_fasta), report=Path(report))
    Path(protein_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(report).parent.mkdir(parents=True, exist_ok=True)
    with open(protein_fasta, "w") as f_prot, open(report, "w") as f_rep:
        f_rep.write("transcript_id\tstrand\tframe\torf_start\torf_end\tprotein_length\tcompleteness\t"
                    "defining_hit\tbitscore\thits\tnote\n")
        for header, seq in sequences.read_fasta(transcripts):
            result.n_transcripts += 1
            tid = header.split()[0] if header else ""
            tx_hits = hits.get(tid)
            if not tx_hits:
                continue
            result.n_with_hits += 1
            best = max(tx_hits, key=lambda h: (h[0], -h[1]))
            qlen = len(seq)
            if best[4] != qlen:
                logger.warning(f"{tid}: the search recorded a length of {best[4]:,} but the transcript has "
                               f"{qlen:,} nucleotides; the hit coordinates may not apply")
            strand, frame, orf_start, orf_end, aa, completeness, note = orf_from_hit(seq, best[2], best[3], table_id)
            if not aa:
                continue
            notes = _frame_notes(tx_hits, qlen) + ([note] if note else [])
            rec = OrfRecord(tid, strand, frame, orf_start, orf_end, len(aa), completeness, best[5], best[0],
                            len(tx_hits), ";".join(notes))
            result.n_translated += 1
            result.completeness[completeness] += 1
            for n in notes:
                for part in n.split(";"):
                    if part != "upstream_stop":
                        result.notes[part] += 1
            if len(aa) < SHORT_PROTEIN:
                result.n_short += 1
            fields = {"orf": f"{orf_start}-{orf_end}", "strand": strand, "frame": frame, "orf_status": completeness}
            sequences._write_fasta_record(f_prot, sequences.format_header(tid, fields), aa)
            f_rep.write(f"{rec.transcript_id}\t{rec.strand}\t{rec.frame}\t{rec.orf_start}\t{rec.orf_end}\t"
                        f"{rec.protein_length}\t{rec.completeness}\t{rec.subject_id}\t{rec.bitscore:g}\t"
                        f"{rec.n_hits}\t{rec.note}\n")
    logger.info(f"  {result.n_translated:,} of {result.n_transcripts:,} transcripts translated "
                f"({result.n_with_hits:,} with hits); " +
                ", ".join(f"{k} {v:,}" for k, v in result.completeness.most_common()) +
                (("; flagged: " + ", ".join(f"{k} {v:,}" for k, v in result.notes.most_common()))
                 if result.notes else ""))
    logger.info(f"  Wrote {protein_fasta}; report {report}")
    return result


def read_report(report: Path) -> dict:
    """transcript_id -> the report row as a dict."""
    rows = {}
    with open(report) as fh:
        header = next(fh, "").rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 7:
                rows[p[0]] = dict(zip(header, p + [""] * (len(header) - len(p))))
    return rows


def orf_header_fields(report: Path, completeness: dict = None) -> dict:
    """
    transcript_id -> ORF header fields (orf, strand, frame, orf_status) read
    back from the report, with completeness= when completeness
    (transcript_id -> (status, evidence)) is given.
    """
    fields = {}
    for tid, r in read_report(report).items():
        fields[tid] = {"orf": f"{r['orf_start']}-{r['orf_end']}", "strand": r["strand"], "frame": r["frame"],
                       "orf_status": r["completeness"]}
        if completeness and tid in completeness:
            fields[tid]["completeness"] = completeness[tid][0]
    return fields


def blastx_fallback(blastp_results: Path, blastx_results: str, report: Path, out: Path) -> tuple:
    """
    Search results for the annotation: the blastp results of the proteins,
    plus the blastx hits of the translated transcripts whose protein found
    no blastp hit (short proteins with hits near the e-value threshold).
    The blastx query coordinates and length are given in residues of the
    transcript's reading frame, so that the coverage filter reads them as
    it reads blastp hits. Returns (path, number of transcripts that use
    their blastx hits).
    """
    with_hits = set()
    with open(blastp_results) as fh:
        for line in fh:
            with_hits.add(line.split("\t", 1)[0])
    translated = set(read_report(report))
    missing = translated - with_hits
    if not missing:
        return Path(blastp_results), 0
    tmp = Path(str(out) + ".tmp")
    used = set()
    with open(tmp, "w") as fo:
        with open(blastp_results) as fh:
            for line in fh:
                fo.write(line)
        with open(blastx_results) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts[0] not in missing or len(parts) < 15 or not parts[12].isdigit():
                    continue
                qs, qe, qlen = int(parts[6]), int(parts[7]), int(parts[12])
                if qs > qe:
                    qs, qe = qlen - qs + 1, qlen - qe + 1
                parts[6], parts[7] = str((qs + 2) // 3), str(qe // 3)
                parts[12] = str(qlen // 3)
                fo.write("\t".join(parts) + "\n")
                used.add(parts[0])
    tmp.replace(out)
    logger.info(f"  {len(used):,} translated transcripts without blastp hits are annotated from their blastx hits")
    return out, len(used)


# Residues by which an alignment may stop short of the reference protein's
# ends and still count as reaching them (N- and C-termini differ between
# orthologs by a few residues)
END_TOLERANCE = 10


def completeness(report: Path, ref_alignments: dict, annotations: dict, tx2gene: dict) -> dict:
    """
    Completeness of each translated transcript's protein:
      complete    both ends present
      5_partial   the N-terminus is missing
      3_partial   the C-terminus is missing
      partial     both ends missing
    The 3' end is present when the ORF ends in a stop codon. The 5' end is
    present when the ORF begins with a methionine and either an in-frame
    stop precedes it in the transcript, or the alignment to the reference
    protein (that of the gene's symbol, else the gene's best reference hit)
    begins within END_TOLERANCE residues of the reference start.
    ref_alignments: transcript -> {reference symbol: (bits, sstart, send, slen)}.
    Returns transcript_id -> (status, evidence); evidence names what supports
    or contradicts each end (5:met+upstream_stop, 5:met+reference_start,
    5:no_met, 5:met_without_support, 3:stop, 3:no_stop) and, when the
    reference alignment reaches the reference end without a stop codon
    in the transcript, 3:reference_end.
    """
    out = {}
    for tid, r in read_report(report).items():
        has_met = r["completeness"] in ("complete", "3_partial")
        has_stop = r["completeness"] in ("complete", "5_partial")
        bounded = "upstream_stop" in r.get("note", "")
        annot = annotations.get(tx2gene.get(tid, tid))
        aligns = ref_alignments.get(tid, {})
        ref = None
        symbol = (annot.consensus_symbol or "").upper() if annot else ""
        if symbol and symbol in aligns:
            ref = aligns[symbol]
        elif aligns:
            ref = max(aligns.values())
        ev = []
        if has_met and bounded:
            five = True
            ev.append("5:met+upstream_stop")
        elif has_met and ref and ref[1] <= END_TOLERANCE + 1:
            five = True
            ev.append("5:met+reference_start")
        elif has_met:
            five = False
            ev.append("5:met_without_support")
        else:
            five = False
            ev.append("5:no_met")
        ev.append("3:stop" if has_stop else "3:no_stop")
        if not has_stop and ref and ref[2] >= ref[3] - END_TOLERANCE:
            ev.append("3:reference_end")
        status = ("complete" if five and has_stop else "3_partial" if five
                  else "5_partial" if has_stop else "partial")
        out[tid] = (status, ",".join(ev))
    return out
