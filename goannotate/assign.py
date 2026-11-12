"""Genome-wide symbol checks against the reference proteome.

The per-gene vote (annotate.py) decides identity from many database hits.
Its weakness is naming: a gene that shares a domain with a well-studied
protein can be named after it, and for paralog pairs the better-studied
sibling tends to win the vote. This module checks the vote against the
complete reference-species proteome:

1. The reference proteome (SwissProt entries of the reference taxon, taken
   from the SwissProt Diamond database or a user-supplied FASTA) is
   searched against the query proteins and the query proteins against it.
   A gene is the reciprocal best hit (RBH) of reference symbol X when X's
   protein has that gene as its best query hit and X is among the gene's
   hits.
2. A voted symbol equal to the RBH symbol is confirmed. An RBH symbol that
   differs but occurs among the resolved symbols of the gene's winning
   cluster replaces the vote. A gene the vote left unnamed takes its RBH
   symbol.
3. Every other voted reference symbol must be supported by the gene's own
   sequence: the reference protein of that symbol must score within a set
   fraction of the gene's best reference hit, and the alignment must cover
   at least half of the reference protein. Unsupported votes are replaced
   by the closest supported reference protein, or removed when no
   reference protein is covered (a shared domain or a fragment). Genes
   without any reference hit keep their vote.
4. Voted symbols outside the reference nomenclature (a zebrafish foxo3a, a
   rodent Vmn1r) are, under the default symbol policy "reference", replaced
   by the closest reference protein that the gene covers, or removed when
   none is covered, so that every symbol in the output belongs to one
   nomenclature. Under the policy "any" they are kept.
5. Paralog designators of other genomes left in voted product names
   ("Forkhead box protein O3-B", "protein kinase C beta a") are removed
   when the name without them is a reference protein name or a name among
   the gene's own hits.

Several genes may carry the same symbol. Paralogs, recent duplicates and
retrocopies are real, and each copy is named for what it is homologous to.
"""

import logging
import os
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

from goannotate import GOAnnotateError
from goannotate.mappings import gene_symbol_of
from goannotate.search import find_diamond


logger = logging.getLogger(__name__)


####################
# REFERENCE PROTEOME
####################

def _getseq_taxon(diamond_bin: str, db: str, taxid: str, fout) -> int:
    """Write the entries of one taxon (OX= in the header) of a Diamond database to fout; returns the count."""
    tag = f"OX={taxid} "
    n = 0
    proc = subprocess.Popen([diamond_bin, "getseq", "--db", db],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    keep = False
    for line in proc.stdout:
        if line.startswith(">"):
            keep = tag in line
            if keep:
                n += 1
        if keep:
            fout.write(line)
    proc.wait()
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.read().strip().splitlines()[-10:])
        raise GOAnnotateError(f"diamond getseq failed on {db} (exit {proc.returncode}):\n{tail}")
    return n


def build_reference_proteome(swissprot_db: str, taxid: str, output: Path, diamond_bin: str = None,
                             trembl_db: str = None, min_proteins: int = 10000) -> Path:
    """
    Write the proteins of the reference taxon to output, extracted with
    ``diamond getseq``: the SwissProt entries, and, when the taxon has fewer
    than min_proteins of them and a TrEMBL database is given, its TrEMBL
    entries as well (Drosophila melanogaster has about 3,900 SwissProt
    entries for about 14,000 protein-coding genes). An existing non-empty
    output is reused.
    """
    if output.exists() and output.stat().st_size > 0:
        logger.info(f"Reusing reference proteome: {output}")
        return output
    if not swissprot_db or not Path(swissprot_db).exists():
        raise GOAnnotateError("Building the reference proteome needs --diamond-swissprot "
                              "(or pass --reference-proteome FILE)")
    diamond_bin = find_diamond(diamond_bin)
    logger.info(f"Extracting taxon {taxid} proteins from {swissprot_db}")
    tmp = output.with_suffix(output.suffix + ".tmp")
    with open(tmp, "w") as fout:
        n = _getseq_taxon(diamond_bin, swissprot_db, taxid, fout)
        logger.info(f"  {n:,} SwissProt proteins")
        if n < min_proteins and trembl_db and Path(trembl_db).exists():
            logger.info(f"  Fewer than {min_proteins:,} SwissProt proteins; adding the TrEMBL entries of taxon "
                        f"{taxid} from {trembl_db}")
            m = _getseq_taxon(diamond_bin, trembl_db, taxid, fout)
            logger.info(f"  {m:,} TrEMBL proteins")
            n += m
    if n == 0:
        tmp.unlink()
        raise GOAnnotateError(f"No proteins with OX={taxid} in {swissprot_db}; "
                              f"pass --reference-proteome FILE")
    os.replace(tmp, output)
    logger.info(f"  {n:,} reference proteins written to {output}")
    return output


def read_reference_symbols(proteome: Path) -> tuple:
    """
    Parse reference proteome headers. Returns (accession -> symbol,
    upper-case symbol -> protein description) for proteins with a GN= symbol.
    """
    symbols = {}
    descriptions = {}
    with open(proteome) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            sid = header.split()[0]
            parts = sid.split("|")
            acc = parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else sid
            m = re.search(r"GN=(\S+)", header)
            if not m:
                continue
            symbols[acc] = gene_symbol_of(m.group(1))
            desc = header[len(sid):].strip()
            desc = desc.split(" OS=")[0]
            desc = re.sub(r"\s*(GN|PE|SV)=\S+", "", desc).strip()
            descriptions.setdefault(m.group(1).upper(), desc)
    logger.info(f"  {len(symbols):,} reference proteins carry a gene symbol")
    return symbols, descriptions


def read_reference_existence(proteome: Path) -> dict:
    """
    Upper-case symbol -> best UniProt protein existence level (PE=1 evidence
    at protein level ... 5 uncertain) among the reference proteins carrying
    it. Used to prefer an established gene over a retrocopy or a recent
    duplicate with weaker evidence that scores almost alike (NANOG PE=1 over
    NANOGP8 PE=2, CD8B over CD8B2 PE=3).
    """
    existence = {}
    with open(proteome) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            g = re.search(r"GN=(\S+)", line)
            pe = re.search(r"PE=(\d)", line)
            if g and pe:
                key = gene_symbol_of(g.group(1)).upper()
                existence[key] = min(existence.get(key, 9), int(pe.group(1)))
    return existence


####################
# SEARCHES
####################

def _makedb(diamond_bin, fasta: Path, db: Path, threads: int) -> None:
    if db.exists() and db.stat().st_mtime >= fasta.stat().st_mtime:
        return
    cmd = [diamond_bin, "makedb", "--in", str(fasta), "--db", str(db), "--threads", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.strip().splitlines()[-10:])
        raise GOAnnotateError(f"diamond makedb failed for {fasta}:\n{tail}")


def _blastp(diamond_bin, query: Path, db: Path, out: Path, threads: int, evalue: float,
            sensitivity_flags: list, force: bool) -> Path:
    if out.exists() and out.stat().st_size > 0 and not force:
        with open(out) as fh:
            n_cols = len(fh.readline().split("\t"))
        if n_cols >= 12:
            logger.info(f"  Reusing {out}")
            return out
        logger.info(f"  {out.name} lacks the length or alignment coordinate columns; rerunning")
    tmp = out.with_suffix(out.suffix + ".tmp")
    cmd = [diamond_bin, "blastp", "--db", str(db), "--query", str(query), "--out", str(tmp),
           "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qlen", "slen",
           "qstart", "qend", "sstart", "send",
           "--max-target-seqs", "25", "--evalue", str(evalue), "--threads", str(threads)] + sensitivity_flags
    start = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.strip().splitlines()[-10:])
        raise GOAnnotateError(f"diamond blastp failed ({query.name} vs {db.name}):\n{tail}")
    os.replace(tmp, out)
    logger.info(f"  {query.name} vs {db.name}: {time.time() - start:.0f}s")
    return out


def run_reciprocal_searches(query_proteins: Path, reference_proteome: Path, output_dir: Path, prefix: str,
                            threads: int, evalue: float, sensitivity_flags: list,
                            diamond_bin: str = None, force: bool = False) -> tuple:
    """Search queries against the reference proteome and the reverse. Returns (forward, reverse) paths."""
    diamond_bin = find_diamond(diamond_bin)
    output_dir = Path(output_dir)
    query_db = output_dir / f"{prefix}_proteins.dmnd"
    ref_db = output_dir / f"{prefix}_reference_proteome.dmnd"
    logger.info("Reciprocal searches against the reference proteome")
    _makedb(diamond_bin, query_proteins, query_db, threads)
    _makedb(diamond_bin, reference_proteome, ref_db, threads)
    forward = _blastp(diamond_bin, query_proteins, ref_db, output_dir / f"{prefix}_vs_reference.tsv",
                      threads, evalue, sensitivity_flags, force)
    reverse = _blastp(diamond_bin, reference_proteome, query_db, output_dir / f"{prefix}_reference_vs_query.tsv",
                      threads, evalue, sensitivity_flags, force)
    return forward, reverse


def best_hits(path: Path) -> dict:
    """query -> (subject, bitscore) of the highest-scoring hit."""
    best = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            q, s, bits = parts[0], parts[1], float(parts[5])
            if q not in best or bits > best[q][1]:
                best[q] = (s, bits)
    return best


def reference_alignments(forward: Path, ref_symbols: dict) -> dict:
    """
    transcript -> {upper-case reference symbol: (bits, sstart, send, slen)}
    of the highest-scoring alignment to each reference symbol's proteins.
    """
    out = defaultdict(dict)
    with open(forward) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            symbol = ref_symbols.get(_accession(parts[1]))
            if not symbol:
                continue
            key = symbol.upper()
            bits = float(parts[5])
            if bits > out[parts[0]].get(key, (-1.0,))[0]:
                out[parts[0]][key] = (bits, int(parts[10]), int(parts[11]), int(parts[7]))
    return out


def reference_scores_by_gene(forward: Path, ref_symbols: dict, tx2gene: dict) -> dict:
    """gene -> {upper-case reference symbol: best bitscore} over the gene's transcripts."""
    out = defaultdict(dict)
    for tx, aligns in reference_alignments(forward, ref_symbols).items():
        g = tx2gene.get(tx, tx)
        for sym, rec in aligns.items():
            if rec[0] > out[g].get(sym, 0.0):
                out[g][sym] = rec[0]
    return dict(out)


def _accession(subject_id: str) -> str:
    parts = subject_id.split("|")
    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
        return parts[1]
    return subject_id


####################
# ASSIGNMENT
####################

def assign_symbols(annotations: dict, universe, forward: Path, reverse: Path, ref_symbols: dict,
                   symbol_case_fn, ref_descriptions: dict = None,
                   support_fraction: float = 0.9, min_coverage: float = 0.5,
                   symbol_policy: str = "ortholog", nomenclature=None,
                   drop_unmatched: bool = False, min_query_coverage: float = 0.8,
                   rbh_margin: float = 0.9, existence: dict = None, established_margin: float = 0.95,
                   name_threshold: float = 0.5, partial_genes: set = None,
                   partial_query_coverage: float = 0.6) -> dict:
    """
    Apply reciprocal-best-hit confirmation and the sequence-support check
    to the annotations in place. Sets annot.symbol_status on every gene with
    hits and returns summary counts. Several genes may share a symbol.

    ref_descriptions: upper-case symbol -> reference protein description,
    used as the product name whenever the final symbol is a reference
    symbol supported by the gene's sequence.
    support_fraction: a voted symbol is supported when the gene's hit to
    that symbol's reference protein scores at least this fraction of the
    gene's best reference hit,
    min_coverage: and the alignment covers at least this fraction of the
    reference protein, or min_query_coverage of the gene itself (partial
    gene models keep the symbol of the protein they are part of).
    partial_genes: genes whose proteins are fragments (translated
    transcripts without both ends); for them partial_query_coverage of the
    gene replaces min_query_coverage, since a fragment's protein runs to
    the transcript ends and its ends rarely align.
    symbol_policy, for voted symbols outside the reference nomenclature:
    "ortholog" replaces one by a reference symbol only with evidence that
    it is the same gene (NCBI synonym, reciprocal best hit; see
    ortholog_evidence) and keeps the rest, or removes them with
    drop_unmatched; "reference" replaces it by the closest reference
    protein the gene covers (at least min_coverage of it), or removes it;
    "any" keeps it.
    nomenclature: mappings.ReferenceNomenclature (official symbols and
    synonyms of the reference species) for the synonym evidence.
    existence: upper-case symbol -> UniProt protein existence level of its
    reference proteins (read_reference_existence). After the checks, a
    gene's reference symbol whose protein existence is weaker than PE=1 is
    replaced by an established (PE=1) reference symbol the gene covers and
    hits at least established_margin as well (NANOGP8 -> NANOG, CD8B2 ->
    CD8B), when the weak symbol is named after it (begins with its symbol);
    the replaced symbol is kept in voted_symbol.
    Reference pseudogenes and readthrough loci (nomenclature.is_excluded)
    never name a gene.
    Product: when the reference symbol confirms the gene's own symbol (rbh,
    rbh_added, vote_supported), the reference
    protein's name replaces the voted name unless SwissProt hits of other
    species carry the voted name and none carries the reference name (names
    matching at name_threshold): a reference species' own gene name
    ("minidiscs", "polychaetoid" for a Drosophila reference) then leaves the
    hits' name; a curated name is otherwise kept (confirm_name). When
    the reference check changes the symbol to another gene, its reference
    name is used.
    """
    ref_descriptions = ref_descriptions or {}
    from goannotate.annotate import name_similarity, normalize_gene_name
    name_stats = Counter()

    def confirm_name(annot, symbol):
        """
        Reference name for a confirmed symbol. It replaces the voted name unless curated
        (SwissProt) entries of other species than the reference name the protein as the hits do
        and none of them uses the reference name: the reference name is then the reference
        species' own gene name (minidiscs, polychaetoid with a Drosophila reference, where other
        species say Y+L amino acid transporter 2, Tight junction protein ZO-2), and the hits'
        name stays. A curated name is kept wherever nothing shows that (Eukaryotic initiation
        factor 4A-III over the automatic family name RNA helicase).
        """
        desc = ref_descriptions.get(symbol)
        if not desc:
            return
        if not annot.consensus_name or desc == annot.consensus_name:
            annot.consensus_name = desc
            return
        curated = getattr(annot, "hit_curated_names", {}) or {}
        if curated:
            ref_norm = normalize_gene_name(desc)
            vote_norm = normalize_gene_name(annot.consensus_name)
            own = symbol.upper()

            def other_reference_gene(syms):
                # an entry of another reference gene (a paralog: UCP2 for UCP3) does not speak for this one
                return any(s in reference_symbols and s != own for s in syms)

            if (not any(name_similarity(ref_norm, h) >= name_threshold for h in curated)
                    and any(name_similarity(vote_norm, h) >= name_threshold and not other_reference_gene(syms)
                            for h, syms in curated.items())):
                name_stats["reference_name_not_supported"] += 1
                return
        annot.consensus_name = desc
        name_stats["reference_name_used"] += 1

    rev = best_hits(reverse)       # ref protein -> (query transcript, bits)
    tx2gene = universe.transcript_to_gene

    # Forward hits per gene: gene -> {reference symbol: best bits}, and the
    # best coverage of the reference protein reached by any of the gene's
    # transcripts: gene -> {reference symbol: aligned length / reference length}
    gene_forward_all = defaultdict(dict)
    gene_ref_coverage = defaultdict(dict)
    gene_query_coverage = defaultdict(dict)   # aligned length / query length
    with open(forward) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            symbol = ref_symbols.get(_accession(parts[1]))
            if not symbol:
                continue
            gene = tx2gene.get(parts[0], parts[0])
            key = symbol.upper()
            bits = float(parts[5])
            if bits > gene_forward_all[gene].get(key, -1.0):
                gene_forward_all[gene][key] = bits
            cov = int(parts[3]) / int(parts[7]) if len(parts) >= 8 and int(parts[7]) > 0 else 1.0
            if cov > gene_ref_coverage[gene].get(key, 0.0):
                gene_ref_coverage[gene][key] = cov
            qcov = int(parts[3]) / int(parts[6]) if len(parts) >= 8 and int(parts[6]) > 0 else 0.0
            if qcov > gene_query_coverage[gene].get(key, 0.0):
                gene_query_coverage[gene][key] = qcov
    # Gene-level forward best: gene -> (reference symbol, bits)
    gene_forward = {g: max(d.items(), key=lambda kv: kv[1]) for g, d in gene_forward_all.items()}

    # Reverse best per reference symbol: symbol -> (gene, bits), best over the symbol's proteins
    symbol_reverse = {}
    for ref_sid, (tx, bits) in rev.items():
        symbol = ref_symbols.get(_accession(ref_sid))
        if not symbol:
            continue
        gene = tx2gene.get(tx, tx)
        key = symbol.upper()
        if key not in symbol_reverse or bits > symbol_reverse[key][1]:
            symbol_reverse[key] = (gene, bits)

    # Reciprocal best hit: the gene is the best query gene of reference
    # symbol X, and X is among the gene's forward hits; the strongest such
    # X is taken. Requiring X to be the gene's single best forward hit
    # would lose paralog pairs where both members hit the better-known
    # sibling first (CALM1 over CALM2).
    # Human pseudogenes (NCBI gene type 'pseudo': CMAHP, VN1R17P, SSPOP) are
    # never given to a gene by the checks below; a functional gene whose human
    # counterpart is a pseudogene keeps its own symbol (CMAH).
    def is_pseudo(symbol):
        if not nomenclature:
            return False
        if hasattr(nomenclature, "is_excluded"):
            return nomenclature.is_excluded(symbol)
        return nomenclature.is_pseudogene(symbol)

    rbh = {}
    for gene, hits in gene_forward_all.items():
        reciprocal = [(bits, sym) for sym, bits in hits.items()
                      if symbol_reverse.get(sym, (None,))[0] == gene and not is_pseudo(sym)]
        if reciprocal:
            rbh[gene] = max(reciprocal)[1]

    stats = Counter()
    reference_symbols = {s.upper() for s in ref_symbols.values()}
    # Step 1 and 2: confirmation and RBH override.
    for gene, annot in annotations.items():
        voted = annot.consensus_symbol.upper() if annot.consensus_symbol else ""
        rbh_symbol = rbh.get(gene)
        if not annot.hits and not rbh_symbol:
            annot.symbol_status = ""
            continue
        if not voted:
            annot.symbol_status = "none" if annot.hits else ""
            if rbh_symbol:
                # a reciprocal best hit names a gene the vote left unnamed
                annot.consensus_symbol = symbol_case_fn(rbh_symbol)
                annot.symbol_status = "rbh_added"
                confirm_name(annot, rbh_symbol)
            continue
        if rbh_symbol is None:
            if annot.symbol_status != "propagated":
                annot.symbol_status = "vote"
        elif rbh_symbol == voted:
            annot.symbol_status = "rbh"
            # The symbol is confirmed by orthology; the reference protein's
            # curated name replaces the voted name, which often comes from an
            # automatically named TrEMBL entry ("Serine/threonine-protein"),
            # unless it is the reference species' own gene name (confirm_name).
            confirm_name(annot, rbh_symbol)
        else:
            cluster_symbols = {s.upper() for s in annot.winning_cluster_symbols if s}
            if rbh_symbol in cluster_symbols:
                annot.consensus_symbol = symbol_case_fn(rbh_symbol)
                annot.symbol_status = "rbh_override"
                if ref_descriptions.get(rbh_symbol):
                    annot.consensus_name = ref_descriptions[rbh_symbol]
            else:
                annot.symbol_status = "vote_rbh_conflict"
    for annot in annotations.values():
        if annot.symbol_status:
            stats[annot.symbol_status] += 1

    # Step 3: sequence support for voted symbols. Several genes may carry the
    # same symbol (paralogs, recent duplicates, retrocopies); that is not
    # checked or changed. What is checked is whether each gene's own
    # sequence supports its voted symbol: the reference protein of that
    # symbol must be among the gene's hits to the reference proteome and
    # score within support_fraction of the gene's best reference hit. A
    # gene whose vote is not supported (typically a gene that shares a
    # domain with a well-studied protein and was named after it) takes the
    # symbol of its closest reference protein instead. Genes without any
    # hit to the reference proteome (lineage-specific families) keep the
    # vote, since nothing in the reference contradicts it.

    partial_genes = partial_genes or set()

    def covers(gene, symbol):
        """
        The gene's alignment to the reference protein covers at least
        min_coverage of the reference protein, or at least
        min_query_coverage of the gene itself (partial_query_coverage for a
        fragment): a partial gene model or a shorter protein is the same
        gene; a gene that shares one domain covers little of either.
        """
        need = partial_query_coverage if gene in partial_genes else min_query_coverage
        return (gene_ref_coverage.get(gene, {}).get(symbol, 0.0) >= min_coverage
                or gene_query_coverage.get(gene, {}).get(symbol, 0.0) >= need)

    def is_supported(gene, symbol, best_bits):
        bits = gene_forward_all.get(gene, {}).get(symbol, 0.0)
        return bits >= support_fraction * best_bits and covers(gene, symbol)

    nomenclature = nomenclature or {}

    def ortholog_evidence(gene, voted):
        """
        (reference symbol, kind of evidence) when a voted non-reference
        symbol names the same gene as a reference protein the gene's
        sequence covers (at least min_coverage of it):
          renamed  the symbol is a unique NCBI synonym of that reference
                   symbol (EIF2C3 -> AGO3)
          rbh      the gene is the reciprocal best hit of that reference symbol
        Comparisons are case-insensitive. ('', '') without evidence.
        """
        current = nomenclature.current(voted).upper() if nomenclature else ""
        if current and current != voted and current in reference_symbols and covers(gene, current) \
                and not is_pseudo(current):
            return current, "renamed"
        r = rbh.get(gene)
        if r and covers(gene, r) and unambiguous(gene, r):
            return r, "rbh"
        return "", ""

    def unambiguous(gene, symbol):
        """
        The reference protein is the gene's best reference hit and the
        second-best reference protein scores below rbh_margin of it. In
        families expanded in one lineage (olfactory and vomeronasal
        receptors, MHC, cytochromes P450, globins) several reference
        paralogs score almost alike, and a reciprocal best hit among them
        is not evidence of orthology.
        """
        if rbh_margin <= 0:
            return True
        hits = sorted(gene_forward_all.get(gene, {}).items(), key=lambda kv: -kv[1])
        if not hits or hits[0][0] != symbol:
            return False
        return len(hits) == 1 or hits[1][1] < rbh_margin * hits[0][1]

    replaced = 0
    supported = 0
    unverified = 0
    no_reference = 0
    foreign_kept = 0
    foreign_replaced = 0
    foreign_removed = 0
    foreign_ortholog = Counter()
    removed = 0
    for gene, annot in annotations.items():
        if annot.symbol_status not in ("vote", "propagated", "vote_rbh_conflict"):
            continue
        voted = annot.consensus_symbol.upper()
        if voted not in reference_symbols:
            # A symbol outside the reference nomenclature (a zebrafish
            # foxo3a, a rodent Vmn1r, an outdated AARS).
            #   any:       kept.
            #   ortholog:  replaced by a reference symbol only with evidence
            #              that it is the same gene (see ortholog_evidence);
            #              otherwise kept (foreign_kept), or removed with
            #              drop_unmatched.
            #   reference: the closest reference protein the gene covers,
            #              or no symbol.
            # A replaced symbol is kept in voted_symbol.
            if symbol_policy == "any":
                foreign_kept += 1
                continue
            if symbol_policy == "ortholog":
                target, how = ortholog_evidence(gene, voted)
                if target:
                    annot.voted_symbol = annot.consensus_symbol
                    annot.consensus_symbol = symbol_case_fn(target)
                    annot.symbol_status = "foreign_ortholog"
                    if ref_descriptions.get(target):
                        annot.consensus_name = ref_descriptions[target]
                    foreign_ortholog[how] += 1
                elif drop_unmatched:
                    annot.voted_symbol = annot.consensus_symbol
                    annot.consensus_symbol = ""
                    annot.symbol_status = "foreign_removed"
                    foreign_removed += 1
                else:
                    annot.symbol_status = "foreign_kept"
                    foreign_kept += 1
                continue
            fwd_hits = gene_forward_all.get(gene) or {}
            best_bits = max(fwd_hits.values()) if fwd_hits else 0.0
            candidates = [(bits, sym) for sym, bits in fwd_hits.items()
                      if is_supported(gene, sym, best_bits) and not is_pseudo(sym)]
            annot.voted_symbol = annot.consensus_symbol
            if candidates:
                best_symbol = max(candidates)[1]
                annot.consensus_symbol = symbol_case_fn(best_symbol)
                annot.symbol_status = "foreign_replaced"
                if ref_descriptions.get(best_symbol):
                    annot.consensus_name = ref_descriptions[best_symbol]
                foreign_replaced += 1
            else:
                annot.consensus_symbol = ""
                annot.symbol_status = "foreign_removed"
                foreign_removed += 1
            continue
        fwd_hits = gene_forward_all.get(gene)
        if not fwd_hits:
            no_reference += 1
            continue
        best_bits = max(fwd_hits.values())
        if is_supported(gene, voted, best_bits):
            supported += 1
            annot.symbol_status = annot.symbol_status if annot.symbol_status != "vote" else "vote_supported"
            confirm_name(annot, voted)
            continue
        annot.voted_symbol = annot.consensus_symbol
        candidates = [(bits, sym) for sym, bits in fwd_hits.items()
                      if is_supported(gene, sym, best_bits) and not is_pseudo(sym)]
        if candidates:
            best_symbol = max(candidates)[1]
            annot.consensus_symbol = symbol_case_fn(best_symbol)
            annot.symbol_status = "reference_best_hit"
            if ref_descriptions.get(best_symbol):
                annot.consensus_name = ref_descriptions[best_symbol]
            replaced += 1
        elif gene_forward_all.get(gene, {}).get(voted, 0.0) >= support_fraction * best_bits:
            # No reference protein covers the gene (a shared domain, a
            # fragment, a lineage-specific gene) but the voted symbol's
            # protein is the gene's best or near-best reference hit: the
            # reference agrees on the identity and cannot confirm the
            # coverage, so the vote is kept unverified with lower confidence.
            # The product name from the vote is kept.
            annot.consensus_symbol = annot.voted_symbol
            annot.voted_symbol = ""
            annot.symbol_status = "vote_unverified"
            unverified += 1
        else:
            # The gene's sequence prefers other reference proteins to the
            # voted symbol's (or does not hit it at all) and none of them
            # covers it: the symbol is removed, the product name kept.
            annot.consensus_symbol = ""
            annot.symbol_status = "unsupported_removed"
            removed += 1

    # Established genes over weakly evidenced near-identical paralogs (retrocopies, recent duplicates)
    existence = existence or {}
    established = 0
    if existence:
        for gene, annot in annotations.items():
            sym = annot.consensus_symbol.upper() if annot.consensus_symbol else ""
            if not sym or existence.get(sym, 1) <= 1:
                continue
            hits = gene_forward_all.get(gene) or {}
            bits = hits.get(sym, 0.0)
            if bits <= 0:
                continue
            # only the gene the weak symbol is named after (NANOG for NANOGP8, CD8B for CD8B2,
            # PRAME for PRAMEF14), not another member of a family (OR4F6 and OR4F3)
            better = [(b, s) for s, b in hits.items() if s != sym and existence.get(s, 9) == 1
                      and sym.startswith(s) and b >= established_margin * bits and covers(gene, s)
                      and not is_pseudo(s)]
            if not better:
                continue
            target = max(better)[1]
            if not annot.voted_symbol:
                annot.voted_symbol = annot.consensus_symbol
            if annot.consensus_symbol not in annot.symbol_aliases:
                annot.symbol_aliases = list(annot.symbol_aliases) + [annot.consensus_symbol]
            annot.consensus_symbol = symbol_case_fn(target)
            if ref_descriptions.get(target):
                annot.consensus_name = ref_descriptions[target]
            established += 1
    stats["established_preferred"] = established
    stats.update(name_stats)
    if name_stats:
        logger.info(f"  Reference names for confirmed symbols: {name_stats['reference_name_used']:,} used, "
                    f"{name_stats['reference_name_not_supported']:,} not supported by the gene's hits "
                    f"(the hits' name kept)")

    stats["vote_supported"] = supported
    stats["vote_no_reference_hit"] = no_reference
    stats["vote_foreign_symbol"] = foreign_kept
    stats["foreign_replaced"] = foreign_replaced
    stats["foreign_removed"] = foreign_removed
    stats["foreign_ortholog"] = sum(foreign_ortholog.values())
    for how in ("renamed", "rbh"):
        stats[f"foreign_ortholog_{how}"] = foreign_ortholog[how]
    stats["reference_best_hit"] = replaced
    stats["vote_unverified"] = unverified
    stats["unsupported_removed"] = removed
    stats["rbh_total"] = len(rbh)
    stats["symbol_policy"] = symbol_policy

    symbols = [a.consensus_symbol.upper() for a in annotations.values() if a.consensus_symbol]
    copies = Counter(symbols)
    shared = sum(1 for n in copies.values() if n > 1)
    logger.info(f"Symbol assignment: {len(rbh):,} reciprocal best hits; "
                f"confirmed={stats['rbh']:,}, rbh override={stats['rbh_override']:,}, "
                f"rbh added={stats['rbh_added']:,}")
    logger.info(f"  Voted reference symbols checked against the reference proteome: {supported:,} supported, "
                f"{replaced:,} replaced by the closest supported reference protein, "
                f"{removed:,} removed (no reference protein covered), "
                f"{no_reference:,} without a reference hit (kept)")
    if symbol_policy == "any":
        logger.info(f"  {foreign_kept:,} symbols outside the reference nomenclature kept (--symbol-policy any)")
    elif symbol_policy == "ortholog":
        logger.info(f"  Symbols outside the reference nomenclature: {sum(foreign_ortholog.values()):,} replaced by "
                    f"the reference ortholog ({foreign_ortholog['renamed']:,} "
                    f"by NCBI synonym, {foreign_ortholog['rbh']:,} by reciprocal best hit), {foreign_kept:,} kept, "
                    f"{foreign_removed:,} removed (--symbol-policy ortholog)")
    else:
        logger.info(f"  Symbols outside the reference nomenclature: {foreign_replaced:,} replaced by the closest "
                    f"covered reference protein, {foreign_removed:,} removed (--symbol-policy reference)")
    if established:
        logger.info(f"  {established:,} symbols with weak protein evidence (PE>1) replaced by an established "
                    f"reference gene scoring at least {established_margin:.2f} as well (NANOGP8 -> NANOG)")
    logger.info(f"  {len(symbols):,} symbols on {len(copies):,} distinct values; "
                f"{shared:,} symbols are carried by more than one gene")
    return dict(stats)


####################
# SECONDARY NOMENCLATURE
####################

def secondary_symbols(annotations: dict, universe, forward: Path, reverse: Path, ref_symbols: dict,
                      symbol_case_fn, nomenclature=None, min_coverage: float = 0.5,
                      min_query_coverage: float = 0.8, rbh_margin: float = 0.9) -> dict:
    """
    Symbol of each gene in a second nomenclature (annot.secondary_symbol),
    in place: the gene must be the reciprocal best hit of a protein of the
    secondary reference proteome, that protein must be the gene's best hit
    in it, the second-best protein must score below rbh_margin of it, and
    the alignment must cover min_coverage of the protein or
    min_query_coverage of the gene. Pseudogenes are not used. These are the
    conditions under which a reciprocal best hit counts as ortholog
    evidence in assign_symbols. Returns counts.
    """
    tx2gene = universe.transcript_to_gene
    fwd = defaultdict(dict)
    cov = defaultdict(dict)
    with open(forward) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            symbol = ref_symbols.get(_accession(parts[1]))
            if not symbol:
                continue
            gene = tx2gene.get(parts[0], parts[0])
            key = symbol.upper()
            bits = float(parts[5])
            if bits > fwd[gene].get(key, -1.0):
                fwd[gene][key] = bits
            length, qlen, slen = int(parts[3]), int(parts[6]), int(parts[7])
            covered = ((slen > 0 and length / slen >= min_coverage)
                       or (qlen > 0 and length / qlen >= min_query_coverage))
            cov[gene][key] = cov[gene].get(key, False) or covered
    best_query = {}
    for ref_sid, (tx, bits) in best_hits(reverse).items():
        symbol = ref_symbols.get(_accession(ref_sid))
        if not symbol:
            continue
        key = symbol.upper()
        if key not in best_query or bits > best_query[key][1]:
            best_query[key] = (tx2gene.get(tx, tx), bits)
    n = 0
    for gene, annot in annotations.items():
        annot.secondary_symbol = ""
        hits = sorted(fwd.get(gene, {}).items(), key=lambda kv: -kv[1])
        if not hits:
            continue
        top, bits = hits[0]
        if best_query.get(top, (None,))[0] != gene or not cov[gene].get(top):
            continue
        if nomenclature and nomenclature.is_pseudogene(top):
            continue
        if rbh_margin > 0 and len(hits) > 1 and hits[1][1] >= rbh_margin * bits:
            continue
        annot.secondary_symbol = symbol_case_fn(top)
        n += 1
    logger.info(f"Secondary symbols: {n:,} genes have an unambiguous reciprocal best hit in the secondary "
                f"reference proteome")
    return {"secondary_symbols": n}


####################
# SYMBOL ALIASES
####################

# Statuses whose final symbol is the gene the vote chose, so the GN= symbols
# of the hits behind the vote name the same gene
_VOTE_KEPT = {"vote", "vote_supported", "rbh", "propagated", "foreign_kept", "foreign_ortholog",
              "vote_rbh_conflict"}


def secondary_names(annotations: dict, secondary_descriptions: dict, ref_descriptions: dict,
                    name_threshold: float = 0.5) -> dict:
    """
    Products from the secondary nomenclature, in place: a gene whose product is
    its primary reference protein's name, where none of the gene's hits uses
    that name (a reference species' own gene name: larval translucida, LD43551p
    without a curated name of another species to check it against), takes the
    name of its secondary reference protein (annot.secondary_symbol, an
    unambiguous reciprocal best hit) when at least one of its hits uses that
    name. Returns counts.
    """
    from goannotate.annotate import name_similarity, normalize_gene_name
    n = 0
    examples = []
    for annot in annotations.values():
        if not (annot.secondary_symbol and annot.consensus_symbol and annot.hit_names):
            continue
        primary = ref_descriptions.get(annot.consensus_symbol.upper())
        if not primary or annot.consensus_name != primary:
            continue
        secondary = secondary_descriptions.get(annot.secondary_symbol.upper())
        if not secondary or secondary == primary:
            continue
        p_norm, s_norm = normalize_gene_name(primary), normalize_gene_name(secondary)
        if any(name_similarity(p_norm, h) >= name_threshold for h in annot.hit_names):
            continue
        if not any(name_similarity(s_norm, h) >= name_threshold for h in annot.hit_names):
            continue
        annot.consensus_name = secondary
        annot.product_note = "secondary_reference"
        n += 1
        if len(examples) < 5:
            examples.append(f"{primary} -> {secondary}")
    if n:
        logger.info(f"  Products from the secondary reference: {n:,} genes whose primary reference name no hit "
                    f"uses (e.g. {'; '.join(examples[:3])})")
    return {"secondary_names": n}


def _alias_key(symbol: str) -> str:
    """Spelling-insensitive key: case and hyphens ignored (IL-15 = il15)."""
    return symbol.lower().replace("-", "")


def add_symbol_aliases(annotations: dict, nomenclature, is_bad_symbol, symbol_case_fn, max_aliases: int = 25) -> dict:
    """
    Other symbols of each gene's final symbol, in place (annot.symbol_aliases):
    the GN= symbols of the hits in the winning symbol group (other species'
    names of the same gene: Foxo3, foxo3a, foxo3b), the voted symbol when a
    reference ortholog replaced it, and the reference species' NCBI synonyms
    of the final symbol. Genome-specific identifiers (symbol patterns) are
    left out; duplicates are removed case-insensitively and every alias is
    written in the --symbol-case convention. Hit symbols are
    used only when the final symbol is the one the vote chose (or its
    reference ortholog). Returns counts.
    """
    n_genes = n_aliases = 0
    for annot in annotations.values():
        final = annot.consensus_symbol
        if not final:
            annot.symbol_aliases = []
            continue
        seen = {_alias_key(final)}
        out = []

        def add(sym):
            if not sym or len(out) >= max_aliases:
                return
            key = _alias_key(sym)
            if key in seen or is_bad_symbol(sym):
                return
            seen.add(key)
            out.append(symbol_case_fn(sym))

        if annot.symbol_status == "foreign_ortholog":
            add(annot.voted_symbol)
        if annot.symbol_status in _VOTE_KEPT:
            for hit in annot.filtered_hits:
                add(hit.gene_symbol)
        if nomenclature:
            for syn in sorted(nomenclature.synonyms.get(final.lower(), ())):
                add(syn)
        annot.symbol_aliases = out
        if out:
            n_genes += 1
            n_aliases += len(out)
    logger.info(f"Symbol aliases: {n_aliases:,} on {n_genes:,} genes (hit symbols of the same gene, replaced "
                f"symbols, reference-species synonyms)")
    return {"genes_with_aliases": n_genes}


####################
# PARALOG DESIGNATORS
####################

# Designators that another genome's nomenclature adds to a name for its own
# duplicates: a trailing " a"/" b" word, a letter glued to a trailing number
# ("O3b", "11ba"), "-A"/"-B" after a number ("O3-B"). Case-sensitive.
_DESIGNATOR_FORMS = (
    re.compile(r"^(.*\S) [ab]$"),
    re.compile(r"^(.*\d)[ab]{1,2}$"),
    re.compile(r"^(.*\d)-[AB]$"),
)


def strip_paralog_designators(annotations: dict, ref_descriptions: dict, normalize_fn) -> dict:
    """
    Remove paralog designators of other genomes from voted product names,
    in place. A designator is removed only when the name without it is a
    reference protein name or a name carried by the gene's own hits
    (winning cluster), so that "Cytochrome b" or "Annexin A" are left
    alone. Products taken from the reference proteome are not touched.
    Returns counts.
    """
    ref_names = {normalize_fn(d) for d in ref_descriptions.values() if d}
    from_reference = {"rbh", "rbh_override", "rbh_added", "vote_supported", "reference_best_hit",
                      "foreign_replaced"}
    stripped = 0
    examples = Counter()
    for annot in annotations.values():
        name = annot.consensus_name
        if not name or (annot.symbol_status in from_reference and normalize_fn(name) in ref_names):
            continue        # a product taken from the reference proteome
        pool = {normalize_fn(n) for n in annot.winning_cluster_names if n}
        for form in _DESIGNATOR_FORMS:
            m = form.match(name)
            if not m:
                continue
            base = m.group(1).rstrip(" ,-")
            norm = normalize_fn(base)
            if norm in ref_names or (norm in pool and norm != normalize_fn(name)):
                examples[f"{name} -> {base}"] += 1
                annot.consensus_name = base
                annot.product_note = "designator_stripped"
                stripped += 1
                break
    if stripped:
        logger.info(f"  Paralog designators removed from {stripped:,} product names, e.g. "
                    + "; ".join(k for k, _ in examples.most_common(5)))
    return {"designators_stripped": stripped}
