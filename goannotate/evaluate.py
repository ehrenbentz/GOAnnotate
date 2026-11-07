"""
Evaluation of the final annotation against each source of evidence.

For every gene with hits or a symbol, the final symbol and product are
scored against each source independently: the gene's own hits, the
reference proteome, the relatives' orthologs, the same-species annotation,
the InterPro entries, and the product against the symbol. Nothing is
changed; the report lists each gene's evidence and the summary counts how
often each source agrees with the final call.
"""

import logging
from collections import Counter, defaultdict
from pathlib import Path

from goannotate.annotate import load_stop_words, normalize_gene_name, symbol_bases

logger = logging.getLogger(__name__)

# Reference evidence: a symbol's reference protein within this fraction of
# the gene's best reference hit counts as agreeing with the sequence
NEAR_BEST = 0.9
# Hit support at or above this fraction of the hits counts as agreement
HIT_AGREE = 0.5

SOURCES = ("hits_symbol", "hits_name", "reference", "relatives", "same_species", "domains", "product_symbol")

COLUMNS = ["gene_id", "symbol", "symbol_status", "orthology_decision", "decided_by_relatives", "product",
           "hits_symbol", "hits_symbol_fraction", "hits_name", "hits_name_fraction",
           "reference", "reference_best_symbol", "reference_ratio",
           "relatives", "relatives_symbol", "same_species", "same_species_symbol",
           "domains", "domain_names", "product_symbol", "product_generic", "n_agree", "n_disagree", "conflict"]


def _tokens(name: str, stop_words) -> set:
    """Informative words of a name: stop words removed, separators split."""
    words = set()
    for w in normalize_gene_name(name).replace("-", " ").replace("/", " ").replace(",", " ").split():
        w = w.strip("()[]")
        if len(w) >= 2 and w not in stop_words:
            words.add(w)
    return words


def _names_agree(a: str, b: str, stop_words) -> bool:
    """Two names agree when they share at least half of the informative words of the shorter one."""
    ta, tb = _tokens(a, stop_words), _tokens(b, stop_words)
    if not ta or not tb:
        return False
    return len(ta & tb) >= 0.5 * min(len(ta), len(tb))


def _reference_evidence(gene_ref: dict, symbol: str, reference_symbol: str, ref_symbol_set: set) -> tuple:
    """
    (class, best symbol, ratio) of the final symbol against the gene's
    alignments to the reference proteome: best, near_best, below, absent
    (a reference symbol the gene does not hit), not_reference (a symbol
    outside the reference nomenclature with no reference form), none (no
    reference hit).
    """
    if not gene_ref:
        return "none", "", ""
    best_sym, best = max(gene_ref.items(), key=lambda kv: kv[1])
    candidates = {symbol.upper()} | symbol_bases(symbol.upper())
    if reference_symbol:
        candidates.add(reference_symbol.upper())
    scores = [gene_ref[s] for s in candidates if s in gene_ref]
    if not scores:
        in_ref = any(s in ref_symbol_set for s in candidates)
        return ("absent" if in_ref else "not_reference"), best_sym, ""
    ratio = max(scores) / best if best else 0.0
    if ratio >= 1.0:
        return "best", best_sym, f"{ratio:.2f}"
    if ratio >= NEAR_BEST:
        return "near_best", best_sym, f"{ratio:.2f}"
    return "below", best_sym, f"{ratio:.2f}"


def evaluate(annotations: dict, universe, ref_alignments: dict, ref_symbol_set: set, ref_descriptions: dict,
             modules_run: list, report: Path, stop_words=None, bad_name_filter=None) -> dict:
    """
    Write the per-gene evaluation report and return the summary statistics:
    for each source, how many genes it agrees with, disagrees with or says
    nothing about; relatives split by whether they decided the symbol; the
    conflicts (genes with a symbol that some source contradicts).

    ref_alignments: transcript -> {reference symbol: (bits, ...)} from the
    forward search against the reference proteome.
    ref_symbol_set: the reference nomenclature's symbols, upper case.
    bad_name_filter: marks family-level products (product_generic), which
    name a family rather than the gene of the symbol.
    """
    stop_words = stop_words or load_stop_words()
    orthology = next((m for m in modules_run if hasattr(m, "ortholog_evidence_for")), None)
    interpro = next((m for m in modules_run if hasattr(m, "gene_entries") and hasattr(m, "entries")), None)
    # reference evidence per gene: the best alignment to each reference symbol over the gene's transcripts
    gene_ref = defaultdict(dict)
    tx2gene = universe.transcript_to_gene if universe is not None else {}
    for tx, aligns in (ref_alignments or {}).items():
        g = tx2gene.get(tx, tx)
        for sym, rec in aligns.items():
            bits = rec[0] if isinstance(rec, tuple) else rec
            if bits > gene_ref[g].get(sym, 0.0):
                gene_ref[g][sym] = bits

    counts = {s: Counter() for s in SOURCES}
    relatives_split = {"decided": Counter(), "independent": Counter()}
    conflicts = Counter()
    n_conflict = 0
    n_rows = 0
    Path(report).parent.mkdir(parents=True, exist_ok=True)
    with open(report, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for gene_id, annot in annotations.items():
            symbol = annot.consensus_symbol or ""
            if not symbol and not annot.hits:
                continue
            n_rows += 1
            decision = orthology.decision.get(gene_id, "") if orthology is not None else ""
            decided = (annot.symbol_source == "relative_override") or decision in ("taken", "refined", "filled")
            row = {c: "" for c in COLUMNS}
            row.update(gene_id=gene_id, symbol=symbol, symbol_status=annot.symbol_status or "",
                       orthology_decision=decision, decided_by_relatives="yes" if decided else "no",
                       product=annot.consensus_name or "")
            verdict = {}

            # hits
            if symbol and annot.hit_symbols:
                frac = annot.symbol_concordance
                verdict["hits_symbol"] = "agree" if frac >= HIT_AGREE else "minority" if frac > 0 else "disagree"
                row["hits_symbol_fraction"] = f"{frac:.2f}"
            else:
                verdict["hits_symbol"] = "none"
            if annot.consensus_name and annot.hit_names:
                frac = annot.name_concordance
                verdict["hits_name"] = "agree" if frac >= HIT_AGREE else "minority" if frac > 0 else "disagree"
                row["hits_name_fraction"] = f"{frac:.2f}"
            else:
                verdict["hits_name"] = "none"

            # reference proteome
            if symbol:
                cls, best_sym, ratio = _reference_evidence(gene_ref.get(gene_id, {}), symbol,
                                                           annot.ortholog_reference_symbol or "", ref_symbol_set)
                verdict["reference"] = cls
                row["reference_best_symbol"] = best_sym
                row["reference_ratio"] = ratio
            else:
                verdict["reference"] = "none"

            # relatives and the same species
            rel_sym = same_sym = ""
            verdict["relatives"] = verdict["same_species"] = "none"
            if orthology is not None:
                ev = orthology.ortholog_evidence_for(gene_id)
                rel_sym, same_sym = ev.get("relatives", ""), ev.get("same_species", "")
                if bad_name_filter is not None:      # an identifier (LOC123) names nothing
                    rel_sym = "" if rel_sym and bad_name_filter.is_bad_symbol(rel_sym) else rel_sym
                    same_sym = "" if same_sym and bad_name_filter.is_bad_symbol(same_sym) else same_sym
                from_relatives = (annot.ortholog_reference_symbol or "").upper() == rel_sym.upper() and rel_sym
                if symbol and rel_sym:
                    agree = from_relatives or orthology._equivalent(symbol, rel_sym)
                    verdict["relatives"] = "agree" if agree else "disagree"
                    relatives_split["decided" if decided else "independent"][verdict["relatives"]] += 1
                if symbol and same_sym:
                    verdict["same_species"] = "agree" if orthology._equivalent(symbol, same_sym) else "disagree"
            row["relatives_symbol"], row["same_species_symbol"] = rel_sym, same_sym

            # domains: the product shares an informative word with an InterPro entry of the gene.
            # Entry descriptions name domains and families, not genes, so a shared word supports
            # the product and the absence of one says nothing against it.
            verdict["domains"] = "none"
            if interpro is not None and annot.consensus_name:
                names = [interpro.entries.get(acc, ("", ""))[0] for _, acc in interpro.gene_entries.get(gene_id, [])]
                names = [n for n in names if n]
                if names:
                    words = _tokens(annot.consensus_name, stop_words)
                    shared = any(words & _tokens(n, stop_words) for n in names)
                    verdict["domains"] = "agree" if shared else "no_shared_word" if words else "none"
                    row["domain_names"] = "|".join(names[:5])

            # product against the symbol: the reference protein's name, else curated hits carrying the symbol
            verdict["product_symbol"] = "none"
            if symbol and annot.consensus_name:
                # the reference protein's name, and the curated hits' names of the gene (the vote's
                # products follow the hits' naming style, the reference's may differ)
                ref_name = ref_descriptions.get((annot.ortholog_reference_symbol or symbol).upper(), "")
                product = annot.consensus_name
                forms = {symbol.upper()} | symbol_bases(symbol.upper())
                if annot.ortholog_reference_symbol:
                    forms.add(annot.ortholog_reference_symbol.upper())
                carriers = [name for name, syms in (annot.hit_curated_names or {}).items()
                            if any(s.upper() in forms for s in syms)]
                names = ([ref_name] if ref_name else []) + carriers
                if names:
                    ok = any(_names_agree(product, name, stop_words) for name in names)
                    verdict["product_symbol"] = "agree" if ok else "disagree"
                if bad_name_filter is not None and bad_name_filter.is_generic_name(product):
                    row["product_generic"] = "yes"

            for s in SOURCES:
                counts[s][verdict[s]] += 1
                row[s] = verdict[s]
            disagree = [s for s in SOURCES if verdict[s] in ("disagree", "below", "absent")]
            agree = [s for s in SOURCES if verdict[s] in ("agree", "best", "near_best")]
            row["n_agree"], row["n_disagree"] = str(len(agree)), str(len(disagree))
            if symbol and disagree:
                n_conflict += 1
                conflicts[",".join(disagree)] += 1
                row["conflict"] = ",".join(disagree)
            fh.write("\t".join(row[c] for c in COLUMNS) + "\n")

    stats = {"rows": n_rows, "counts": counts, "relatives_split": relatives_split,
             "conflicts": n_conflict, "conflict_kinds": conflicts, "report": str(report)}
    logger.info(f"Evaluation of the final annotation: {n_rows:,} genes; {n_conflict:,} with a symbol that a source "
                f"contradicts; {report}")
    return stats


def summary_lines(stats: dict) -> list:
    """Lines for the summary file."""
    if not stats:
        return []
    c = stats["counts"]
    lines = ["Evaluation of the final annotation against each source (reports: evaluation.tsv):",
             "-" * 70,
             f"  Genes evaluated:              {stats['rows']:,} (with hits or a symbol)",
             "  Source                        agree    disagree   partial  no evidence",
             "  (partial: a minority of the hits, or no word shared with the InterPro entries)"]

    def line(label, cnt, agree_keys=("agree",), disagree_keys=("disagree",), minority_keys=("minority",)):
        agree = sum(cnt[k] for k in agree_keys)
        disagree = sum(cnt[k] for k in disagree_keys)
        minority = sum(cnt[k] for k in minority_keys)
        none = sum(v for k, v in cnt.items() if k not in agree_keys + disagree_keys + minority_keys)
        return f"  {label:<30}{agree:>7,}{disagree:>11,}{minority:>10,}{none:>13,}"

    lines.append(line("Hits: symbol", c["hits_symbol"]))
    lines.append(line("Hits: product name", c["hits_name"]))
    lines.append(line("Reference proteome", c["reference"], ("best", "near_best"), ("below", "absent"), ()))
    lines.append(line("Relatives' orthologs", c["relatives"]))
    rs = stats["relatives_split"]
    for key, label in (("independent", "  relatives did not decide"), ("decided", "  relatives decided")):
        if rs[key]:
            lines.append(line(label, rs[key]))
    lines.append(line("Same-species annotation", c["same_species"]))
    lines.append(line("InterPro entries (words)", c["domains"], ("agree",), (), ("no_shared_word",)))
    lines.append(line("Product names the symbol", c["product_symbol"]))
    lines.append(f"  Conflicts:                    {stats['conflicts']:,} genes with a symbol that a source contradicts "
                 f"(column conflict; product_generic marks family-level products)")
    for kind, n in stats["conflict_kinds"].most_common(6):
        lines.append(f"    {n:>7,}  {kind}")
    lines.append("")
    return lines
