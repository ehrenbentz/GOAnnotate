"""Pipeline orchestration: inputs -> sequences -> Diamond -> annotation -> comparison -> outputs."""

import logging
import os
import shutil
from pathlib import Path

from goannotate import GOAnnotateError
from goannotate import annotate as annotate_mod
from goannotate import (assign, blast, compare, gff3, mappings, modules, output, refgo, reference, search,
                        sequences, taxon)
from goannotate.filters import BadNameFilter, NameAudit, write_name_audit
from goannotate.go import GOHierarchy, resolve_go_obo

logger = logging.getLogger(__name__)


def _coverage_filter(args) -> blast.HitCoverageFilter:
    return blast.HitCoverageFilter(
        mode=getattr(args, "hit_coverage", "permissive"),
        min_query_cov=getattr(args, "min_query_cov", 0.5),
        min_subject_cov=getattr(args, "min_subject_cov", 0.5),
        min_identity=getattr(args, "min_identity", 40.0),
        min_align_length=getattr(args, "min_align_length", 100),
        fallback=not getattr(args, "no_coverage_fallback", False))


def _taxon_constraints(args, go_hierarchy) -> dict:
    """
    GO terms excluded for the query organism (GO ID -> reason), or {} when
    the constraints are off, no --taxid was given, or no constraint data
    is available (SQLite tables from the current build_databases.py, or the two
    OBO files).
    """
    if getattr(args, "taxon_constraints", "on") == "off":
        logger.info("GO taxon constraints disabled (--taxon-constraints off)")
        return {}
    taxid = getattr(args, "taxid", None)
    if not taxid:
        logger.info("GO taxon constraints not applied: give --taxid TAXID (the query organism's NCBI taxon) "
                    "to remove GO terms that do not apply to it")
        return {}
    constraints = unions = None
    if getattr(args, "taxon_constraints_obo", None):
        constraints = taxon.parse_taxon_constraints_obo(args.taxon_constraints_obo)
        unions = (taxon.parse_taxon_groupings_obo(args.taxon_groupings_obo)
                  if getattr(args, "taxon_groupings_obo", None) else {})
    elif args.db:
        constraints, unions = taxon.load_taxon_constraints_sqlite(args.db)
    if not constraints:
        logger.warning("GO taxon constraints not applied: the SQLite database has no go_taxon_constraints table "
                       "(build_databases.py --tables taxonomy --update) and no --taxon-constraints-obo was given")
        return {}
    lineage = []
    if getattr(args, "taxdump", None):
        lineage = taxon.lineage_from_taxdump(args.taxdump, taxid)
    if not lineage and args.db:
        lineage = taxon.lineage_sqlite(args.db, taxid)
    if not lineage:
        logger.warning(f"GO taxon constraints not applied: no lineage for taxon {taxid} (the SQLite database "
                       f"has no taxonomy_nodes table and no --taxdump directory was given)")
        return {}
    name = taxon.taxon_name_sqlite(args.db, lineage[0]) if args.db else ""
    forbidden = taxon.forbidden_go_terms(constraints, unions, lineage, go_hierarchy)
    logger.info(f"GO taxon constraints for taxon {taxid}{' (' + name + ')' if name else ''}: lineage of "
                f"{len(lineage)} taxa; {len(forbidden):,} GO terms excluded")
    return forbidden


def _require_file(path, label):
    if path and not Path(path).exists():
        raise GOAnnotateError(f"{label} not found: {path}")


def run(args) -> None:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or output_dir.resolve().name
    # Annotation files at the top level; reports and intermediate files in
    # their own directories.
    reports_dir = output_dir / "reports"
    work_dir = output_dir / "intermediate"
    reports_dir.mkdir(exist_ok=True)
    work_dir.mkdir(exist_ok=True)
    _migrate_old_layout(output_dir, work_dir, prefix)

    for path, label in ((args.gff3, "GFF3 file"), (args.genome, "Genome FASTA"),
                        (getattr(args, "proteins", None), "Protein FASTA"),
                        (args.transcripts, "Transcripts FASTA"), (args.blast_results, "BLAST results"),
                        (args.db, "SQLite database"), (args.go_mapping, "GO mapping file"),
                        (args.ncbi_idmapping, "NCBI ID mapping file"), (args.ncbi_geneinfo, "NCBI gene_info file"),
                        (args.bad_names, "Bad names file"), (args.tx2gene, "Transcript-to-gene file")):
        _require_file(path, label)

    args.go_obo = resolve_go_obo(args.go_obo)
    reference.resolve_references(args)

    ############################################################
    # Stage 1: gene universe and query sequences
    ############################################################
    genes = None
    gff_stats = None
    extraction = None

    protein_check = None
    if args.gff3 and (args.genome or getattr(args, "proteins", None)):
        genes, gff_stats = gff3.parse_gff3(args.gff3)
        gff3.assign_keys(genes, args.gene_id_attr, args.transcript_id_attr)
        universe = gff3.universe_from_genes(genes)
        if args.genome:
            genome = sequences.load_genome(args.genome)
            extraction = sequences.extract_sequences(
                genes, genome, work_dir, prefix,
                default_table=args.codon_table, mito_seqids=args.mito_seqid, report_dir=reports_dir)
            del genome
        if getattr(args, "proteins", None):
            if not sequences.fasta_is_protein(args.proteins):
                raise GOAnnotateError(f"{args.proteins} does not look like a protein FASTA; "
                                      f"use --transcripts for nucleotide sequences")
            protein_check = sequences.check_protein_fasta(
                args.proteins, genes, reports_dir, prefix,
                extracted_fasta=extraction.protein_fasta if extraction else None)
            query_fasta = Path(args.proteins)
            logger.info(f"Using the supplied proteins for the search: {args.proteins}")
        else:
            query_fasta = extraction.protein_fasta
        search_mode = "blastp"
        annotated_fasta_source = query_fasta
    elif args.transcripts:
        universe = sequences.universe_from_fasta(args.transcripts)
        if args.tx2gene:
            _apply_tx2gene(universe, args.tx2gene)
        query_fasta = Path(args.transcripts)
        search_mode = "blastp" if sequences.fasta_is_protein(args.transcripts) else "blastx"
        if search_mode == "blastp":
            logger.info("Input FASTA contains protein sequences; Diamond blastp will be used")
        annotated_fasta_source = query_fasta
        if args.gff3:
            genes, gff_stats = gff3.parse_gff3(args.gff3)
            gff3.assign_keys(genes, args.gene_id_attr, args.transcript_id_attr)
    else:
        raise GOAnnotateError("Provide --gff3 with --genome, or --transcripts")

    ############################################################
    # Stage 2: Diamond search (or precomputed results)
    ############################################################
    coverage = _coverage_filter(args)
    if args.blast_results:
        results_path = Path(args.blast_results)
    else:
        results_path = search.run_searches(
            mode=search_mode, query=query_fasta, output_dir=work_dir, prefix=prefix,
            swissprot_db=args.diamond_swissprot, trembl_db=args.diamond_trembl,
            threads=args.threads or 0, evalue=args.evalue, max_target_seqs=args.max_target_seqs,
            sensitivity=args.sensitivity, extra_args=args.diamond_args,
            diamond_bin=args.diamond_bin, force=args.force_search, need_lengths=coverage.active)

    ############################################################
    # Stage 3: annotation
    ############################################################
    transcript_hits, blast_stats = blast.parse_blast_results(str(results_path), evalue_threshold=args.evalue)

    unknown = [q for q in transcript_hits if q not in universe.transcript_to_gene]
    if unknown:
        logger.warning(f"{len(unknown):,} query IDs in the BLAST results are not in the gene universe "
                       f"(first: {unknown[0]}); their hits are kept under their own IDs")

    extra = getattr(args, "bad_names_extra", None) or []
    for path in extra:
        _require_file(path, "Bad names file")
    bad_name_filter = BadNameFilter(args.bad_names, extra_files=extra)
    output.configure_products(bad_name_filter, cleanup=not getattr(args, "no_product_cleanup", False))
    output.configure_existing(bad_name_filter, keep_identifiers=getattr(args, "keep_existing_identifiers", False))
    name_audit = NameAudit()
    if coverage.active and blast_stats.hits_with_lengths == 0 and blast_stats.hits_after_evalue_filter:
        logger.warning("The search results carry no qlen/slen columns; the hit coverage filter tests identity "
                       "only. Rerun with --force-search (or pass --hit-coverage off) to apply it fully")
    filtered_hits, coverage_status = blast.filter_and_select_hits(
        transcript_hits, bad_name_filter, top_n=args.top_n, audit=name_audit, coverage=coverage)

    accessions = mappings.collect_accessions_from_blast(filtered_hits)
    symbol_pairs = None
    if not args.no_resolve_orthologs:
        symbol_pairs = blast.collect_symbol_pairs(filtered_hits, str(args.reference_taxid))
    go_mapping, resolver, nomenclature = mappings.load_reference_data(
        args.db, args.go_mapping, args.ncbi_idmapping, args.ncbi_geneinfo, accessions,
        resolve_orthologs=not args.no_resolve_orthologs,
        reference_taxid=str(args.reference_taxid), anchor_taxid=str(args.ortholog_anchor_taxid),
        symbol_pairs=symbol_pairs, go_evidence_policy=getattr(args, "go_evidence", "all"))
    resolved_hits = blast.resolve_symbols(filtered_hits, resolver, bad_name_filter, audit=name_audit)
    go_hierarchy = GOHierarchy(args.go_obo)
    forbidden_go = _taxon_constraints(args, go_hierarchy)

    word_weights = None
    scheme = getattr(args, "word_weighting", "idf")
    if scheme != "off":
        word_weights = annotate_mod.WordWeights(
            scheme=scheme, stop_words=annotate_mod.load_stop_words(getattr(args, "stop_words", None)))

    coding_gene_ids = None
    if universe.coding:
        coding_gene_ids = {g for g, c in universe.coding.items() if c}

    ############################################################
    # Stage 2b: modules that run before the vote (orthology)
    ############################################################
    # The orthology module runs OrthoFinder now, so that the relatives'
    # ortholog symbols can steer the vote; it is integrated after Stage 3b.
    selected = modules.select_modules(args)
    modules_auto = (args.modules or "").strip().lower() == "auto"
    early_modules = []
    ortholog_prior = {}
    ctx = None
    if selected:
        ctx = modules.ModuleContext(
            query_fasta=Path(query_fasta), search_mode=search_mode, work_dir=work_dir,
            reports_dir=reports_dir, prefix=prefix, threads=args.threads or 0,
            force=getattr(args, "force_modules", False), universe=universe, genes=genes,
            annotations={}, go_hierarchy=go_hierarchy, forbidden_go=set(forbidden_go),
            namespaces=list(args.namespace), bad_name_filter=bad_name_filter,
            symbol_case_fn=lambda s: annotate_mod.format_symbol(s, args.symbol_case))
        early_modules = modules.start_modules([m for m in selected if m.runs_before_vote], ctx, auto=modules_auto)
        for m in early_modules:
            ortholog_prior.update(m.prior())

    annotations = annotate_mod.annotate_queries(
        resolved_transcript_hits=resolved_hits,
        transcript_to_gene=universe.transcript_to_gene,
        go_mapping=go_mapping,
        go_hierarchy=go_hierarchy,
        bad_name_filter=bad_name_filter,
        all_gene_ids=universe.gene_ids,
        consensus_threshold=args.consensus_threshold,
        namespaces=args.namespace,
        resolver=resolver,
        threads=args.threads,
        symbol_case=args.symbol_case,
        coding_gene_ids=coding_gene_ids,
        symbol_score=args.symbol_score,
        max_gene_hits=getattr(args, "max_gene_hits", 300),
        word_weights=word_weights,
        forbidden_go=set(forbidden_go),
        go_iea_min_hits=getattr(args, "go_iea_min_hits", 1),
        coverage_status=coverage_status if coverage.active else None,
        ortholog_prior=ortholog_prior,
        reference_taxid=str(args.reference_taxid),
    )

    ############################################################
    # Stage 3b: genome-wide symbol assignment (reciprocal best hits)
    ############################################################
    assignment = None
    ref_symbols, ref_descriptions = {}, {}
    if args.no_reciprocal:
        logger.info("Reciprocal best-hit assignment disabled (--no-reciprocal)")
    elif search_mode != "blastp":
        logger.info("Reciprocal best-hit assignment needs protein queries; skipped for nucleotide input")
    else:
        if args.reference_proteome:
            ref_proteome = Path(args.reference_proteome)
            _require_file(ref_proteome, "Reference proteome")
        else:
            ref_proteome = assign.build_reference_proteome(
                args.diamond_swissprot, str(args.reference_taxid),
                work_dir / f"{prefix}_reference_proteome.faa", args.diamond_bin,
                trembl_db=args.diamond_trembl)
        ref_symbols, ref_descriptions = assign.read_reference_symbols(ref_proteome)
        # Reference symbols and names pass through the same bad-name filter
        # as database hits: reference proteins with uninformative symbols
        # (C1orf122, KIAA1191) are not used to name genes, and uninformative
        # reference names are not used as products.
        n_before = len(ref_symbols)
        ref_symbols = {acc: sym for acc, sym in ref_symbols.items()
                       if not bad_name_filter.is_bad_symbol(sym)}
        clean_descriptions = {}
        for sym, desc in ref_descriptions.items():
            cleaned = bad_name_filter.clean_name(desc)
            if cleaned and not bad_name_filter.is_bad_name(cleaned):
                clean_descriptions[sym] = cleaned
        logger.info(f"  Bad-name filter on the reference proteome: {n_before - len(ref_symbols):,} proteins "
                    f"with uninformative symbols excluded, {len(ref_descriptions) - len(clean_descriptions):,} "
                    f"uninformative names not used as products")
        ref_descriptions = clean_descriptions
        forward, reverse = assign.run_reciprocal_searches(
            Path(query_fasta), ref_proteome, work_dir, prefix, threads=args.threads or 0,
            evalue=args.evalue, sensitivity_flags=search.SENSITIVITY_FLAGS.get(args.sensitivity, []),
            diamond_bin=args.diamond_bin, force=args.force_search)
        assignment = assign.assign_symbols(
            annotations, universe, forward, reverse, ref_symbols,
            symbol_case_fn=lambda s: annotate_mod.format_symbol(s, args.symbol_case),
            ref_descriptions=ref_descriptions, support_fraction=args.symbol_support,
            min_coverage=args.symbol_min_coverage, symbol_policy=getattr(args, "symbol_policy", "ortholog"),
            nomenclature=nomenclature, drop_unmatched=getattr(args, "drop_unmatched_symbols", False),
            min_query_coverage=getattr(args, "symbol_min_query_coverage", 0.8),
            rbh_margin=getattr(args, "rbh_margin", 0.9), existence=assign.read_reference_existence(ref_proteome),
            name_threshold=args.consensus_threshold)
        if not getattr(args, "keep_paralog_designators", False):
            assignment.update(assign.strip_paralog_designators(
                annotations, ref_descriptions, bad_name_filter.normalize_for_comparison))
        for annot in annotations.values():
            annot.confidence, annot.confidence_score = annotate_mod.compute_confidence(annot)

        # Second nomenclature (human for an insect annotated with Drosophila
        # symbols): unambiguous reciprocal best hits in that species' proteome
        sec = args.secondary_reference_taxid
        if sec:
            if getattr(args, "secondary_reference_proteome", None):
                sec_proteome = Path(args.secondary_reference_proteome)
                _require_file(sec_proteome, "Secondary reference proteome")
            else:
                sec_proteome = assign.build_reference_proteome(
                    args.diamond_swissprot, str(sec), work_dir / f"{prefix}_secondary_{sec}_proteome.faa",
                    args.diamond_bin, trembl_db=args.diamond_trembl)
            sec_symbols, sec_descriptions = assign.read_reference_symbols(sec_proteome)
            sec_descriptions = {sym: bad_name_filter.clean_name(desc) for sym, desc in sec_descriptions.items()}
            sec_descriptions = {sym: d for sym, d in sec_descriptions.items() if d and not bad_name_filter.is_bad_name(d)}
            sec_symbols = {acc: s for acc, s in sec_symbols.items() if not bad_name_filter.is_bad_symbol(s)}
            sec_forward, sec_reverse = assign.run_reciprocal_searches(
                Path(query_fasta), sec_proteome, work_dir, f"{prefix}_secondary_{sec}",
                threads=args.threads or 0, evalue=args.evalue,
                sensitivity_flags=search.SENSITIVITY_FLAGS.get(args.sensitivity, []),
                diamond_bin=args.diamond_bin, force=args.force_search)
            sec_nomenclature = mappings.load_reference_nomenclature(args.db, str(sec)) if args.db else None
            assignment.update(assign.secondary_symbols(
                annotations, universe, sec_forward, sec_reverse, sec_symbols,
                symbol_case_fn=lambda s: annotate_mod.format_symbol(s, args.symbol_case),
                nomenclature=sec_nomenclature, min_coverage=args.symbol_min_coverage,
                min_query_coverage=getattr(args, "symbol_min_query_coverage", 0.8),
                rbh_margin=getattr(args, "rbh_margin", 0.9)))
            assignment["secondary_reference"] = f"{reference.label(sec)} ({sec})"
            assignment.update(assign.secondary_names(annotations, sec_descriptions, ref_descriptions,
                                                     name_threshold=args.consensus_threshold))

    assign.add_symbol_aliases(annotations, nomenclature, bad_name_filter.is_bad_symbol,
                              symbol_case_fn=lambda s: annotate_mod.format_symbol(s, args.symbol_case))

    ############################################################
    # Stage 3c: ortholog symbols, reference-gene GO terms, other modules
    ############################################################
    modules_run = []
    if ctx is not None:
        ctx.annotations = annotations
        ctx.reference_descriptions = ref_descriptions
        modules_run = modules.integrate_modules(early_modules, ctx)
    reference_go = refgo.add_reference_go(
        annotations, ref_symbols, args.db, go_hierarchy, set(forbidden_go), list(args.namespace),
        evidence_policy=getattr(args, "go_evidence", "all"), mode=getattr(args, "reference_go", "add"))
    if ctx is not None:
        modules_run += modules.run_modules([m for m in selected if not m.runs_before_vote], ctx, auto=modules_auto)

    # Concordance of the final symbols and products with each gene's hits, and
    # the confidence recomputed from it (with the orthology evidence)
    annotate_mod.final_concordance(annotations, consensus_threshold=args.consensus_threshold)

    ############################################################
    # Stage 4: comparison with existing annotation
    ############################################################
    comparison = None
    if genes is not None and not args.no_compare:
        has_existing = any(compare.existing_symbol(g, args.existing_symbol_attr)
                           or compare.existing_product(g, args.existing_product_attr) for g in genes)
        if has_existing:
            comparison = compare.compare_with_existing(
                annotations, genes, reports_dir, prefix,
                symbol_attrs=args.existing_symbol_attr, product_attrs=args.existing_product_attr,
                gene_strata=args.strata_gene, transcript_strata=args.strata_transcript,
                is_bad_symbol=bad_name_filter.is_bad_symbol)
        else:
            logger.info("No existing functional annotation found in the GFF3; comparison skipped")

    ############################################################
    # Stage 5: outputs
    ############################################################
    write_name_audit(name_audit, bad_name_filter, reports_dir / f"{prefix}_name_audit.txt",
                     annotations=annotations)
    output.write_go_evidence(annotations, go_hierarchy, reports_dir / f"{prefix}_go_evidence.tsv",
                             forbidden_go=forbidden_go)
    product_audit = output.product_cleanup_audit(annotations)
    if product_audit:
        logger.info("Product naming rules applied: " + "; ".join(
            f"{n:,}x {rule[:60]}" for rule, n in sorted(product_audit.items(), key=lambda kv: -kv[1])))
    protein_lengths = None
    if search_mode == "blastp":
        protein_lengths = {}
        for header, seq in sequences.read_fasta(str(query_fasta)):
            protein_lengths[header.split()[0] if header else ""] = len(seq.rstrip("*"))
    output.write_annotation_tables(output_dir, prefix, annotations, universe, genes=genes,
                                   fasta_headers=str(annotated_fasta_source), mode=args.annotation_mode,
                                   protein_lengths=protein_lengths, modules=modules_run)
    output.write_mapping_files(annotations, output_dir, prefix, universe=universe, evidence_dir=reports_dir)

    # Annotated sequence files, named by content: amino acid (AA), coding
    # sequence (CDS) and spliced transcript sequences.
    fasta_outputs = []
    if search_mode == "blastp":
        fasta_outputs.append((annotated_fasta_source, f"{prefix}_annotated_AA.faa"))
    else:
        fasta_outputs.append((annotated_fasta_source, f"{prefix}_annotated_transcripts.fna"))
    if extraction is not None:
        fasta_outputs.append((extraction.cds_fasta, f"{prefix}_annotated_CDS.fna"))
        fasta_outputs.append((extraction.transcript_fasta, f"{prefix}_annotated_transcripts.fna"))
    for source, name in fasta_outputs:
        output.write_annotated_fasta(str(source), annotations, output_dir / name,
                                     transcript_to_gene=universe.transcript_to_gene,
                                     mode=args.annotation_mode)
    if genes is not None:
        output.write_annotated_gff3(args.gff3, output_dir / f"{prefix}_annotated.gff3",
                                    annotations, genes, mode=args.annotation_mode,
                                    hypothetical=args.hypothetical, modules=modules_run)
    output.write_summary(output_dir, prefix, annotations, args, blast_stats=blast_stats,
                         universe=universe, extraction=extraction, gff_stats=gff_stats,
                         comparison=comparison, search_results=str(results_path),
                         assignment=assignment, protein_check=protein_check,
                         coverage=coverage, forbidden_go=forbidden_go, product_audit=product_audit,
                         word_weights=word_weights, modules=selected, reference_go=reference_go)


# Files of the earlier output layout (everything at the top level)
_OLD_INTERMEDIATES = ("diamond_results.tsv", "diamond_SwissProt.tsv", "diamond_TrEMBL.tsv",
                      "reference_proteome.faa", "reference_proteome.dmnd", "proteins.dmnd",
                      "vs_reference.tsv", "reference_vs_query.tsv")
_OLD_OUTPUTS = ("proteins.faa", "cds.fna", "transcripts.fna", "annotated.fasta", "extraction_report.tsv",
                "protein_input_check.tsv", "comparison.tsv", "comparison_summary.txt", "review.tsv",
                "annotation_evidence.tsv", "name_audit.txt")


def _migrate_old_layout(output_dir: Path, work_dir: Path, prefix: str) -> None:
    """
    Move files of the earlier single-directory layout out of the top level:
    search results and reference files into intermediate/ (so they are
    reused), outdated sequence and report files into intermediate/old_layout/.
    Nothing is deleted. Symbolic links are recreated with absolute targets.
    """
    moved = []
    for suffix in _OLD_INTERMEDIATES + _OLD_OUTPUTS:
        src = output_dir / f"{prefix}_{suffix}"
        if not (src.exists() or src.is_symlink()):
            continue
        dest_dir = work_dir if suffix in _OLD_INTERMEDIATES else work_dir / "old_layout"
        dest_dir.mkdir(exist_ok=True)
        dest = dest_dir / src.name
        if dest.exists() or dest.is_symlink():
            continue
        if src.is_symlink():
            target = (src.parent / os.readlink(src)).resolve()
            dest.symlink_to(target)
            src.unlink()
        else:
            shutil.move(str(src), str(dest))
        moved.append(src.name)
    if moved:
        logger.info(f"Moved {len(moved)} files of the earlier output layout into {work_dir}: {', '.join(moved)}")


def _apply_tx2gene(universe, path: str) -> None:
    """Override the transcript-to-gene mapping with a two-column file (transcript_id, gene_id)."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0].startswith("#"):
                continue
            mapping[parts[0]] = parts[1]
    if not mapping:
        raise GOAnnotateError(f"No transcript-to-gene pairs read from {path}")

    gene_ids = []
    seen = set()
    gene_to_transcripts = {}
    n_mapped = 0
    for tx in list(universe.transcript_to_gene):
        gene_id = mapping.get(tx)
        if gene_id is None:
            gene_id = tx
        else:
            n_mapped += 1
        universe.transcript_to_gene[tx] = gene_id
        gene_to_transcripts.setdefault(gene_id, []).append(tx)
        if gene_id not in seen:
            seen.add(gene_id)
            gene_ids.append(gene_id)
    universe.gene_ids = gene_ids
    universe.gene_to_transcripts = gene_to_transcripts
    logger.info(f"Applied transcript-to-gene mapping from {path}: {n_mapped:,} transcripts mapped, "
                f"{len(gene_ids):,} genes")
