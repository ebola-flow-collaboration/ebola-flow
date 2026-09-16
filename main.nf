#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHAC-NML/ebola-flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/PHAC-NML/ebola-flow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PATHOPLEXUS_FETCH    } from './modules/nf-core/pathoplexus_fetch/main'
include { PATHOPLEXUS_FETCH as PATHOPLEXUS_FETCH_OUTGROUP } from './modules/nf-core/pathoplexus_fetch/main'
include { CHECK_NEW_SUBMISSIONS } from './modules/nf-core/check_new_submissions/main'
include { CLEAN_METADATA       } from './modules/nf-core/clean_metadata/main'
include { CLEAN_OUTGROUP       } from './modules/nf-core/clean_outgroup/main'
include { SEQKIT_GREP          } from './modules/nf-core/seqkit_grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_OUTGROUP } from './modules/nf-core/seqkit_grep/main'
include { TRIM_SEQUENCES       } from './modules/nf-core/trim_sequences/main'
include { TRIM_SEQUENCES as TRIM_SEQUENCES_OUTGROUP } from './modules/nf-core/trim_sequences/main'
include { NEXTCLADE_DATASETGET } from './modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN        } from './modules/nf-core/nextclade/run/main'
include { IQTREE               } from './modules/nf-core/iqtree/main'
include { IQTREE as IQTREE_POST_FILTER } from './modules/nf-core/iqtree/main'
include { ITERATIVE_TREETIME  } from './modules/nf-core/iterative_treetime/main'
include { ITERATIVE_TREETIME as ITERATIVE_TREETIME_POST } from './modules/nf-core/iterative_treetime/main'
include { EMERGENCE_ANALYSIS   } from './modules/nf-core/emergence_analysis/main'
include { COMBINE_REMOVED_SEQUENCES } from './modules/nf-core/combine_removed_sequences/main'
include { DOWNSAMPLE           } from './modules/nf-core/downsample/main'
include { BEAST2_XML_GEN       } from './modules/nf-core/beast2_xml_gen/main'
include { BEAST2               } from './modules/nf-core/beast2/main'
include { BEAST_PYPE_STATIC_DIAGNOSE } from './modules/nf-core/beast_pype_static_diagnose/main'
include { BEAST_PYPE_PARAMETERS_REPORT } from './modules/nf-core/beast_pype_parameters_report/main'
include { TREEANNOTATOR              } from './modules/nf-core/treeannotator/main'
include { BEAST_PYPE_SUMMARY_TREE_REPORT } from './modules/nf-core/beast_pype_summary_tree_report/main'
include { DATA_AND_RUN_STATS_REPORT } from './modules/nf-core/data_and_run_stats_report/main'
include { TEMPEST_METADATA                   } from './modules/nf-core/tempest_metadata/main'
include { TEMPEST_METADATA as TEMPEST_METADATA_NO_REF } from './modules/nf-core/tempest_metadata/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_ebola-flow_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_ebola-flow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
workflow EBOLAFLOW {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // MODULE: Optionally check for new submissions today before proceeding
    //
    if (params.check_new_submissions && !params.fasta) {
        ch_check_meta = channel.of(
            [id: params.pathoplexus_organism, organism: params.pathoplexus_organism]
        )
        CHECK_NEW_SUBMISSIONS(ch_check_meta)

        // Gate: only proceed if has_submissions.txt was emitted (optional output)
        ch_proceed = CHECK_NEW_SUBMISSIONS.out.proceed
            .map { meta, txt -> meta }
    }

    //
    // MODULE: Fetch sequences and metadata from Pathoplexus (if --fasta not provided)
    //
    if (params.fasta) {
        ch_fasta    = channel.fromPath(params.fasta, checkIfExists: true)
            .map { fasta -> tuple([id: fasta.baseName], fasta) }
        ch_metadata = channel.fromPath(params.metadata, checkIfExists: true)
        // No Pathoplexus fetch when a FASTA is supplied directly.
        ch_raw_metadata          = channel.fromPath(params.metadata, checkIfExists: true)
        ch_cleaned_metadata      = channel.fromPath(params.metadata, checkIfExists: true)
        ch_pathoplexus_timestamp = channel.of(file('NO_FILE'))
    } else {
        ch_pathoplexus_meta = channel.of(
            [id: params.pathoplexus_organism, organism: params.pathoplexus_organism, date_from: params.pathoplexus_date_from ?: '']
        )

        // If check_new_submissions is enabled, gate PATHOPLEXUS_FETCH on the check result
        if (params.check_new_submissions) {
            ch_pathoplexus_meta = ch_proceed
                .map { meta ->
                    [id: params.pathoplexus_organism, organism: params.pathoplexus_organism, date_from: params.pathoplexus_date_from ?: '']
                }
        }

        PATHOPLEXUS_FETCH(ch_pathoplexus_meta)

        // Capture the Pathoplexus fetch date/time for the data and run stats report
        ch_pathoplexus_timestamp = PATHOPLEXUS_FETCH.out.timestamp.map { meta, ts -> ts }

        //
        // MODULE: Clean metadata (JSON→TSV, filter excluded IDs, filter invalid dates)
        //
        ch_exclude_ids = params.exclude_ids
            ? channel.fromPath(params.exclude_ids, checkIfExists: true)
            : channel.of(file('NO_FILE'))

        ch_clean_meta_input = PATHOPLEXUS_FETCH.out.metadata
            .combine(ch_exclude_ids)
            .map { meta, json_file, ids -> tuple(meta, json_file, ids) }

        CLEAN_METADATA(ch_clean_meta_input)
        ch_metadata = CLEAN_METADATA.out.metadata.map { meta, tsv -> tsv }
        // Raw Pathoplexus pull (unfiltered) and post-clean_metadata data
        ch_raw_metadata     = CLEAN_METADATA.out.metadata_raw.map { meta, tsv -> tsv }
        ch_cleaned_metadata = CLEAN_METADATA.out.metadata_full.map { meta, tsv -> tsv }

        //
        // MODULE: Remove sequences listed in ids_to_remove using seqkit
        //
        ch_seqkit_input = PATHOPLEXUS_FETCH.out.fasta
            .combine(CLEAN_METADATA.out.ids_to_remove.map { meta, ids -> ids })
            .map { meta, fasta, ids -> tuple(meta, fasta, ids) }

        SEQKIT_GREP(ch_seqkit_input)
        ch_fasta = SEQKIT_GREP.out.fasta

        //
        // MODULE: Optionally trim sequences to a maximum length
        //
        if (params.trim_max_length) {
            TRIM_SEQUENCES(ch_fasta)
            ch_fasta = TRIM_SEQUENCES.out.fasta
        }
    }

    //
    // MODULE: Download nextstrain/ebola/bdbv dataset
    // When check_new_submissions is enabled, gate on ch_proceed so this waits
    //
    if (params.check_new_submissions && !params.fasta) {
        ch_nextclade_trigger = ch_proceed.map { it -> params.nextclade_dataset }
    } else {
        ch_nextclade_trigger = channel.of(params.nextclade_dataset)
    }

    NEXTCLADE_DATASETGET(
        ch_nextclade_trigger,
        params.nextclade_dataset_tag ?: []
    )

    //
    // Extract reference ID from Nextclade dataset
    //
    ch_ref_id = NEXTCLADE_DATASETGET.out.dataset
        .map { dataset_dir ->
            def ref_file = file("${dataset_dir}/reference.fasta")
            def first_line = ref_file.text?.readLines()?.getAt(0)
            first_line ? first_line.replaceFirst('^>', '').split('\\s+')[0] : 'REFERENCE'
        }

    //
    // MODULE: Optionally fetch outgroup sequences, clean them, and merge with main FASTA before alignment
    //
    def has_outgroup = (params.outgroup_date_from && params.outgroup_date_from != 'null') ||
                       (params.outgroup_date_to && params.outgroup_date_to != 'null')
    def build_tree = params.build_initial_tree && params.build_initial_tree.toString() != 'false'
    def build_with_ref = params.build_initial_tree_with_reference && params.build_initial_tree_with_reference.toString() != 'false'

    if ((build_tree || build_with_ref) && has_outgroup) {
        ch_outgroup_meta = channel.of(
            [id: "${params.pathoplexus_organism}_outgroup", organism: params.pathoplexus_organism]
        )
        PATHOPLEXUS_FETCH_OUTGROUP(ch_outgroup_meta)

        // Clean outgroup: identify reference duplicates and round dates
        CLEAN_OUTGROUP(
            PATHOPLEXUS_FETCH_OUTGROUP.out.metadata,
            ch_ref_id
        )

        // Remove reference duplicate from outgroup FASTA
        ch_seqkit_outgroup_input = PATHOPLEXUS_FETCH_OUTGROUP.out.fasta
            .combine(CLEAN_OUTGROUP.out.ids_to_remove.map { meta, ids -> ids })
            .map { meta, fasta, ids -> tuple(meta, fasta, ids) }

        SEQKIT_GREP_OUTGROUP(ch_seqkit_outgroup_input)

        // Optionally trim outgroup sequences to max length
        if (params.trim_max_length) {
            TRIM_SEQUENCES_OUTGROUP(SEQKIT_GREP_OUTGROUP.out.fasta)
            ch_outgroup_fasta_clean = TRIM_SEQUENCES_OUTGROUP.out.fasta
        } else {
            ch_outgroup_fasta_clean = SEQKIT_GREP_OUTGROUP.out.fasta
        }

        // Merge main FASTA + cleaned outgroup FASTA before alignment
        ch_fasta_for_alignment = ch_fasta
            .combine(ch_outgroup_fasta_clean.map { meta, fasta -> fasta })
            .map { meta, main_fasta, outgroup_fasta ->
                tuple(meta, [main_fasta, outgroup_fasta])
            }
    } else {
        ch_fasta_for_alignment = ch_fasta
    }

    //
    // MODULE: Align sequences using Nextclade
    //
    NEXTCLADE_RUN(
        ch_fasta_for_alignment,
        NEXTCLADE_DATASETGET.out.dataset
    )

    //
    // MODULE: Optionally build a genetic distance tree with IQ-TREE
    //
    if (build_tree || build_with_ref) {
        if (build_with_ref && has_outgroup) {
            // Root at MRCA of reference + outgroup sequences (comma-separated for IQ-TREE -o)
            ch_outgroup_ids_for_iqtree = SEQKIT_GREP_OUTGROUP.out.fasta
                .map { meta, fasta ->
                    fasta.readLines().findAll { it.startsWith('>') }.collect { it.replaceFirst('^>', '').split('\\s+')[0] }
                }
            ch_outgroup = ch_ref_id
                .combine(ch_outgroup_ids_for_iqtree)
                .map { items -> items.flatten().join(',') }
        } else if (build_with_ref) {
            // Reference only (no outgroup)
            ch_outgroup = ch_ref_id
        } else if (has_outgroup) {
            // Root at first outgroup sequence when using outgroup without reference
            ch_outgroup = SEQKIT_GREP_OUTGROUP.out.fasta
                .map { meta, fasta ->
                    fasta.readLines().find { it.startsWith('>') }?.replaceFirst('^>', '')?.split('\\s+')[0] ?: ''
                }
        } else {
            ch_outgroup = channel.of('')
        }

        ch_iqtree_input = NEXTCLADE_RUN.out.fasta_aligned
            .map { meta, fasta -> tuple(meta, fasta, []) }

        IQTREE(
            ch_iqtree_input,
            [], // tree_te
            [], // lmclust
            [], // mdef
            [], // partitions_equal
            [], // partitions_proportional
            [], // partitions_unlinked
            [], // guide_tree
            [], // sitefreq_in
            [], // constraint_tree
            [], // trees_z
            [], // suptree
            [], // trees_rf
            ch_outgroup
        )

        if (build_with_ref) {
            // Generate TempEst metadata with reference included (for iqtree output)
            ch_tempest_meta_input = CLEAN_METADATA.out.metadata

            // Get outgroup dates channel (if outgroup was fetched)
            if (has_outgroup) {
                ch_outgroup_dates = CLEAN_OUTGROUP.out.dates.map { meta, tsv -> tsv }
            } else {
                ch_outgroup_dates = channel.of(file('NO_FILE'))
            }

            TEMPEST_METADATA(
                ch_tempest_meta_input,
                ch_ref_id,
                channel.of(params.reference_date ?: ''),
                channel.of(true),
                ch_outgroup_dates
            )

            // Collect outgroup + ref IDs for rerooting/pruning in TreeTime
            if (has_outgroup) {
                ch_outgroup_ids = SEQKIT_GREP_OUTGROUP.out.fasta
                    .map { meta, fasta ->
                        fasta.readLines().findAll { it.startsWith('>') }.collect { it.replaceFirst('^>', '').split('\\s+')[0] }
                    }
                    .flatten()
                    .collect()

                ch_reroot_ids = ch_ref_id
                    .combine(ch_outgroup_ids)
                    .map { items -> items.flatten() }
            } else {
                ch_reroot_ids = ch_ref_id.map { id -> [id] }
            }

            //
            // MODULE: ITERATIVE_TREETIME — iterative clock filtering, stats, plots, and tree outputs
            // Removes outliers until convergence, prunes outgroup/reference (--remove-root),
            // produces temporal tree, genetic distance tree, stats, and plots.
            //
            ch_treetime_input = IQTREE.out.phylogeny
                .combine(NEXTCLADE_RUN.out.fasta_aligned.map { meta, fasta -> fasta })
                .combine(ch_metadata)
                .map { meta, tree, fasta, metadata ->
                    tuple(meta, tree, fasta, metadata)
                }

            ITERATIVE_TREETIME(
                ch_treetime_input,
                ch_reroot_ids,
                channel.of(params.treetime_clock_filter),
                channel.of(params.treetime_clock_filter_method ?: 'local'),
                channel.of(params.treetime_clock_model ?: 'strict'),
                channel.of(params.treetime_remove_future_tips),
                channel.of(params.beast_seed)
            )

            //
            // MODULE: EMERGENCE_ANALYSIS — variant emergence surveillance report
            // Combines the raw diversity tree (IQ-TREE), the clock-cleaned genetic
            // distance tree, the temporal tree, and metadata into an R Markdown report.
            //
            ch_emergence_input = IQTREE.out.phylogeny
                .join(ITERATIVE_TREETIME.out.genetic_distance_tree)
                .join(ITERATIVE_TREETIME.out.temporal_tree)
                .join(ITERATIVE_TREETIME.out.metadata)
                .map { meta, div_raw, div_clean, time_tree, metadata ->
                    tuple(meta, div_raw, div_clean, time_tree, metadata)
                }

            EMERGENCE_ANALYSIS(ch_emergence_input)

            //
            // MODULE: COMBINE_REMOVED_SEQUENCES — consolidate all removed sequences
            // Combines clean_metadata's excluded_samples.csv (exclude-list hits and
            // invalid-date removals) with iterative_treetime's all_outliers.tsv
            // (clock outliers) into one CSV documenting what was removed and why.
            //
            ch_removed_input = CLEAN_METADATA.out.excluded_samples
                .join(ITERATIVE_TREETIME.out.outliers)
                .map { meta, excluded, outliers ->
                    tuple(meta, excluded, outliers)
                }

            COMBINE_REMOVED_SEQUENCES(ch_removed_input)

            //
            // MODULE: DOWNSAMPLE (optional) — residual-based temporal downsampling
            // IQTREE_POST_FILTER and ITERATIVE_TREETIME_POST only run if downsampling took place
            //
            if (params.downsample_to) {
                ch_downsample_input = ITERATIVE_TREETIME.out.genetic_distance_tree
                    .join(ITERATIVE_TREETIME.out.fasta)
                    .join(ITERATIVE_TREETIME.out.metadata)
                    .map { meta, tree, fasta, metadata ->
                        tuple(meta, tree, fasta, metadata)
                    }

                DOWNSAMPLE(
                    ch_downsample_input,
                    ITERATIVE_TREETIME.out.clock_stats.map { meta, yml -> yml },
                    channel.of(params.downsample_to),
                    channel.of(params.beast_seed)
                )

                // Rebuild tree on downsampled data
                ch_iqtree_post_input = DOWNSAMPLE.out.fasta
                    .map { meta, fasta -> tuple(meta, fasta, []) }

                IQTREE_POST_FILTER(
                    ch_iqtree_post_input,
                    [], [], [], [], [], [], [], [], [], [], [], [],
                    channel.of('')
                )

                // Re-run ITERATIVE_TREETIME on downsampled data (no clock filtering)
                ch_treetime_post_input = IQTREE_POST_FILTER.out.phylogeny
                    .combine(DOWNSAMPLE.out.fasta.map { meta, fasta -> fasta })
                    .combine(DOWNSAMPLE.out.metadata.map { meta, csv -> csv })
                    .map { meta, tree, fasta, metadata ->
                        tuple(meta, tree, fasta, metadata)
                    }

                ITERATIVE_TREETIME_POST(
                    ch_treetime_post_input,
                    channel.of([]),         // no reroot
                    channel.of('0'),        // clock-filter=0 (no outlier removal)
                    channel.of('local'),
                    channel.of(params.treetime_clock_model ?: 'strict'),
                    channel.of(params.treetime_remove_future_tips),
                    channel.of(params.beast_seed)
                )

                // Outputs for BEAST2_XML_GEN
                ch_aligned_fasta = DOWNSAMPLE.out.fasta
                ch_metadata = DOWNSAMPLE.out.metadata.map { meta, csv -> csv }
                ch_initial_tree = ITERATIVE_TREETIME_POST.out.temporal_tree.map { meta, tree -> tree }
            } else {
                // No downsampling — use ITERATIVE_TREETIME outputs directly
                ch_aligned_fasta = ITERATIVE_TREETIME.out.fasta
                ch_metadata = ITERATIVE_TREETIME.out.metadata.map { meta, tsv -> tsv }
                ch_initial_tree = ITERATIVE_TREETIME.out.temporal_tree.map { meta, tree -> tree }
            }

            // Generate TempEst metadata without reference/outgroup (for treetime output)
            TEMPEST_METADATA_NO_REF(
                ch_tempest_meta_input,
                ch_ref_id,
                channel.of(params.reference_date ?: ''),
                channel.of(false),
                channel.of(file('NO_FILE'))
            )
        } else {
            // build_initial_tree (without reference) — may still have outgroup to prune
            if (has_outgroup) {
                ch_outgroup_ids_only = SEQKIT_GREP_OUTGROUP.out.fasta
                    .map { meta, fasta ->
                        fasta.readLines().findAll { it.startsWith('>') }.collect { it.replaceFirst('^>', '').split('\\s+')[0] }
                    }
                    .flatten()
                    .collect()
                    .map { ids -> ids }

                ch_reroot_ids_no_ref = ch_outgroup_ids_only
            } else {
                ch_reroot_ids_no_ref = channel.of([])
            }

            def do_remove = has_outgroup

            // ITERATIVE_TREETIME (no reference — use outgroup or default rerooting)
            ch_treetime_input_noref = IQTREE.out.phylogeny
                .combine(NEXTCLADE_RUN.out.fasta_aligned.map { meta, fasta -> fasta })
                .combine(ch_metadata)
                .map { meta, tree, fasta, metadata -> tuple(meta, tree, fasta, metadata) }

            ITERATIVE_TREETIME(
                ch_treetime_input_noref,
                ch_reroot_ids_no_ref,
                channel.of(params.treetime_clock_filter),
                channel.of(params.treetime_clock_filter_method ?: 'local'),
                channel.of(params.treetime_clock_model ?: 'strict'),
                channel.of(params.treetime_remove_future_tips),
                channel.of(params.beast_seed)
            )

            //
            // MODULE: EMERGENCE_ANALYSIS — variant emergence surveillance report
            //
            ch_emergence_input_noref = IQTREE.out.phylogeny
                .join(ITERATIVE_TREETIME.out.genetic_distance_tree)
                .join(ITERATIVE_TREETIME.out.temporal_tree)
                .join(ITERATIVE_TREETIME.out.metadata)
                .map { meta, div_raw, div_clean, time_tree, metadata ->
                    tuple(meta, div_raw, div_clean, time_tree, metadata)
                }

            EMERGENCE_ANALYSIS(ch_emergence_input_noref)

            //
            // MODULE: COMBINE_REMOVED_SEQUENCES — consolidate all removed sequences
            //
            ch_removed_input_noref = CLEAN_METADATA.out.excluded_samples
                .join(ITERATIVE_TREETIME.out.outliers)
                .map { meta, excluded, outliers ->
                    tuple(meta, excluded, outliers)
                }

            COMBINE_REMOVED_SEQUENCES(ch_removed_input_noref)

            // FILTER_SEQUENCES
            if (params.downsample_to) {
                ch_downsample_input_noref = ITERATIVE_TREETIME.out.genetic_distance_tree
                    .join(ITERATIVE_TREETIME.out.fasta)
                    .join(ITERATIVE_TREETIME.out.metadata)
                    .map { meta, tree, fasta, metadata ->
                        tuple(meta, tree, fasta, metadata)
                    }

                DOWNSAMPLE(
                    ch_downsample_input_noref,
                    ITERATIVE_TREETIME.out.clock_stats.map { meta, yml -> yml },
                    channel.of(params.downsample_to),
                    channel.of(params.beast_seed)
                )

                // Rebuild tree on downsampled data
                ch_iqtree_post_input_noref = DOWNSAMPLE.out.fasta
                    .map { meta, fasta -> tuple(meta, fasta, []) }

                IQTREE_POST_FILTER(
                    ch_iqtree_post_input_noref,
                    [], [], [], [], [], [], [], [], [], [], [], [],
                    channel.of('')
                )

                // Re-run ITERATIVE_TREETIME on downsampled data (no clock filtering)
                ch_treetime_post_input_noref = IQTREE_POST_FILTER.out.phylogeny
                    .combine(DOWNSAMPLE.out.fasta.map { meta, fasta -> fasta })
                    .combine(DOWNSAMPLE.out.metadata.map { meta, csv -> csv })
                    .map { meta, tree, fasta, metadata -> tuple(meta, tree, fasta, metadata) }

                ITERATIVE_TREETIME_POST(
                    ch_treetime_post_input_noref,
                    channel.of([]),         // no reroot
                    channel.of('0'),        // clock-filter=0
                    channel.of('local'),
                    channel.of(params.treetime_clock_model ?: 'strict'),
                    channel.of(params.treetime_remove_future_tips),
                    channel.of(params.beast_seed)
                )

                ch_aligned_fasta = DOWNSAMPLE.out.fasta
                ch_metadata = DOWNSAMPLE.out.metadata.map { meta, csv -> csv }
                ch_initial_tree = ITERATIVE_TREETIME_POST.out.temporal_tree.map { meta, tree -> tree }
            } else {
                // No downsampling — use ITERATIVE_TREETIME outputs directly
                ch_aligned_fasta = ITERATIVE_TREETIME.out.fasta
                ch_metadata = ITERATIVE_TREETIME.out.metadata.map { meta, tsv -> tsv }
                ch_initial_tree = ITERATIVE_TREETIME.out.temporal_tree.map { meta, tree -> tree }
            }
        }
    } else if (params.initial_tree) {
        ch_initial_tree = channel.fromPath(params.initial_tree, checkIfExists: true)
        ch_aligned_fasta = NEXTCLADE_RUN.out.fasta_aligned
    } else {
        ch_aligned_fasta = NEXTCLADE_RUN.out.fasta_aligned
    }

    //
    // MODULE: Generate BEAST2 XML from template
    //
    if (params.stop_after != 'treetime' && params.stop_after != 'iqtree' && params.stop_after != 'nextclade') {
    ch_template     = channel.fromPath(params.template_xml, checkIfExists: true)
    ch_xml_gen_input = ch_aligned_fasta
        .combine(ch_template)
        .combine(ch_metadata)
        .combine(ch_initial_tree)
        .map { meta, aligned_fasta, template, metadata, tree ->
            def new_meta = [id: params.run_name ?: template.baseName]
            tuple(new_meta, template, aligned_fasta, metadata, tree)
        }

    BEAST2_XML_GEN(ch_xml_gen_input)

    //
    // MODULE: Run BEAST2 with generated XML
    //
    if (params.stop_after != 'beast2_xml_gen') {
    def rng = new Random(params.beast_seed != null ? params.beast_seed as long : System.nanoTime())
    ch_seeds = channel.fromList(
        (1..params.beast_repeats).collect { i ->
            [i, rng.nextInt(999999) + 1]
        }
    )

    ch_beast_input = BEAST2_XML_GEN.out.xml
        .combine(ch_seeds)
        .map { meta, xml, rep_num, seed ->
            def new_meta = meta + [id: "rep_${rep_num}", seed: seed, template_id: meta.id, rep: rep_num]
            tuple(new_meta, xml)
        }

    BEAST2(ch_beast_input)

    //
    // MODULE: Generate the data and run stats report
    //   - Pathoplexus fetch date/time
    //   - metadata report on the raw Pathoplexus pull
    //   - metadata report after the clean_metadata step
    //   - metadata report on the outlier-removed (iterative_treetime) data
    //   - BEAST run times per model
    //
    // Build aligned model-name / .out-file lists from the BEAST2 screen logs.
    // Collect model names and out-files into two parallel, order-stable lists.
    // Wrapping each collected List in a single-element list ([it]) prevents
    // `combine` from unwrapping it into separate tuple elements.
    ch_beast_models = BEAST2.out.screen_log
        .map { meta, out -> meta.template_id }
        .collect()
        .map { models -> [models] }

    ch_beast_outs = BEAST2.out.screen_log
        .map { meta, out -> out }
        .collect()
        .map { outs -> [outs] }

    ch_report_input = ch_raw_metadata
        .combine(ch_cleaned_metadata)               // after clean_metadata
        .combine(ch_metadata)                       // outlier-removed metadata
        .combine(ch_pathoplexus_timestamp)
        .combine(ch_beast_models)
        .combine(ch_beast_outs)
        .map { raw_meta, cleaned_meta, filtered_meta, timestamp, models, outs ->
            def new_meta = [id: params.run_name ?: 'beast2']
            tuple(new_meta, raw_meta, cleaned_meta, filtered_meta, timestamp, models, outs)
        }

    DATA_AND_RUN_STATS_REPORT(ch_report_input)

    //
    // MODULE: Diagnose and merge BEAST2 outputs (per template)
    //
    ch_beast_pype_input = BEAST2.out.log
        .map { meta, log -> tuple(meta.template_id, log) }
        .mix(
            BEAST2.out.trees.map { meta, trees -> tuple(meta.template_id, trees) }
        )
        .groupTuple()
        .map { template_id, files ->
            def meta = [id: params.run_name ?: template_id]
            tuple(meta, files)
        }

    BEAST_PYPE_STATIC_DIAGNOSE(ch_beast_pype_input)

    //
    // MODULE: Summarise posterior trees with TreeAnnotator
    //
    ch_treeannotator_input = BEAST_PYPE_STATIC_DIAGNOSE.out.merged_trees
        .map { meta, trees ->
            def new_meta = [id: meta.id ?: 'beast2_mcc']
            tuple(new_meta, trees)
        }

    TREEANNOTATOR(ch_treeannotator_input)

    //
    // MODULE: Generate parameters report from merged log(s) and BEAST2 XML(s)
    //
    ch_all_merged_logs = BEAST_PYPE_STATIC_DIAGNOSE.out.merged_log
        .map { meta, merged_log -> merged_log }
        .collect()

    ch_all_xmls = BEAST2_XML_GEN.out.xml
        .map { meta, xml -> xml }
        .collect()

    ch_params_report_input = ch_all_merged_logs
        .combine(ch_all_xmls)
        .map { items ->
            def merged_logs = items.findAll { it.name.endsWith('.csv') }
            def xmls = items.findAll { it.name.endsWith('.xml') }
            def meta = [id: params.run_name ?: 'beast2_diagnostics']
            tuple(meta, merged_logs, xmls)
        }

    BEAST_PYPE_PARAMETERS_REPORT(ch_params_report_input)

    //
    // MODULE: Generate summary tree report from TreeAnnotator output
    // Join by template ID so each summary tree is paired with its corresponding XML
    //
    ch_summary_tree_report_input = TREEANNOTATOR.out.tree
        .map { meta, tree -> tuple(meta.id, meta, tree) }
        .join(
            BEAST2_XML_GEN.out.xml.map { meta, xml -> tuple(meta.id, xml) }
        )
        .map { id, meta, tree, xml ->
            tuple(meta, tree, xml)
        }

    BEAST_PYPE_SUMMARY_TREE_REPORT(ch_summary_tree_report_input)
    } // end stop_after != 'beast2_xml_gen'
    } // end stop_after != 'treetime'/'iqtree'/'nextclade'
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    EBOLAFLOW (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
