/*
 * ITERATIVE_TREETIME
 *
 * Runs iterative TreeTime clock filtering, generates stats, plots, and
 * temporal/genetic distance trees. Uses the TreeTime Python API directly.
 * Replaces the old ITERATIVE_TREETIME_FILTER + TREETIME_STATS_AND_PLOTS modules.
 */

process ITERATIVE_TREETIME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-treetime:latest' :
        'quay.io/biocontainers/mulled-v2-treetime:latest' }"

    input:
    tuple val(meta), path(tree), path(fasta), path(metadata)
    val outgroup_ids
    val clock_filter
    val clock_filter_method
    val clock_model
    val remove_future_tips
    val seed

    output:
    tuple val(meta), path("temporal_tree.nwk"),          emit: temporal_tree
    tuple val(meta), path("genetic_distance_tree.nwk"), emit: genetic_distance_tree
    tuple val(meta), path("final.alignment.fasta"),     emit: fasta
    tuple val(meta), path("final.metadata.tsv"),        emit: metadata
    tuple val(meta), path("final.dates.tsv"),           emit: dates
    tuple val(meta), path("clock_model_stats.yml"),     emit: clock_stats
    tuple val(meta), path("all_outliers.tsv"),          emit: outliers
    tuple val(meta), path("root_to_tip.png"),           emit: rtt_plot
    tuple val(meta), path("timetree.png"),              emit: timetree_plot
    tuple val(meta), path("future_tips_warning.txt"),   emit: future_tips_warning, optional: true
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outgroup_list = outgroup_ids instanceof List ? outgroup_ids : (outgroup_ids ? [outgroup_ids] : [])
    def reroot_args   = outgroup_list ? "--reroot ${outgroup_list.join(' ')}" : ""
    def seed_arg      = seed ? "--seed ${seed}" : ""
    def cf            = clock_filter ?: 3.0
    def cf_method     = clock_filter_method ?: 'local'
    def model_arg     = clock_model ?: 'strict'
    def future_arg    = (remove_future_tips != null && remove_future_tips.toString() != 'null' && remove_future_tips.toString() != '-1') ? "--remove-future-tips ${remove_future_tips}" : "--remove-future-tips -1"
    def remove_root_arg = outgroup_list ? "--remove-root" : "--no-remove-root"
    def sample_id_field  = params.sample_id_field ?: (params.fasta ? 'strain' : 'accessionVersion')
    def date_field       = params.collection_date_field ?: (params.fasta ? 'date' : 'sampleCollectionDate')
    def sample_id_arg    = "--sample-id-field ${sample_id_field}"
    def date_field_arg   = "--date-field ${date_field}"
    """
    python3 ${moduleDir}/resources/run_iterative_treetime.py \\
        --tree ${tree} \\
        --aln ${fasta} \\
        --dates ${metadata} \\
        --outdir . \\
        --prefix '' \\
        ${reroot_args} \\
        --clock-filter ${cf} \\
        --clock-filter-method ${cf_method} \\
        --clock-model ${model_arg} \\
        ${future_arg} \\
        ${remove_root_arg} \\
        ${sample_id_arg} \\
        ${date_field_arg} \\
        ${seed_arg}

    # Write the genetic distance tree (tree BEFORE branch_length_to_years)
    # The last pruned tree from the iteration is the genetic distance tree
    # It's already written by the script as part of iteration output
    # We need to save it separately — the script outputs temporal tree only
    # So we use the last iteration's pruned tree or the input tree if no iterations needed
    if [ -f *_iter*.nwk ]; then
        cp \$(ls -t *_iter*.nwk | head -1) genetic_distance_tree.nwk
    else
        cp ${tree} genetic_distance_tree.nwk
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treetime: \$(python -c "from importlib.metadata import version; print(version('treetime'))")
        biopython: \$(python -c "from importlib.metadata import version; print(version('biopython'))")
        matplotlib: \$(python -c "from importlib.metadata import version; print(version('matplotlib'))")
        pandas: \$(python -c "from importlib.metadata import version; print(version('pandas'))")
    END_VERSIONS
    """

    stub:
    """
    touch temporal_tree.nwk
    touch genetic_distance_tree.nwk
    touch final.alignment.fasta
    touch final.metadata.tsv
    touch final.dates.tsv
    touch clock_model_stats.yml
    touch all_outliers.tsv
    touch root_to_tip.png
    touch timetree.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treetime: \$(python -c "from importlib.metadata import version; print(version('treetime'))")
        biopython: \$(python -c "from importlib.metadata import version; print(version('biopython'))")
        matplotlib: \$(python -c "from importlib.metadata import version; print(version('matplotlib'))")
        pandas: \$(python -c "from importlib.metadata import version; print(version('pandas'))")
    END_VERSIONS
    """
}
