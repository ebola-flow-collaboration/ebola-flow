/*
 * EMERGENCE_ANALYSIS
 *
 * Runs the variant emergence surveillance analytics (root-to-tip regression,
 * clock-outlier detection, and highlighted tree plots) as a parameterised
 * R Markdown report. Emits both the raw .Rmd source and the rendered HTML.
 */

process EMERGENCE_ANALYSIS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-treetime:latest' :
        'quay.io/biocontainers/mulled-v2-treetime:latest' }"

    input:
    tuple val(meta), path(div_tree_raw), path(div_tree_clean), path(time_tree), path(metadata)

    output:
    tuple val(meta), path("*_emergence_analysis.html"), emit: html
    tuple val(meta), path("*_emergence_analysis.Rmd") , emit: rmd
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def name_col  = task.ext.name_col ?: (params.sample_id_field ?: 'accessionVersion')
    def date_col  = task.ext.date_col ?: (params.collection_date_field ?: 'sampleCollectionDate')
    def y_cut_val = task.ext.y_cut_val != null ? task.ext.y_cut_val : 0.00081
    def n_sd      = task.ext.n_sd != null ? task.ext.n_sd : 3
    """
    # Copy the report source and RTT helper into the work dir so paths resolve
    cp ${moduleDir}/resources/emergence_analysis.Rmd ${prefix}_emergence_analysis.Rmd
    cp ${moduleDir}/resources/rtt_fit_plot.R rtt_fit_plot.R

    Rscript -e "rmarkdown::render(
        '${prefix}_emergence_analysis.Rmd',
        output_file = '${prefix}_emergence_analysis.html',
        params = list(
            metadata       = '${metadata}',
            div_tree_raw   = '${div_tree_raw}',
            div_tree_clean = '${div_tree_clean}',
            time_tree      = '${time_tree}',
            rtt_script     = 'rtt_fit_plot.R',
            name_col       = '${name_col}',
            date_col       = '${date_col}',
            y_cut_val      = ${y_cut_val},
            n_sd           = ${n_sd}
        ),
        knit_root_dir = getwd()
    )" ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
        ggtree: \$(Rscript -e "cat(as.character(packageVersion('ggtree')))")
        ape: \$(Rscript -e "cat(as.character(packageVersion('ape')))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_emergence_analysis.html
    touch ${prefix}_emergence_analysis.Rmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
        ggtree: \$(Rscript -e "cat(as.character(packageVersion('ggtree')))")
        ape: \$(Rscript -e "cat(as.character(packageVersion('ape')))")
    END_VERSIONS
    """
}