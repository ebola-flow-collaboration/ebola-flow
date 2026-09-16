process BEAST_PYPE_SUMMARY_TREE_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.0' :
        'community.wave.seqera.io/library/beast_pype:0.7.0' }"

    input:
    tuple val(meta), path(summary_tree), path(beast_xml)

    output:
    tuple val(meta), path("*_summary_tree_report.html") , emit: html
    tuple val(meta), path("*_summary_tree_report.ipynb"), emit: notebook
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    prefix                    = task.ext.prefix ?: "${meta.id}"
    def kernel_name           = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow_R'
    def collection_date_field = task.ext.collection_date_field ? "-c '${task.ext.collection_date_field}'" : ''
    def plot_width            = task.ext.plot_width ? "--plot_width ${task.ext.plot_width}" : ''
    def plot_height           = task.ext.plot_height ? "--plot_height ${task.ext.plot_height}" : ''
    def plot_res              = task.ext.plot_res ? "--plot_res ${task.ext.plot_res}" : ''
    """
    # Register the R kernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    Rscript -e "IRkernel::installspec(name = 'ebola_flow_R', displayname = 'ebola_flow_R', prefix = '\$(pwd)/.local')"

    beast_pype \\
        summary-tree-report \\
        -s ${summary_tree} \\
        -x ${beast_xml} \\
        -o ${prefix}_summary_tree_report.ipynb \\
        ${kernel_name} \\
        ${collection_date_field} \\
        ${plot_width} \\
        ${plot_height} \\
        ${plot_res} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary_tree_report.html
    touch ${prefix}_summary_tree_report.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}
