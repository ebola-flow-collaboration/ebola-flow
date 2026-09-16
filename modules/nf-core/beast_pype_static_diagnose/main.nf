process BEAST_PYPE_STATIC_DIAGNOSE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.0' :
        'community.wave.seqera.io/library/beast_pype:0.7.0' }"

    input:
    tuple val(meta), path(beast_outputs, stageAs: 'beast_outputs/*')

    output:
    tuple val(meta), path("${prefix}_*")            , emit: results
    tuple val(meta), path("${prefix}_*.html")       , emit: notebook, optional: true
    tuple val(meta), path("${prefix}_*.csv")        , emit: merged_log, optional: true
    tuple val(meta), path("${prefix}_*.trees")      , emit: merged_trees, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    prefix                 = task.ext.prefix ?: "${meta.id}"
    def burnin             = task.ext.burnin ? "-b ${task.ext.burnin}" : ''
    def params_per_section = task.ext.parameters_per_section ? "-n ${task.ext.parameters_per_section}" : ''
    def kernel_name        = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow'
    """
    # Register the ipykernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    python -m ipykernel install --prefix="\$(pwd)/.local" --name ebola_flow --display-name "ebola_flow"

    beast_pype \\
        static-diagnose-and-merge \\
        ${burnin} \\
        -o ${prefix}_ \\
        ${params_per_section} \\
        ${kernel_name} \\
        ${args} \\
        beast_outputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.log
    touch ${prefix}_merged.trees
    touch ${prefix}_diagnostics.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}
