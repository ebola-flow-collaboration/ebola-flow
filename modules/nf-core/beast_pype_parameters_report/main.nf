process BEAST_PYPE_PARAMETERS_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.1' :
        'community.wave.seqera.io/library/beast_pype:0.7.1' }"

    input:
    tuple val(meta), path(merged_logs), path(beast_xmls)

    output:
    tuple val(meta), path("*_parameters_report.html") , emit: html
    tuple val(meta), path("*_parameters_report.ipynb"), emit: notebook
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    prefix              = task.ext.prefix ?: "${meta.id}"
    def log_files       = merged_logs instanceof List ? merged_logs : [merged_logs]
    def xml_files       = beast_xmls instanceof List ? beast_xmls : [beast_xmls]
    def is_comparative  = xml_files.size() > 1
    def report_template = task.ext.report_template ?: (is_comparative ? 'Generic-Comparative' : 'Generic')
    def kernel_name     = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow'
    def xml_set_label   = task.ext.xml_set_label ? "--xml_set_label '${task.ext.xml_set_label}'" : (is_comparative ? "--xml_set_label 'model'" : '')
    def log_args        = log_files.collect { "-l ${it}" }.join(' ')
    def xml_args        = xml_files.collect { "-x ${it}" }.join(' ')
    """
    # Register the ipykernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    python -m ipykernel install --prefix="\$(pwd)/.local" --name ebola_flow --display-name "ebola_flow"

    beast_pype \\
        parameters-report \\
        -o ${prefix}_parameters_report.ipynb \\
        ${log_args} \\
        ${xml_args} \\
        ${kernel_name} \\
        ${xml_set_label} \\
        ${args} \\
        ${report_template}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_parameters_report.html
    touch ${prefix}_parameters_report.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}
