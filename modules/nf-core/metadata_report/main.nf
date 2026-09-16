process METADATA_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.1' :
        'community.wave.seqera.io/library/beast_pype:0.7.1' }"

    input:
    tuple val(meta), path(metadata)

    output:
    tuple val(meta), path("*_metadata_report.html") , emit: html
    tuple val(meta), path("*_metadata_report.ipynb"), emit: notebook
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    prefix                    = task.ext.prefix ?: "${meta.id}"
    def meta_files            = metadata instanceof List ? metadata : [metadata]
    def is_comparative        = meta_files.size() > 1
    def kernel_name           = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow'
    def collection_date_field = task.ext.collection_date_field ? "-c '${task.ext.collection_date_field}'" : ''
    def xml_set_label         = task.ext.xml_set_label ? "--xml_set_label '${task.ext.xml_set_label}'" : ''
    def meta_args             = meta_files.collect { "-m ${it}" }.join(' ')
    def set_name_args         = is_comparative ? meta_files.withIndex().collect { f, i -> "-n set_${i + 1}" }.join(' ') : ''
    """
    # Register the ipykernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    python -m ipykernel install --prefix="\$(pwd)/.local" --name ebola_flow --display-name "ebola_flow"

    # Build the metadata report directly from the inlined gen_metadata_report source
    # (see resources/gen_metadata_report.py) rather than the beast_pype CLI command.
    python3 ${moduleDir}/resources/gen_metadata_report.py \\
        -o ${prefix}_metadata_report.ipynb \\
        ${meta_args} \\
        ${set_name_args} \\
        ${kernel_name} \\
        ${collection_date_field} \\
        ${xml_set_label} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metadata_report.html
    touch ${prefix}_metadata_report.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}
