process METADATA_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.1' :
        'community.wave.seqera.io/library/beast_pype:0.7.1' }"

    input:
    tuple val(meta), path(raw_metadata), path(cleaned_metadata), path(filtered_metadata), path(pathoplexus_timestamp)

    output:
    tuple val(meta), path("*_metadata_report.html") , emit: html
    tuple val(meta), path("*_metadata_report.ipynb"), emit: notebook
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    prefix                    = task.ext.prefix ?: "${meta.id}"
    def kernel_name           = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow'
    def collection_date_field = task.ext.collection_date_field ? "-c '${task.ext.collection_date_field}'" : ''
    """
    # Register the ipykernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    python -m ipykernel install --prefix="\$(pwd)/.local" --name ebola_flow --display-name "ebola_flow"

    # Build the metadata report directly from the inlined
    # gen_metadata_report source (see resources/).
    python3 ${moduleDir}/resources/gen_metadata_report.py \\
        -o ${prefix}_metadata_report.ipynb \\
        -r ${raw_metadata} \\
        -m ${cleaned_metadata} \\
        -f ${filtered_metadata} \\
        -t ${pathoplexus_timestamp} \\
        ${kernel_name} \\
        ${collection_date_field} \\
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
