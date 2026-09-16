process COMBINE_REMOVED_SEQUENCES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas:2.2.1' :
        'community.wave.seqera.io/library/pandas:2.2.1' }"

    input:
    tuple val(meta), path(excluded_samples), path(all_outliers)

    output:
    tuple val(meta), path("*_removed_sequences.csv"), emit: csv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${moduleDir}/resources/combine_removed_sequences.py \\
        --excluded ${excluded_samples} \\
        --outliers ${all_outliers} \\
        -o ${prefix}_removed_sequences.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_removed_sequences.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
