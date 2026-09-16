process DATA_AND_RUN_STATS_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beast_pype:0.7.1' :
        'community.wave.seqera.io/library/beast_pype:0.7.1' }"

    input:
    tuple val(meta), path(raw_metadata), path(cleaned_metadata), path(filtered_metadata), path(pathoplexus_timestamp), val(beast_out_models), path(beast_out_files)

    output:
    tuple val(meta), path("*_data_and_run_stats_report.html") , emit: html
    tuple val(meta), path("*_data_and_run_stats_report.ipynb"), emit: notebook
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    prefix                    = task.ext.prefix ?: "${meta.id}"
    def kernel_name           = task.ext.kernel_name ? "-k ${task.ext.kernel_name}" : '-k ebola_flow'
    def collection_date_field = task.ext.collection_date_field ? "-c '${task.ext.collection_date_field}'" : ''
    // beast_out_models and beast_out_files are aligned lists (one entry per
    // BEAST .out file). Group the .out files back per model to build one
    // --beast_model / --beast_out pair per model.
    def model_list            = beast_out_models instanceof List ? beast_out_models : [beast_out_models]
    def out_list              = beast_out_files  instanceof List ? beast_out_files  : [beast_out_files]
    def grouped               = [:]
    model_list.eachWithIndex { model, i ->
        grouped.get(model, []) << out_list[i].name
    }
    def beast_args            = grouped.collect { model, files ->
        "--beast_model '${model}' --beast_out '${files.join(',')}'"
    }.join(' ')
    """
    # Register the ipykernel in a task-local directory to avoid race conditions on SLURM
    export JUPYTER_DATA_DIR="\$(pwd)/.local/share/jupyter"
    python -m ipykernel install --prefix="\$(pwd)/.local" --name ebola_flow --display-name "ebola_flow"

    # Build the data and run stats report directly from the inlined
    # gen_data_and_run_stats_report source (see resources/).
    python3 ${moduleDir}/resources/gen_data_and_run_stats_report.py \\
        -o ${prefix}_data_and_run_stats_report.ipynb \\
        -r ${raw_metadata} \\
        -m ${cleaned_metadata} \\
        -f ${filtered_metadata} \\
        -t ${pathoplexus_timestamp} \\
        ${beast_args} \\
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
    touch ${prefix}_data_and_run_stats_report.html
    touch ${prefix}_data_and_run_stats_report.ipynb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast_pype: \$(pip show beast_pype | grep '^Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}
