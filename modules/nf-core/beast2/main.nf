process BEAST2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beast2:2.7.7--haf2fa61_0' :
        'quay.io/biocontainers/beast2:2.7.7--haf2fa61_0' }"

    input:
    tuple val(meta), path(xml)

    output:
    tuple val(meta), path("*.log")       , emit: log
    tuple val(meta), path("*.trees")     , emit: trees     , optional: true
    tuple val(meta), path("*.state.xml") , emit: state     , optional: true
    tuple val(meta), path("*.xml.state") , emit: xml_state , optional: true
    tuple val(meta), path("*.out")       , emit: screen_log
    tuple val("${task.process}"), val('beast2'), eval("beast -version 2>&1 | head -n1 | sed 's/.*v//'"), emit: versions_beast2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seed   = task.ext.seed   ?: '1'
    def packages = task.ext.packages ?: ''
    def install_packages = packages ? "packagemanager -add ${packages} || true" : ''
    """
    ${install_packages}

    beast \\
        -threads $task.cpus \\
        -seed $seed \\
        -prefix $prefix \\
        -statefile ${prefix}.state.xml \\
        $xml \\
        | tee ${prefix}.out
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.log
    touch ${prefix}.trees
    touch ${prefix}.state.xml
    touch ${prefix}.out
    """
}