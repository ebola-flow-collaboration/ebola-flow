process PATHOPLEXUS_FETCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.fasta")            , emit: fasta
    tuple val(meta), path("metadata.json")      , emit: metadata
    tuple val(meta), path("fetch_timestamp.txt"), emit: timestamp
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: [:]
    def prefix                = task.ext.prefix ?: "${meta.id}"
    def organism              = args.organism ?: meta.organism ?: 'ebola-bdbv'
    def date_from             = (args.date_from && args.date_from != 'null') ? args.date_from : ((meta.date_from && meta.date_from != 'null') ? meta.date_from : '')
    def date_to               = (args.date_to && args.date_to != 'null') ? args.date_to : ((meta.date_to && meta.date_to != 'null') ? meta.date_to : '')
    def submitted_date        = (args.submitted_date && args.submitted_date != 'null') ? args.submitted_date : ((meta.submitted_date && meta.submitted_date != 'null') ? meta.submitted_date : '')
    def length_from           = (args.length_from && args.length_from != 'null') ? args.length_from : ((meta.length_from && meta.length_from != 'null') ? meta.length_from : '')
    def date_from_param       = date_from ? "sampleCollectionDateRangeLowerFrom=${date_from}&" : ''
    def date_to_param         = date_to ? "sampleCollectionDateRangeUpperTo=${date_to}&" : ''
    def submitted_date_param  = submitted_date ? "submittedDate=${submitted_date}&" : ''
    def length_from_param     = length_from ? "lengthFrom=${length_from}&" : ''
    """
    # Record the date and time this Pathoplexus data was fetched
    date +%s > fetch_timestamp.txt

    # Download unaligned nucleotide sequences
    curl -f -sS \\
        "https://lapis.pathoplexus.org/${organism}/sample/unalignedNucleotideSequences?${date_from_param}${date_to_param}${submitted_date_param}${length_from_param}versionStatus=LATEST_VERSION" \\
        -o ${prefix}.fasta

    # Download metadata as JSON
    curl -f -sS \\
        "https://lapis.pathoplexus.org/${organism}/sample/details?${date_from_param}${date_to_param}${submitted_date_param}${length_from_param}versionStatus=LATEST_VERSION" \\
        -o metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    date +%s > fetch_timestamp.txt
    touch ${prefix}.fasta
    echo '{"data":[]}' > metadata.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
}
