process BEAST2_XML_GEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-beast2xml:latest' :
        'quay.io/biocontainers/mulled-v2-beast2xml:latest' }"

    input:
    tuple val(meta), path(template_xml), path(fasta), path(metadata), path(initial_tree)

    output:
    tuple val(meta), path("*.xml"), emit: xml
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: [:]
    def prefix                = task.ext.prefix ?: "${meta.id}"
    def sample_id_field       = args.sample_id_field ?: 'strain'
    def collection_date_field = args.collection_date_field ?: 'date'
    def log_file_basename     = "log_file_basename='${args.log_file_basename ?: ''}',"
    def chain_length          = args.chain_length ? "chain_length=${args.chain_length}," : ''
    def trace_log_every       = args.trace_log_every ? "trace_log_every=${args.trace_log_every}," : ''
    def tree_log_every        = args.tree_log_every ? "tree_log_every=${args.tree_log_every}," : ''
    def screen_log_every      = args.screen_log_every ? "screen_log_every=${args.screen_log_every}," : ''
    def store_state_every     = args.store_state_every ? "store_state_every=${args.store_state_every}," : ''
    def use_initial_tree      = initial_tree.name != 'NO_FILE' ? 'True' : 'False'
    """
    #!/usr/bin/env python

    from dark.fasta import FastaReads
    from beast2xml.beast2 import BEAST2XML
    from importlib.metadata import version as _pkg_version

    template_xml_path = '${template_xml}'
    fasta_path = '${fasta}'
    metadata_path = '${metadata}'
    output_path = '${prefix}_generated.xml'

    if metadata_path.endswith('.tsv'):
        delimiter = '\\t'
    elif metadata_path.endswith('.csv'):
        delimiter = ','
    else:
        raise TypeError(
            f"metadata_path must be a csv or tsv file. Value given is {metadata_path}")

    beast2xml = BEAST2XML(template=template_xml_path)
    seqs = FastaReads([fasta_path])
    beast2xml.add_sequences(seqs)
    beast2xml.add_dates(
        date_data=metadata_path,
        seperator=delimiter,
        sample_id_field='${sample_id_field}',
        collection_date_field='${collection_date_field}',
    )

    if ${use_initial_tree}:
        beast2xml.add_initial_tree('${initial_tree}')

    beast2xml.to_xml(
        output_path,
        ${log_file_basename}
        ${chain_length}
        ${trace_log_every}
        ${tree_log_every}
        ${screen_log_every}
        ${store_state_every}
    )

    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"BEAST2_XML_GEN":\\n')
        f.write(f'  dark-matter: {_pkg_version("dark-matter")}\\n')
        f.write(f'  beast2-xml: {_pkg_version("beast2-xml")}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_generated.xml
    cat <<-END_VERSIONS > versions.yml
    "BEAST2_XML_GEN":
        dark-matter: \$(python -c "from importlib.metadata import version; print(version('dark-matter'))")
        beast2-xml: \$(python -c "from importlib.metadata import version; print(version('beast2-xml'))")
    END_VERSIONS
    """
}
