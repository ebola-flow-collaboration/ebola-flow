process CLEAN_OUTGROUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(metadata_json)
    val(ref_id)

    output:
    tuple val(meta), path("*_ids_to_remove.txt")   , emit: ids_to_remove
    tuple val(meta), path("*_outgroup_dates.tsv")  , emit: dates
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def date_field = task.ext.args?.date_field ?: 'sampleCollectionDate'
    def id_field = task.ext.args?.id_field ?: 'accessionVersion'
    """
    #!/usr/bin/env python3

    import json
    import re
    import platform

    ref_id = '${ref_id}'
    ref_base = ref_id.split('.')[0]  # e.g. NC_014373

    with open('${metadata_json}') as f:
        data = json.load(f)

    ids_to_remove = []
    outgroup_dates = []

    for record in data.get('data', []):
        insdc_full = record.get('insdcAccessionFull', '') or ''
        insdc_base = record.get('insdcAccessionBase', '') or ''
        accession = record.get('${id_field}', '')
        collection_date = record.get('${date_field}', '') or ''

        # Check if this record is the reference genome
        if ref_base in insdc_full or ref_base in insdc_base:
            ids_to_remove.append(accession)
            continue

        # Round date to mid-month/mid-year if incomplete
        if re.match(r'^\\d{4}-\\d{2}-\\d{2}\$', collection_date):
            rounded_date = collection_date
        elif re.match(r'^\\d{4}-\\d{2}\$', collection_date):
            rounded_date = collection_date + '-15'
        elif re.match(r'^\\d{4}\$', collection_date):
            rounded_date = collection_date + '-07-01'
        else:
            continue  # skip sequences with no valid date

        outgroup_dates.append((accession, rounded_date))

    # Write IDs to remove (reference's Pathoplexus accessionVersion)
    with open('${prefix}_ids_to_remove.txt', 'w') as f:
        for seq_id in ids_to_remove:
            f.write(seq_id + '\\n')

    # Write outgroup dates TSV (for TempEst metadata)
    with open('${prefix}_outgroup_dates.tsv', 'w') as f:
        f.write('name\\tdate\\n')
        for name, d in outgroup_dates:
            f.write(f'{name}\\t{d}\\n')

    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {platform.python_version()}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ids_to_remove.txt
    touch ${prefix}_outgroup_dates.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
