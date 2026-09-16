process TEMPEST_METADATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(metadata_tsv)
    val(reference_id)
    val(reference_date)
    val(include_reference)
    path(outgroup_dates)

    output:
    tuple val(meta), path("*_metadata_for_TempEst.tsv"), emit: tempest_metadata
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args              = task.ext.args ?: [:]
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def sample_id_field   = args.sample_id_field ?: 'accessionVersion'
    def date_field        = args.date_field ?: 'sampleCollectionDate'
    def include_ref       = include_reference ? 'True' : 'False'
    def has_outgroup      = outgroup_dates.name != 'NO_FILE' ? 'True' : 'False'
    """
    #!/usr/bin/env python3

    import re
    import platform
    from importlib.metadata import version as _pkg_version

    import pandas as pd

    metadata_tsv = '${metadata_tsv}'
    sample_id_field = '${sample_id_field}'
    date_field = '${date_field}'
    reference_id = '${reference_id}'
    reference_date_raw = '${reference_date}'
    include_ref = ${include_ref}
    has_outgroup = ${has_outgroup}
    prefix = '${prefix}'

    df = pd.read_csv(metadata_tsv, sep='\\t')

    # Select name and date columns
    tempest_df = df[[sample_id_field, date_field]].copy()
    tempest_df.columns = ['name', 'date']

    # Process reference date - if not YYYY-MM-DD, round to mid-month
    if reference_date_raw:
        # Match YYYY-MM-DD
        if re.fullmatch(r'\\d{4}-\\d{2}-\\d{2}', reference_date_raw):
            ref_date = reference_date_raw
        # Match YYYY-MM (round to 15th)
        elif re.fullmatch(r'\\d{4}-\\d{2}', reference_date_raw):
            ref_date = reference_date_raw + '-15'
        # Match YYYY (round to July 1st)
        elif re.fullmatch(r'\\d{4}', reference_date_raw):
            ref_date = reference_date_raw + '-07-01'
        else:
            ref_date = reference_date_raw
    else:
        ref_date = ''

    # Add reference row if requested
    if include_ref and reference_id and ref_date:
        ref_row = pd.DataFrame({'name': [reference_id], 'date': [ref_date]})
        tempest_df = pd.concat([tempest_df, ref_row], ignore_index=True)

    # Add outgroup dates if provided
    if has_outgroup:
        outgroup_df = pd.read_csv('${outgroup_dates}', sep='\\t')
        tempest_df = pd.concat([tempest_df, outgroup_df], ignore_index=True)

    tempest_df.to_csv(f'{prefix}_metadata_for_TempEst.tsv', sep='\\t', index=False)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {platform.python_version()}\\n')
        f.write(f'    pandas: {_pkg_version("pandas")}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metadata_for_TempEst.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$(python -c "from importlib.metadata import version; print(version('pandas'))")
    END_VERSIONS
    """
}
