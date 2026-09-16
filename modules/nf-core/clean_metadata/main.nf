process CLEAN_METADATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(metadata_json), path(exclude_ids)

    output:
    tuple val(meta), path("*_metadata_raw.tsv")   , emit: metadata_raw
    tuple val(meta), path("*_metadata_full.tsv")  , emit: metadata_full
    tuple val(meta), path("*_metadata_trimmed.tsv")       , emit: metadata
    tuple val(meta), path("*_ids_to_remove.txt")  , emit: ids_to_remove
    tuple val(meta), path("*_excluded_samples.csv"), emit: excluded_samples
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args              = task.ext.args ?: [:]
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def sample_id_field   = args.sample_id_field ?: 'accessionVersion'
    def date_field        = args.date_field ?: 'sampleCollectionDate'
    def filter_dates      = args.filter_invalid_dates ?: false
    """
    #!/usr/bin/env python

    import re
    import json
    import pandas as pd
    from pathlib import Path
    from importlib.metadata import version as _pkg_version

    # --- Convert JSON to DataFrame ---
    with open('${metadata_json}') as f:
        data = json.load(f)
    df = pd.DataFrame(data['data'])

    # --- Write raw (unfiltered) metadata straight from the Pathoplexus pull ---
    df.to_csv('${prefix}_metadata_raw.tsv', sep='\\t', index=False)

    ids_to_remove = set()

    # Provenance records for each excluded sample. Maps accessionVersion (or the
    # raw accession when it can't be resolved) -> dict with reason/source fields.
    excluded_records = {}

    def record_exclusion(seq_id, accession, reason, source):
        rec = excluded_records.get(seq_id)
        if rec is None:
            excluded_records[seq_id] = {
                'accessionVersion': seq_id,
                'accession': accession,
                'reason_for_removal': reason,
                'source_for_reason': source,
            }
        else:
            if reason and reason not in rec['reason_for_removal'].split('; '):
                rec['reason_for_removal'] = '; '.join(
                    [r for r in [rec['reason_for_removal'], reason] if r]
                )
            if source and source not in rec['source_for_reason'].split('; '):
                rec['source_for_reason'] = '; '.join(
                    [s for s in [rec['source_for_reason'], source] if s]
                )

    # --- Filter by exclude list ---
    # The exclude list is a CSV with a header and the columns:
    #   accession, reason_for_removal, source_for_reason
    # The 'accession' column supports both accessionVersion (e.g. PP_006X5RB.1)
    # and accession (e.g. PP_006X5RB). If an accession without a version is
    # given, it is looked up in the metadata to resolve the full accessionVersion.
    exclude_file = Path('${exclude_ids}')
    if exclude_file.exists() and exclude_file.stat().st_size > 0:
        accession_to_version = {}
        if 'accession' in df.columns and '${sample_id_field}' in df.columns:
            accession_to_version = dict(zip(df['accession'].astype(str), df['${sample_id_field}'].astype(str)))

        exclude_df = pd.read_csv(exclude_file, dtype=str).fillna('')
        if 'accession' not in exclude_df.columns:
            raise ValueError(
                "Exclude list '${exclude_ids}' must be a CSV with an 'accession' "
                "column (expected header: "
                "accession,reason_for_removal,source_for_reason)."
            )

        for _, ex_row in exclude_df.iterrows():
            accession = str(ex_row['accession']).strip()
            reason = str(ex_row.get('reason_for_removal', '')).strip()
            source = str(ex_row.get('source_for_reason', '')).strip()
            if accession:
                if accession in df['${sample_id_field}'].values:
                    seq_id = accession
                elif accession in accession_to_version:
                    seq_id = accession_to_version[accession]
                else:
                    seq_id = accession
                ids_to_remove.add(seq_id)
                record_exclusion(seq_id, accession, reason, source)

    # --- Filter invalid dates (not YYYY-MM-DD) ---
    if ${filter_dates ? 'True' : 'False'}:
        date_pattern = re.compile(r'^\\d{4}-\\d{2}-\\d{2}\$')
        invalid_date_mask = ~df['${date_field}'].astype(str).apply(
            lambda x: bool(date_pattern.match(x))
        )
        for _, dt_row in df.loc[invalid_date_mask].iterrows():
            seq_id = str(dt_row['${sample_id_field}'])
            accession = str(dt_row['accession']) if 'accession' in df.columns else ''
            bad_date = str(dt_row['${date_field}'])
            ids_to_remove.add(seq_id)
            record_exclusion(
                seq_id,
                accession,
                'Invalid collection date (not YYYY-MM-DD): ' + bad_date,
                'clean_metadata date filter',
            )

    # --- Apply filtering ---
    df_clean = df[~df['${sample_id_field}'].isin(ids_to_remove)]

    # --- Write full cleaned metadata ---
    df_clean.to_csv('${prefix}_metadata_full.tsv', sep='\\t', index=False)

    # --- Write trimmed cleaned metadata (subset of fields for downstream use) ---
    trimmed_fields = [
        'accessionVersion',
        'sampleCollectionDate',
        'geoLocAdmin1',
        'geoLocAdmin2',
        'geoLocCity',
        'geoLocCountry',
        'hostAge',
        'submittedAtTimestamp',
    ]
    # Only keep fields that are present in the metadata
    df_trimmed = df_clean[[f for f in trimmed_fields if f in df_clean.columns]]
    df_trimmed.to_csv('${prefix}_metadata_trimmed.tsv', sep='\\t', index=False)

    # --- Write IDs to remove (for seqkit grep) ---
    with open('${prefix}_ids_to_remove.txt', 'w') as f:
        for seq_id in sorted(ids_to_remove):
            f.write(seq_id + '\\n')

    # --- Write excluded samples provenance CSV ---
    # Combines exclude-list entries and invalid-date removals, with the reason
    # for removal (and its source) for each excluded sample.
    excluded_cols = [
        'accessionVersion',
        'accession',
        'reason_for_removal',
        'source_for_reason',
    ]
    excluded_out = pd.DataFrame(
        list(excluded_records.values()), columns=excluded_cols
    ).sort_values('accessionVersion')
    excluded_out.to_csv('${prefix}_excluded_samples.csv', index=False)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'  python: {__import__("platform").python_version()}\\n')
        f.write(f'  pandas: {_pkg_version("pandas")}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metadata_raw.tsv
    touch ${prefix}_metadata_full.tsv
    touch ${prefix}_metadata_trimmed.tsv
    touch ${prefix}_ids_to_remove.txt
    touch ${prefix}_excluded_samples.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$(python -c "from importlib.metadata import version; print(version('pandas'))")
    END_VERSIONS
    """
}
