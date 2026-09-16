process DOWNSAMPLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-treetime:latest' :
        'quay.io/biocontainers/mulled-v2-treetime:latest' }"

    input:
    tuple val(meta), path(tree), path(fasta), path(metadata)
    path clock_stats_yml
    val sample_size
    val seed

    output:
    tuple val(meta), path("*_downsampled_sequences.fasta") , emit: fasta
    tuple val(meta), path("*_downsampled_metadata.csv")    , emit: metadata
    tuple val(meta), path("*_downsampled_tree.nwk")        , emit: tree
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: [:]
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def sample_id_field = args.sample_id_field ?: 'accessionVersion'
    def seed_val        = seed != null && seed != '' ? seed : 'None'
    """
    #!/usr/bin/env python3

    from Bio import Phylo, SeqIO
    import pandas as pd
    import numpy as np
    import yaml
    from copy import deepcopy
    from importlib.metadata import version as _pkg_version

    prefix = '${prefix}'
    fasta_path = '${fasta}'
    metadata_path = '${metadata}'
    tree_path = '${tree}'
    sample_id_field = '${sample_id_field}'
    sample_size = ${sample_size}
    seed = ${seed_val}

    # Read tree
    tree = Phylo.read(tree_path, 'newick')
    tips = tree.get_terminals()
    tip_names = [t.name for t in tips]
    n_tips = len(tips)

    if sample_size >= n_tips:
        # No downsampling needed — copy inputs as outputs
        import shutil
        shutil.copy(fasta_path, f'{prefix}_downsampled_sequences.fasta')
        sep = '\\t' if metadata_path.endswith('.tsv') else ','
        metadata_df = pd.read_csv(metadata_path, sep=sep)
        metadata_df.to_csv(f'{prefix}_downsampled_metadata.csv', index=False)
        Phylo.write(tree, f'{prefix}_downsampled_tree.nwk', 'newick')
    else:
        dist2root = np.array([tree.distance(tree.root, t) for t in tips])

        # Read clock model parameters from YAML
        with open('${clock_stats_yml}') as f:
            clock_stats = yaml.safe_load(f)

        slope = clock_stats.get('clock_rate', None)
        t_mrca = clock_stats.get('TMRCA_decimal_year', None)
        if slope is None:
            raise ValueError("clock_model_stats.yml missing 'clock_rate'")

        intercept = clock_stats.get('intercept', None)
        if intercept is None and t_mrca is not None:
            intercept = -slope * t_mrca
        if intercept is None:
            raise ValueError("Cannot determine intercept from clock_model_stats.yml")

        # Get decimal year dates from metadata
        import csv
        from datetime import datetime
        with open(metadata_path) as f:
            sample = f.read(4096)
        dialect = csv.Sniffer().sniff(sample)
        sep = dialect.delimiter
        metadata_df = pd.read_csv(metadata_path, sep=sep)

        date_candidates = [c for c in metadata_df.columns if 'date' in c.lower()]
        date_col = date_candidates[0] if date_candidates else metadata_df.columns[1]

        def to_decimal_year(d):
            try:
                dt = pd.to_datetime(d)
                year_start = datetime(dt.year, 1, 1)
                year_end = datetime(dt.year + 1, 1, 1)
                return dt.year + (dt - year_start).total_seconds() / (year_end - year_start).total_seconds()
            except Exception:
                return np.nan

        metadata_df['_decimal_year'] = metadata_df[date_col].apply(to_decimal_year)
        date_lookup = dict(zip(metadata_df[sample_id_field].astype(str), metadata_df['_decimal_year']))

        # Only use tips that have valid dates
        valid_data = [(i, n) for i, n in enumerate(tip_names) if not np.isnan(date_lookup.get(n, np.nan))]
        valid_indices = [i for i, n in valid_data]
        valid_names = [n for i, n in valid_data]
        valid_dist = dist2root[valid_indices]
        x_year_decimal = np.array([date_lookup[n] for n in valid_names])

        # Compute residuals (same logic as temporal_pruning_sampler)
        root_to_tip_expected = slope * x_year_decimal + intercept
        abs_residuals = np.absolute(valid_dist - root_to_tip_expected)
        prune_prob = abs_residuals / abs_residuals.sum()

        # Probabilistically remove tips proportional to their residual
        to_prune = len(valid_names) - sample_size
        rng = np.random.default_rng(seed=seed)
        to_prune_indices = rng.choice(len(valid_names), size=to_prune, p=prune_prob, replace=False)
        sampled_ids = set(
            name for i, name in enumerate(valid_names) if i not in to_prune_indices
        )

        # Write downsampled outputs
        sep = '\\t' if metadata_path.endswith('.tsv') else ','
        metadata_df = pd.read_csv(metadata_path, sep=sep)
        selected_metadata = metadata_df[metadata_df[sample_id_field].isin(sampled_ids)]
        selected_seqs = [rec for rec in SeqIO.parse(fasta_path, 'fasta') if rec.id in sampled_ids]

        selected_metadata.to_csv(f'{prefix}_downsampled_metadata.csv', index=False)
        with open(f'{prefix}_downsampled_sequences.fasta', 'w') as handle:
            SeqIO.write(selected_seqs, handle, 'fasta')

        # Prune tree to sampled tips
        ds_tree = deepcopy(tree)
        tips_to_remove = [t for t in ds_tree.get_terminals() if t.name not in sampled_ids]
        for t in tips_to_remove:
            ds_tree.prune(t)
        Phylo.write(ds_tree, f'{prefix}_downsampled_tree.nwk', 'newick')

    # Write versions
    with open('versions.yml', 'w') as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    biopython: {_pkg_version("biopython")}\\n')
        f.write(f'    numpy: {_pkg_version("numpy")}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_downsampled_sequences.fasta
    touch ${prefix}_downsampled_metadata.csv
    touch ${prefix}_downsampled_tree.nwk
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "from importlib.metadata import version; print(version('biopython'))")
        numpy: \$(python -c "from importlib.metadata import version; print(version('numpy'))")
    END_VERSIONS
    """
}
