#!/usr/bin/env python3
"""
run_iterative_treetime.py

Thin CLI wrapper around tree_time_scale.py for the ITERATIVE_TREETIME
Nextflow module. All TreeTime logic lives in tree_time_scale.py.

Outputs:
    - temporal_tree.nwk           : Final tree with branch lengths in years
    - genetic_distance_tree.nwk   : Final tree with branch lengths in subs/site
    - final.dates.tsv             : Dates for all retained tips
    - clock_model_stats.yml       : Clock model stats (rate, R^2, TMRCA, intercept)
    - final.alignment.fasta       : Final filtered alignment
    - final.metadata.tsv          : Final filtered metadata
    - all_outliers.tsv            : Combined outlier table across iterations
    - root_to_tip.png             : Root-to-tip regression plot
    - timetree.png                : Temporal tree plot
    - future_tips_warning.txt     : Warning about tips placed in future (if applicable)
"""

import argparse
import copy
import os
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from Bio import Phylo, SeqIO

# Local imports (same resources directory)
from tree_time_scale import (
    iterative_timescale,
    plot_root_to_tip,
    plot_temporal_tree,
)
from date_utilities import decimal_to_date


def parse_args():
    parser = argparse.ArgumentParser(
        description='Iterative TreeTime clock filtering and analysis'
    )
    parser.add_argument('--tree', required=True, help='Input newick tree')
    parser.add_argument('--aln', required=True, help='Input FASTA alignment')
    parser.add_argument('--dates', required=True, help='Input metadata (CSV or TSV)')
    parser.add_argument('--outdir', default='.', help='Output directory')
    parser.add_argument('--prefix', default='', help='Output file prefix')

    # TreeTime options
    parser.add_argument('--reroot', nargs='*', default=None,
                        help='Reroot method or tip names for rooting')
    parser.add_argument('--clock-filter', type=float, default=3.0,
                        help='Clock filter threshold (z-score or n_iqd)')
    parser.add_argument('--clock-filter-method', default='local',
                        choices=['local', 'residual'],
                        help='Clock filter method')
    parser.add_argument('--clock-model', default='strict',
                        choices=['strict', 'relaxed'],
                        help='Clock model type')
    parser.add_argument('--clock-rate', type=float, default=None,
                        help='Fixed clock rate')
    parser.add_argument('--clock-std', type=float, default=None,
                        help='Clock rate std (for relaxed clock)')
    parser.add_argument('--coalescent', default='opt',
                        help='Coalescent time constant')

    # Column names
    parser.add_argument('--sample-id-field', default=None,
                        help='Column name for sample IDs in dates file')
    parser.add_argument('--date-field', default='date',
                        help='Column name for dates in dates file')

    # Filtering options
    parser.add_argument('--remove-future-tips', type=float, default=0,
                        help='Buffer in days beyond youngest date for future tip removal. 0=no buffer. -1=disable.')
    parser.add_argument('--remove-root', action='store_true', default=True,
                        help='Remove root/outgroup sequences after filtering')
    parser.add_argument('--no-remove-root', dest='remove_root',
                        action='store_false')
    parser.add_argument('--max-iterations', type=int, default=50,
                        help='Maximum clock filter iterations')
    parser.add_argument('--negative-tolerance', type=float, default=0.001,
                        help='Tolerance for negative branch lengths')

    # Other
    parser.add_argument('--seed', type=int, default=None, help='Random seed')

    return parser.parse_args()


def fname(outdir, prefix, name):
    """Build output filepath, handling empty prefix."""
    if prefix:
        return os.path.join(outdir, f'{prefix}_{name}')
    return os.path.join(outdir, name)


def write_clock_stats(time_tree, filepath):
    """Write clock model statistics to YAML."""
    rate = time_tree.clock_model.get('slope', 0)
    intercept = time_tree.clock_model.get('intercept', 0)
    r_val = time_tree.clock_model.get('r_val', 0)
    r_sq = r_val ** 2 if r_val else 0
    t_mrca = -intercept / rate if rate and rate != 0 else None

    stats = {
        'clock_rate': float(rate) if rate else 0.0,
        'intercept': float(intercept) if intercept else 0.0,
        'r_squared': float(r_sq),
        'correlation_coefficient_r': float(r_val) if r_val else 0.0,
    }

    if t_mrca is not None:
        stats['TMRCA_decimal_year'] = float(t_mrca)
        tmrca_date = decimal_to_date(t_mrca)
        stats['TMRCA_date'] = tmrca_date.strftime('%Y-%m-%d')

    with open(filepath, 'w') as f:
        yaml.safe_dump(stats, f, default_flow_style=False)

    return stats


def write_future_tips_warning(time_tree, outliers_df, filepath, remove_future_tips):
    """Write future tips warning file (only if remove_future_tips is False)."""
    if remove_future_tips is not None and remove_future_tips >= 0:
        return

    tips = time_tree.tree.get_terminals()
    raw_dates = [t.raw_date_constraint for t in tips
                 if hasattr(t, 'raw_date_constraint')
                 and t.raw_date_constraint is not None]
    if not raw_dates:
        return

    youngest_date = max(raw_dates)
    tip_dates = {t.name: t.numdate for t in tips
                 if hasattr(t, 'numdate') and t.numdate is not None}
    future_tips = [(name, numdate, numdate - youngest_date)
                   for name, numdate in tip_dates.items()
                   if numdate > youngest_date]

    with open(filepath, 'w') as wf:
        if future_tips:
            wf.write(f'WARNING: The following tips are projected beyond the '
                     f'youngest sample date ({youngest_date:.4f}):\n')
            wf.write('This may indicate local rate variation or issues with '
                     'the molecular clock model.\n')
            wf.write('\nSUGGESTION: Consider using a relaxed clock model to '
                     'account for rate variation\namong lineages. Set '
                     '--treetime_clock_model relaxed in the pipeline '
                     'parameters.\n\n')
            wf.write(f'{"Tip":<60s} {"Projected_Date":>16s} '
                     f'{"Years_Beyond":>16s}\n')
            wf.write('-' * 96 + '\n')
            for name, proj_date, excess in sorted(future_tips, key=lambda x: -x[2]):
                wf.write(f'{name:<60s} {proj_date:>16.4f} {excess:>16.4f}\n')
            wf.write(f'\nTotal: {len(future_tips)} / {len(tips)} tips '
                     f'projected into the future.\n')
        else:
            wf.write(f'No tips are projected beyond the youngest sample '
                     f'date ({youngest_date:.4f}).\n')


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Determine reroot strategy
    if args.reroot is None or len(args.reroot) == 0:
        reroot = 'least-squares'
    elif (len(args.reroot) == 1
          and args.reroot[0] in ('best', 'least-squares', 'min_dev', 'oldest')):
        reroot = args.reroot[0]
    else:
        reroot = args.reroot

    # Determine clock_std for relaxed clock
    clock_std = args.clock_std
    if args.clock_model == 'relaxed' and clock_std is None:
        clock_std = 0.2

    # ---- Run iterative timescale (from tree_time_scale.py) ----
    time_tree, all_outliers_df, final_fasta, final_metadata = iterative_timescale(
        ftree=args.tree,
        falignment=args.aln,
        fdates=args.dates,
        reroot=reroot,
        clock_rate=args.clock_rate,
        clock_std=clock_std if args.clock_model == 'relaxed' else None,
        clock_filter=args.clock_filter,
        clock_filter_method=args.clock_filter_method,
        remove_future_tips=None if args.remove_future_tips < 0 else args.remove_future_tips,
        remove_root=args.remove_root,
        coalescent_tc=args.coalescent,
        sample_id_field=args.sample_id_field,
        collection_date_field=args.date_field,
        rng_seed=args.seed,
        negative_tolerance=args.negative_tolerance,
        max_iterations=args.max_iterations,
    )

    prefix = args.prefix
    _fname = lambda name: fname(args.outdir, prefix, name)

    # ---- Write outputs ----

    # Temporal tree (branch lengths in years)
    Phylo.write(time_tree.tree, _fname('temporal_tree.nwk'),
                format='newick', format_branch_length='%1.8f')

    # Genetic distance tree — reconstruct by multiplying branch lengths by rate
    rate = time_tree.clock_model.get('slope', None)
    if rate and rate > 0:
        gd_tree = copy.deepcopy(time_tree.tree)
        for clade in gd_tree.find_clades():
            if clade.branch_length is not None:
                clade.branch_length = clade.branch_length * rate
        Phylo.write(gd_tree, _fname('genetic_distance_tree.nwk'),
                    format='newick', format_branch_length='%1.8f')
    else:
        Phylo.write(time_tree.tree, _fname('genetic_distance_tree.nwk'),
                    format='newick', format_branch_length='%1.8f')

    # Final alignment
    shutil.copy2(final_fasta, _fname('final.alignment.fasta'))

    # Final metadata (ensure TSV)
    if final_metadata.endswith('.tsv'):
        shutil.copy2(final_metadata, _fname('final.metadata.tsv'))
    else:
        df = pd.read_csv(final_metadata)
        df.to_csv(_fname('final.metadata.tsv'), sep='\t', index=False)

    # Dates TSV
    tips = time_tree.tree.get_terminals()
    dates_records = [
        {'name': tip.name, 'date': tip.raw_date_constraint}
        for tip in tips
        if hasattr(tip, 'raw_date_constraint')
        and tip.raw_date_constraint is not None
    ]
    pd.DataFrame(dates_records).to_csv(
        _fname('final.dates.tsv'), sep='\t', index=False
    )

    # Clock model stats YAML
    write_clock_stats(time_tree, _fname('clock_model_stats.yml'))

    # All outliers TSV
    all_outliers_df.to_csv(_fname('all_outliers.tsv'), sep='\t', index=False)

    # Root-to-tip plot (from tree_time_scale.py)
    fig, ax = plot_root_to_tip(time_tree, all_outliers_df)
    fig.savefig(_fname('root_to_tip.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Temporal tree plot (from tree_time_scale.py)
    fig, ax = plot_temporal_tree(time_tree)
    fig.savefig(_fname('timetree.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Future tips warning (only if not removing them)
    write_future_tips_warning(time_tree, all_outliers_df,
                               _fname('future_tips_warning.txt'),
                               args.remove_future_tips)

    print(f"\nDone. Outputs written to: {args.outdir}/")
    print(f"  Retained tips: {len(tips)}")
    print(f"  Total outliers removed: {len(all_outliers_df)}")


if __name__ == '__main__':
    main()
