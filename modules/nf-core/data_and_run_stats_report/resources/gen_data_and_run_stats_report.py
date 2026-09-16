#!/usr/bin/env python3
"""Standalone data-and-run-stats report generator.

This script builds a Jupyter notebook that summarises both the *data* and the
*run statistics* of an ebola-flow pipeline run, then executes it and exports it
to HTML (mirroring the behaviour of the old ``metadata_report`` module).

The report is organised into the following sections:

1. **Pathoplexus Run Date and Time** - the date and time the Pathoplexus data
   was fetched for this run.
2. **Metadata Report - Raw Pathoplexus Data** - the collection-date summaries
   (plus out-of-hours submission check) run on the *raw* Pathoplexus metadata.
3. **Metadata Report - After Outlier Removal** - the same collection-date
   summaries run on the metadata *after* outliers are removed in the
   ``modules/nf-core/iterative_treetime`` step.
4. **BEAST Run Times** - one subsection per BEAST model, detailing the start,
   end and duration of each BEAST run. Run times are extracted from the BEAST
   ``.out`` screen-log files (see the ``runtime.py`` helper logic that is
   inlined into the notebook).

The notebook cells import helper functions (``read_metadata``,
``describe_collection_dates``, ``plot_collection_date_histogram`` etc.) from the
installed ``beast_pype`` package at notebook-execution time.
"""
import argparse

import nbformat as nbf
from nbconvert import HTMLExporter
from beast_pype.nb_utils import execute_notebook, make_kernelspec


# --------------------------------------------------------------------------- #
# Section 1: Pathoplexus run date and time
# --------------------------------------------------------------------------- #
def _pathoplexus_datetime_cells(timestamp_path):
    """Build notebook cells reporting the Pathoplexus fetch date and time.

    ``timestamp_path`` is a file whose contents are an ISO-8601 timestamp (or
    a Unix epoch in seconds). If the file is missing or cannot be parsed the
    modification time of the file is used as a fallback.
    """
    code = (
        "import os\n"
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "\n"
        "timestamp_path = " + repr(timestamp_path) + "\n"
        "fetched_at = None\n"
        "if timestamp_path and os.path.exists(timestamp_path):\n"
        "    with open(timestamp_path) as _f:\n"
        "        _raw = _f.read().strip()\n"
        "    if _raw:\n"
        "        try:\n"
        "            fetched_at = pd.to_datetime(float(_raw), unit='s')\n"
        "        except ValueError:\n"
        "            try:\n"
        "                fetched_at = pd.to_datetime(_raw)\n"
        "            except (ValueError, TypeError):\n"
        "                fetched_at = None\n"
        "    if fetched_at is None:\n"
        "        fetched_at = pd.to_datetime(os.path.getmtime(timestamp_path), unit='s')\n"
        "\n"
        "if fetched_at is None:\n"
        "    display(Markdown('The Pathoplexus fetch date and time is unavailable.'))\n"
        "else:\n"
        "    display(Markdown(\n"
        "        'The Pathoplexus data for this run was fetched on '\n"
        "        f'**{fetched_at.strftime(\"%Y-%m-%d %H:%M:%S\")}**.'))"
    )
    return [
        nbf.v4.new_markdown_cell(
            "## Pathoplexus Run Date and Time\n\n"
            "This section reports the date and time the Pathoplexus data was "
            "fetched for this pipeline run."),
        nbf.v4.new_code_cell(code),
    ]


# --------------------------------------------------------------------------- #
# Section 2 & 3: Metadata report (run on raw and outlier-removed data)
# --------------------------------------------------------------------------- #
def _metadata_report_cells(metadata_path, section_title, collection_date_field,
                           var_suffix):
    """Build the collection-date summary + out-of-hours cells for one dataset.

    ``var_suffix`` keeps notebook variable names unique between the two
    sections that share the same logic (raw vs. outlier-removed data).
    """
    md_var = "metadata_" + var_suffix
    path_var = "metadata_path_" + var_suffix
    sub_var = "submission_df_" + var_suffix

    load_code = (
        "import warnings\n"
        "warnings.filterwarnings('ignore')\n"
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "from beast_pype.outputs import (read_metadata, describe_collection_dates,\n"
        "    plot_collection_date_histogram)\n"
        "\n"
        "collection_date_field = '" + collection_date_field + "'\n"
        + path_var + " = " + repr(metadata_path) + "\n"
        + md_var + " = read_metadata(" + path_var + ", collection_date_field)"
    )

    # Out-of-hours (weekend / Friday evening) submission check, in DRC UTC+2.
    weekend_code = (
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "\n"
        "_sep = '\\t' if str(" + path_var + ").endswith(('.tsv', '.txt')) else ','\n"
        + sub_var + " = pd.read_csv(" + path_var + ", sep=_sep)\n"
        "\n"
        "# DRC time zone is UTC+2 (no daylight saving). 'Etc/GMT-2' == UTC+2.\n"
        "DRC_TZ = 'Etc/GMT-2'\n"
        "timestamp_field = 'submittedAtTimestamp'\n"
        "\n"
        "if timestamp_field not in " + sub_var + ".columns:\n"
        "    display(Markdown(\n"
        "        f'No `{timestamp_field}` field found in the metadata, so '\n"
        "        'submission times could not be assessed.'))\n"
        "else:\n"
        "    submitted_utc = pd.to_datetime(\n"
        "        " + sub_var + "[timestamp_field], unit='s', utc=True)\n"
        "    submitted_local = submitted_utc.dt.tz_convert(DRC_TZ)\n"
        "    dayofweek = submitted_local.dt.dayofweek  # Mon=0 .. Sun=6\n"
        "    hour = submitted_local.dt.hour\n"
        "\n"
        "    friday_evening = (dayofweek == 4) & (hour >= 18)\n"
        "    saturday = dayofweek == 5\n"
        "    sunday = dayofweek == 6\n"
        "    out_of_hours = friday_evening | saturday | sunday\n"
        "\n"
        "    report = " + sub_var + ".copy()\n"
        "    report['Submitted (UTC+2)'] = submitted_local.dt.strftime('%Y-%m-%d %H:%M:%S')\n"
        "    report['Submission day'] = submitted_local.dt.day_name()\n"
        "\n"
        "    n_total = len(report)\n"
        "    n_flagged = int(out_of_hours.sum())\n"
        "    display(Markdown(\n"
        "        f'**{n_flagged}** of **{n_total}** sequences were submitted on a '\n"
        "        'Friday evening (>= 18:00), Saturday or Sunday in the DRC time '\n"
        "        'zone (UTC+2):\\n\\n'\n"
        "        f'- Friday evening (>= 18:00): {int(friday_evening.sum())}\\n'\n"
        "        f'- Saturday: {int(saturday.sum())}\\n'\n"
        "        f'- Sunday: {int(sunday.sum())}'))\n"
        "\n"
        "    if n_flagged:\n"
        "        id_cols = [c for c in ('accessionVersion', 'accession', 'strain')\n"
        "                   if c in report.columns]\n"
        "        show_cols = id_cols + ['Submitted (UTC+2)', 'Submission day']\n"
        "        display(report.loc[out_of_hours, show_cols]\n"
        "                .sort_values('Submitted (UTC+2)')\n"
        "                .reset_index(drop=True))\n"
        "    else:\n"
        "        display(Markdown('No out-of-hours submissions were found.'))"
    )

    cells = [
        nbf.v4.new_markdown_cell("## " + section_title),
        nbf.v4.new_code_cell(load_code),
        nbf.v4.new_markdown_cell("### Summary statistics of collection dates"),
        nbf.v4.new_code_cell(
            "describe_collection_dates(" + md_var + ", collection_date_field)"),
        nbf.v4.new_markdown_cell("### Histogram of collection dates"),
        nbf.v4.new_code_cell(
            "fig, ax = plot_collection_date_histogram(" + md_var + ", collection_date_field)"),
        nbf.v4.new_markdown_cell(
            "### Out-of-hours submissions (DRC time zone, UTC+2)\n\n"
            "Whether any sequences were submitted on a **Friday evening "
            "(>= 18:00), Saturday or Sunday** in the Democratic Republic of "
            "the Congo time zone (UTC+2)."),
        nbf.v4.new_code_cell(weekend_code),
    ]
    cells += _geo_location_cells(path_var, var_suffix)
    return cells


# --------------------------------------------------------------------------- #
# Geographic location breakdowns (bar charts + raw count tables)
# --------------------------------------------------------------------------- #
def _geo_location_cells(path_var, var_suffix):
    """Build geographic-location bar-chart + count-table cells for one dataset.

    Produces, for each of ``geoLocCountry``, ``geoLocAdmin1``, ``geoLocAdmin2``
    and ``geoLocCity``, a bar chart and a raw count table (including a count of
    missing/unspecified entries). ``geoLocAdmin2`` values are normalised to
    title case before counting so that block-caps (e.g. ``BUTEMBO``) and
    first-letter-caps (e.g. ``Butembo``) variants are collapsed into a single
    location.

    ``path_var`` is the name of the notebook variable already holding the path
    to this dataset's metadata file (set in ``_metadata_report_cells``). This
    reuses that loaded path rather than re-reading it from disk.
    """
    geo_var = "geo_df_" + var_suffix

    code = (
        "import matplotlib.pyplot as plt\n"
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "\n"
        "_sep = '\\t' if str(" + path_var + ").endswith(('.tsv', '.txt')) else ','\n"
        + geo_var + " = pd.read_csv(" + path_var + ", sep=_sep, dtype=str)\n"
        "\n"
        "# Fields to summarise, in order of geographic granularity. The bool\n"
        "# flag marks fields whose values should be normalised to title case\n"
        "# before counting (to collapse block-caps vs first-letter-caps).\n"
        "geo_fields = [\n"
        "    ('geoLocCountry', False),\n"
        "    ('geoLocAdmin1', False),\n"
        "    ('geoLocAdmin2', True),\n"
        "    ('geoLocCity', False),\n"
        "]\n"
        "\n"
        "for field, normalise_case in geo_fields:\n"
        "    display(Markdown(f'### {field}'))\n"
        "    if field not in " + geo_var + ".columns:\n"
        "        display(Markdown(f'No `{field}` field found in the metadata.'))\n"
        "        continue\n"
        "    values = " + geo_var + "[field].copy()\n"
        "    # Treat blank/whitespace-only strings as missing.\n"
        "    values = values.where(~values.fillna('').str.strip().eq(''), other=pd.NA)\n"
        "    # Treat sentinel 'no data' strings (case-insensitive) as missing too.\n"
        "    sentinels = {'not provided', 'unknown', 'unspecified', 'na', 'n/a',\n"
        "                 'none', 'null', 'not applicable', 'missing'}\n"
        "    _norm = values.fillna('').str.strip().str.lower()\n"
        "    values = values.where(~_norm.isin(sentinels), other=pd.NA)\n"
        "    if normalise_case:\n"
        "        values = values.where(values.isna(), values.str.strip().str.title())\n"
        "    counts = values.value_counts(dropna=False)\n"
        "    counts.index = [('(unspecified)' if pd.isna(i) else i) for i in counts.index]\n"
        "    count_table = (counts.rename_axis(field).reset_index(name='count'))\n"
        "    display(count_table)\n"
        "    # Bar chart of counts, including the '(unspecified)' bucket as its\n"
        "    # own bar so missing/blank entries are visible in the chart too.\n"
        "    plot_counts = counts.sort_values(ascending=False)\n"
        "    if plot_counts.empty:\n"
        "        display(Markdown('No values to plot for this field.'))\n"
        "        continue\n"
        "    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(plot_counts)), 4))\n"
        "    ax.bar(plot_counts.index.astype(str), plot_counts.values)\n"
        "    ax.set_xlabel(field)\n"
        "    ax.set_ylabel('Number of sequences')\n"
        "    ax.set_title(f'Sequences by {field}')\n"
        "    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')\n"
        "    fig.tight_layout()\n"
        "    plt.show()"
    )

    return [
        nbf.v4.new_markdown_cell(
            "### Geographic location breakdown\n\n"
            "Bar charts and raw count tables of sequences by "
            "`geoLocCountry`, `geoLocAdmin1`, `geoLocAdmin2` and "
            "`geoLocCity`. `geoLocAdmin2` values are normalised to title "
            "case so that block-caps and first-letter-caps variants of the "
            "same location are counted together. Missing or blank entries, "
            "along with sentinel 'no data' strings such as `Not Provided` and "
            "`Unknown`, are reported as `(unspecified)` in both the count "
            "tables and the bar charts."),
        nbf.v4.new_code_cell(code),
    ]


# --------------------------------------------------------------------------- #
# Section 4: BEAST run times
# --------------------------------------------------------------------------- #
def _beast_runtimes_cells(beast_out_paths):
    """Build notebook cells reporting per-model BEAST run start/end/duration.

    ``beast_out_paths`` is a mapping of ``{model_name: [out_file_path, ...]}``.

    Run duration is extracted from the ``Total calculation time:`` line of each
    BEAST ``.out`` screen-log file (see ``runtime.py``). The run *end* time is
    taken from the ``.out`` file modification time, and the run *start* time is
    derived as ``end - duration``.
    """
    code = (
        "import os\n"
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "\n"
        "beast_out_paths = " + repr(beast_out_paths) + "\n"
        "\n"
        "\n"
        "def _extract_calc_seconds(out_path):\n"
        "    \"\"\"Extract 'Total calculation time:' (seconds) from a BEAST .out file.\"\"\"\n"
        "    val_name = 'Total calculation time:'\n"
        "    with open(out_path) as _f:\n"
        "        for line in _f:\n"
        "            if line.strip().startswith(val_name):\n"
        "                cleaned = line.split(':', 1)[1]\n"
        "                number = ''.join(c for c in cleaned if c.isdigit() or c == '.')\n"
        "                if number:\n"
        "                    return float(number)\n"
        "    return None\n"
        "\n"
        "\n"
        "for model, out_files in beast_out_paths.items():\n"
        "    display(Markdown(f'### {model}'))\n"
        "    records = []\n"
        "    for out_file in out_files:\n"
        "        run_name = os.path.basename(out_file).replace('.out', '')\n"
        "        seconds = _extract_calc_seconds(out_file) if os.path.exists(out_file) else None\n"
        "        if seconds is None:\n"
        "            records.append({'BEAST run': run_name, 'Start': 'NA',\n"
        "                            'End': 'NA', 'Duration (D H:M:S)': 'NA'})\n"
        "            continue\n"
        "        duration = pd.to_timedelta(seconds, unit='s')\n"
        "        end = pd.to_datetime(os.path.getmtime(out_file), unit='s')\n"
        "        start = end - duration\n"
        "        records.append({\n"
        "            'BEAST run': run_name,\n"
        "            'Start': start.strftime('%Y-%m-%d %H:%M:%S'),\n"
        "            'End': end.strftime('%Y-%m-%d %H:%M:%S'),\n"
        "            'Duration (D H:M:S)': str(duration),\n"
        "        })\n"
        "    display(pd.DataFrame(records).sort_values('BEAST run').reset_index(drop=True))"
    )
    return [
        nbf.v4.new_markdown_cell(
            "## BEAST Run Times\n\n"
            "Each subsection below corresponds to a BEAST model and details "
            "the start, end and duration of each BEAST run for that model. "
            "Durations are extracted from the `Total calculation time:` line "
            "of the BEAST `.out` screen-log files."),
        nbf.v4.new_code_cell(code),
    ]


# --------------------------------------------------------------------------- #
# Notebook assembly
# --------------------------------------------------------------------------- #
def gen_data_and_run_stats_report(output_report_path,
                                  raw_metadata_path,
                                  cleaned_metadata_path,
                                  filtered_metadata_path,
                                  pathoplexus_timestamp_path,
                                  beast_out_paths,
                                  collection_date_field='collection_date',
                                  kernel_name='beast_pype',
                                  as_version=4):
    """Generate a notebook report summarising data and run statistics.

    Parameters
    ----------
    output_report_path : str
        Path (including filename) where the report notebook will be saved.
    raw_metadata_path : str
        Path to the raw Pathoplexus metadata file (.csv or .tsv).
    cleaned_metadata_path : str
        Path to the metadata file after the ``clean_metadata`` step
        (excluded IDs and invalid dates removed) (.csv or .tsv).
    filtered_metadata_path : str
        Path to the metadata file after outliers are removed in the
        ``iterative_treetime`` step (.csv or .tsv).
    pathoplexus_timestamp_path : str
        Path to a file containing the Pathoplexus fetch timestamp.
    beast_out_paths : dict of {str: list of str}
        Mapping of BEAST model name to a list of its BEAST ``.out`` file paths.
    collection_date_field : str, default 'collection_date'
        Name of the field in the metadata containing collection dates
        (formatted YYYY-MM-DD).
    kernel_name : str, default 'beast_pype'
        Name of the Jupyter (python) kernel to use when executing the notebook.
    as_version : int, default 4
        Jupyter notebook version.

    Returns
    -------
    str
        Path to the generated notebook.
    """
    notebook = nbf.v4.new_notebook()
    notebook['metadata']['kernelspec'] = make_kernelspec(kernel_name)

    notebook['cells'] = [
        nbf.v4.new_markdown_cell(
            "# Data and Run Stats Report\n\n"
            "This report summarises both the data used in this analysis "
            "(sample collection dates and submission times) and the run "
            "statistics of the BEAST models used."),
    ]

    # Section 1: Pathoplexus run date and time
    notebook['cells'] += _pathoplexus_datetime_cells(pathoplexus_timestamp_path)

    # Section 2: metadata report on the raw Pathoplexus pull
    notebook['cells'] += _metadata_report_cells(
        raw_metadata_path,
        "Metadata Report - Raw Pathoplexus Pull",
        collection_date_field,
        var_suffix='raw')

    # Section 3: metadata report after the clean_metadata step
    notebook['cells'] += _metadata_report_cells(
        cleaned_metadata_path,
        "Metadata Report - After clean_metadata Step",
        collection_date_field,
        var_suffix='cleaned')

    # Section 4: metadata report after outlier removal (iterative_treetime)
    notebook['cells'] += _metadata_report_cells(
        filtered_metadata_path,
        "Metadata Report - After Outlier Removal (iterative_treetime)",
        collection_date_field,
        var_suffix='filtered')

    # Section 5: BEAST run times (per model)
    notebook['cells'] += _beast_runtimes_cells(beast_out_paths)

    with open(output_report_path, 'w') as f:
        nbf.write(notebook, f)

    return output_report_path


def _parse_beast_out(model_args, out_args):
    """Pair each --beast_model with its (comma-separated) --beast_out paths."""
    beast_out_paths = {}
    for model, outs in zip(model_args, out_args):
        paths = [p for p in outs.split(',') if p]
        beast_out_paths[model] = paths
    return beast_out_paths


def main():
    parser = argparse.ArgumentParser(
        description="Generate a data-and-run-stats report notebook, then "
                    "execute it and export to HTML.")
    parser.add_argument(
        '-r', '--raw_metadata', required=True,
        help='Path to the raw Pathoplexus metadata file (.csv or .tsv).')
    parser.add_argument(
        '-m', '--cleaned_metadata', required=True,
        help='Path to the metadata file after the clean_metadata step '
             '(.csv or .tsv).')
    parser.add_argument(
        '-f', '--filtered_metadata', required=True,
        help='Path to the outlier-removed metadata file (.csv or .tsv).')
    parser.add_argument(
        '-t', '--pathoplexus_timestamp', required=True,
        help='Path to a file containing the Pathoplexus fetch timestamp.')
    parser.add_argument(
        '--beast_model', action='append', default=[],
        help='Name of a BEAST model. Pair with --beast_out in order.')
    parser.add_argument(
        '--beast_out', action='append', default=[],
        help='Comma-separated list of BEAST .out files for the paired '
             '--beast_model.')
    parser.add_argument(
        '-o', '--output', required=True,
        help='Path to save the output report notebook (.ipynb).')
    parser.add_argument(
        '-c', '--collection_date_field', default='collection_date',
        help='Name of field in metadata containing collection dates '
             '(YYYY-MM-DD). Default: "collection_date".')
    parser.add_argument(
        '-k', '--kernel_name', default='beast_pype',
        help='Name of the Jupyter kernel to use when executing the notebook. '
             'Default: "beast_pype".')
    args = parser.parse_args()

    if len(args.beast_model) != len(args.beast_out):
        parser.error(
            'Each --beast_model must be paired with a --beast_out (same '
            'number, same order).')

    beast_out_paths = _parse_beast_out(args.beast_model, args.beast_out)

    notebook_path = gen_data_and_run_stats_report(
        output_report_path=args.output,
        raw_metadata_path=args.raw_metadata,
        cleaned_metadata_path=args.cleaned_metadata,
        filtered_metadata_path=args.filtered_metadata,
        pathoplexus_timestamp_path=args.pathoplexus_timestamp,
        beast_out_paths=beast_out_paths,
        collection_date_field=args.collection_date_field,
        kernel_name=args.kernel_name,
    )

    execute_notebook(
        input_path=notebook_path,
        output_path=notebook_path,
        kernel_name=args.kernel_name,
        progress_bar=True,
    )

    executed_nb = nbf.read(notebook_path, as_version=4)
    html_exporter = HTMLExporter(exclude_input=True)
    html_body, _ = html_exporter.from_notebook_node(executed_nb)
    html_path = notebook_path.replace('.ipynb', '.html')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_body)

    print(f"Notebook: {notebook_path}")
    print(f"Notebook HTML: {html_path}")


if __name__ == '__main__':
    main()
