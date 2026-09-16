#!/usr/bin/env python3
"""Standalone sample-metadata report generator.

This script inlines ``beast_pype.report_gen.gen_metadata_report`` so the
ebola-flow pipeline can build the sample-metadata report directly from the
underlying source, rather than dispatching through the ``beast_pype
metadata-report`` CLI command. After the notebook is built it is executed and
exported to HTML, mirroring the behaviour of the ``metadata-report`` CLI
command (see ``beast_pype/cli.py``).

The notebook cells generated here import helper functions
(``read_metadata``, ``describe_collection_dates``,
``plot_collection_date_histogram`` etc.) from the installed ``beast_pype``
package at notebook-execution time.

In addition to the collection-date summaries produced by beast_pype, this
version appends a section reporting whether any sequences were submitted on a
Friday evening, Saturday or Sunday in the DRC time zone (UTC+2).
"""
import argparse
import os

import nbformat as nbf
from nbconvert import HTMLExporter
from beast_pype.nb_utils import execute_notebook, make_kernelspec


def _weekend_submission_cells(metadata_paths, xml_set_comparisons):
    """Build notebook cells reporting weekend/Friday-evening submissions.

    Submissions are assessed in the DRC time zone (UTC+2). A submission counts
    as "out of hours" if it was made on a Saturday, a Sunday, or on a Friday at
    or after 18:00 (local UTC+2 time).

    The submission timestamp is read from the ``submittedAtTimestamp`` field
    (Unix epoch seconds) when present. If that field is absent the section
    reports that submission-time information is unavailable.
    """
    if xml_set_comparisons:
        paths_repr = repr(metadata_paths)
        load_line = (
            f"metadata_paths = {paths_repr}\n"
            "submission_frames = []\n"
            "for _xml_set, _path in metadata_paths.items():\n"
            "    _sep = '\\t' if str(_path).endswith(('.tsv', '.txt')) else ','\n"
            "    _df = pd.read_csv(_path, sep=_sep)\n"
            "    _df['xml set'] = _xml_set\n"
            "    submission_frames.append(_df)\n"
            "submission_df = pd.concat(submission_frames, ignore_index=True)"
        )
    else:
        paths_repr = repr(metadata_paths)
        load_line = (
            f"metadata_path = {paths_repr}\n"
            "_sep = '\\t' if str(metadata_path).endswith(('.tsv', '.txt')) else ','\n"
            "submission_df = pd.read_csv(metadata_path, sep=_sep)"
        )

    analysis_code = (
        "import pandas as pd\n"
        "from IPython.display import display, Markdown\n"
        "\n"
        f"{load_line}\n"
        "\n"
        "# DRC time zone is UTC+2 (no daylight saving). 'Etc/GMT-2' == UTC+2.\n"
        "DRC_TZ = 'Etc/GMT-2'\n"
        "timestamp_field = 'submittedAtTimestamp'\n"
        "\n"
        "if timestamp_field not in submission_df.columns:\n"
        "    display(Markdown(\n"
        "        f'No `{timestamp_field}` field found in the metadata, so '\n"
        "        'submission times could not be assessed.'))\n"
        "else:\n"
        "    submitted_utc = pd.to_datetime(\n"
        "        submission_df[timestamp_field], unit='s', utc=True)\n"
        "    submitted_local = submitted_utc.dt.tz_convert(DRC_TZ)\n"
        "    dayofweek = submitted_local.dt.dayofweek  # Mon=0 .. Sun=6\n"
        "    hour = submitted_local.dt.hour\n"
        "\n"
        "    friday_evening = (dayofweek == 4) & (hour >= 18)\n"
        "    saturday = dayofweek == 5\n"
        "    sunday = dayofweek == 6\n"
        "    out_of_hours = friday_evening | saturday | sunday\n"
        "\n"
        "    report = submission_df.copy()\n"
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
        "        if 'xml set' in report.columns:\n"
        "            show_cols = ['xml set'] + show_cols\n"
        "        display(report.loc[out_of_hours, show_cols]\n"
        "                .sort_values('Submitted (UTC+2)')\n"
        "                .reset_index(drop=True))\n"
        "    else:\n"
        "        display(Markdown('No out-of-hours submissions were found.'))"
    )

    return [
        nbf.v4.new_markdown_cell(
            "## Out-of-hours Submissions (DRC time zone, UTC+2)\n\n"
            "This section reports whether any sequences were submitted on a "
            "**Friday evening (>= 18:00), Saturday or Sunday** in the "
            "Democratic Republic of the Congo time zone (UTC+2)."),
        nbf.v4.new_code_cell(analysis_code),
    ]


def gen_metadata_report(output_report_path,
                        metadata_paths,
                        collection_date_field='collection_date',
                        xml_set_comparisons=False,
                        xml_set_label='xml set',
                        kernel_name='beast_pype',
                        as_version=4):
    """Generate a notebook report summarising sample metadata.

    Copied from ``beast_pype.report_gen.gen_metadata_report`` so the
    report-generation logic lives with the pipeline module, then extended with
    an out-of-hours (weekend / Friday evening) submission section.

    For simple workflows (``xml_set_comparisons=False``) the report summarises a
    single metadata file via :func:`pandas.DataFrame.describe` and a histogram.

    For comparative workflows (``xml_set_comparisons=True``) the report starts
    with a combined summary and histogram of all sequences, then adds a section
    comparing each xml set via :func:`pandas.DataFrame.describe` and a stacked
    histogram.

    Parameters
    ----------
    output_report_path : str
        Path (including filename) where the report notebook will be saved.
    metadata_paths : str or dict of {str: str}
        If ``xml_set_comparisons`` is False, a path to a single metadata file
        (``.csv`` or ``.tsv``). If ``xml_set_comparisons`` is True, a mapping of
        xml set name to its metadata file path.
    collection_date_field : str, default 'collection_date'
        Name of the field in the metadata containing collection dates
        (formatted YYYY-MM-DD).
    xml_set_comparisons : bool, default False
        Whether the report should compare xml sets.
    xml_set_label : str, default 'xml set'
        Label used for the xml set grouping variable in comparative reports.
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
            "# Sample Collection Date Report\n\n" +
            "This report summarises the sample collection date metadata used in "
            "this analysis."
        ),
        nbf.v4.new_code_cell(
            "import warnings\n" +
            "warnings.filterwarnings('ignore')\n" +
            "import pandas as pd\n" +
            "from beast_pype.outputs import (read_metadata, describe_collection_dates,\n" +
            "    plot_collection_date_histogram)"
        ),
    ]

    if xml_set_comparisons:
        if not isinstance(metadata_paths, dict):
            raise ValueError(
                "metadata_paths must be a dict of {xml_set: path} when "
                "xml_set_comparisons is True.")
        notebook['cells'].append(
            nbf.v4.new_code_cell(
                "from beast_pype.outputs import (\n" +
                "    describe_collection_dates_by_xml_set,\n" +
                "    plot_stacked_collection_date_histogram)"
            ))
        notebook['cells'].append(
            nbf.v4.new_code_cell(
                f"collection_date_field = '{collection_date_field}'\n" +
                f"xml_set_label = '{xml_set_label}'\n" +
                f"metadata_paths = {metadata_paths!r}\n" +
                "metadata_dict = {xml_set: read_metadata(path, collection_date_field)\n" +
                "    for xml_set, path in metadata_paths.items()}\n" +
                "combined_metadata = pd.concat(metadata_dict.values(), ignore_index=True)"
            ))
        notebook['cells'] += [
            nbf.v4.new_markdown_cell(
                "## All Sequences\n\n" +
                "### Summary statistics of collection dates"),
            nbf.v4.new_code_cell(
                "describe_collection_dates(combined_metadata, collection_date_field)"),
            nbf.v4.new_markdown_cell("### Histogram of collection dates"),
            nbf.v4.new_code_cell(
                "fig, ax = plot_collection_date_histogram(combined_metadata, collection_date_field)\n"),
            nbf.v4.new_markdown_cell(
                "## Comparison of XML Sets\n\n" +
                "### Summary statistics of collection dates per xml set"),
            nbf.v4.new_code_cell(
                "describe_collection_dates_by_xml_set(metadata_dict, collection_date_field, xml_set_label=xml_set_label)"),
            nbf.v4.new_markdown_cell("### Stacked histogram of collection dates per xml set"),
            nbf.v4.new_code_cell(
                "fig, ax = plot_stacked_collection_date_histogram(metadata_dict, collection_date_field, xml_set_label=xml_set_label)\n"),
        ]
    else:
        if isinstance(metadata_paths, dict):
            raise ValueError(
                "metadata_paths must be a single path (str) when "
                "xml_set_comparisons is False.")
        notebook['cells'].append(
            nbf.v4.new_code_cell(
                f"collection_date_field = '{collection_date_field}'\n" +
                f"metadata_path = '{metadata_paths}'\n" +
                "metadata = read_metadata(metadata_path, collection_date_field)"
            ))
        notebook['cells'] += [
            nbf.v4.new_markdown_cell(
                "## Summary statistics of collection dates"),
            nbf.v4.new_code_cell(
                "describe_collection_dates(metadata, collection_date_field)"),
            nbf.v4.new_markdown_cell("## Histogram of collection dates"),
            nbf.v4.new_code_cell(
                "fig, ax = plot_collection_date_histogram(metadata, collection_date_field)\n" +
                "fig"),
        ]

    # --- Out-of-hours (weekend / Friday evening) submission section ---
    notebook['cells'] += _weekend_submission_cells(
        metadata_paths, xml_set_comparisons)

    with open(output_report_path, 'w') as f:
        nbf.write(notebook, f)

    return output_report_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate a report notebook summarising sample metadata, "
                    "then execute it and export to HTML.")
    parser.add_argument(
        '-m', '--metadata', action='append', required=True,
        help='Path to a metadata file (.csv or .tsv). Can be specified multiple '
             'times. When more than one is provided, comparative (xml-set) '
             'visualisations are used.')
    parser.add_argument(
        '-n', '--xml_set_name', action='append', default=[],
        help='Name for an xml set, paired with each --metadata in order. Only '
             'used (and required) when more than one --metadata is provided.')
    parser.add_argument(
        '-o', '--output', required=True,
        help='Path to save the output report notebook (.ipynb).')
    parser.add_argument(
        '-c', '--collection_date_field', default='collection_date',
        help='Name of field in metadata containing collection dates '
             '(YYYY-MM-DD). Default: "collection_date".')
    parser.add_argument(
        '--xml_set_label', default='xml set',
        help='Label for the xml set grouping variable in comparative reports. '
             'Default: "xml set".')
    parser.add_argument(
        '-k', '--kernel_name', default='beast_pype',
        help='Name of the Jupyter kernel to use when executing the notebook. '
             'Default: "beast_pype".')
    args = parser.parse_args()

    xml_set_comparisons = len(args.metadata) > 1
    if xml_set_comparisons:
        if len(args.xml_set_name) != len(args.metadata):
            parser.error(
                'When more than one --metadata is provided, a matching '
                '--xml_set_name must be given for each (same number, same order).')
        metadata_paths = dict(zip(args.xml_set_name, args.metadata))
    else:
        metadata_paths = args.metadata[0]

    notebook_path = gen_metadata_report(
        output_report_path=args.output,
        metadata_paths=metadata_paths,
        collection_date_field=args.collection_date_field,
        xml_set_comparisons=xml_set_comparisons,
        xml_set_label=args.xml_set_label,
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
