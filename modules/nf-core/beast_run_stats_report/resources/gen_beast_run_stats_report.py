#!/usr/bin/env python3
"""Standalone BEAST run stats report generator.

This script builds a Jupyter notebook that summarises the *run statistics* of
the BEAST models used in an ebola-flow pipeline run, then executes it and
exports it to HTML.

The report contains a single section:

1. **BEAST Run Times** - one subsection per BEAST model, detailing the start,
   end and duration of each BEAST run. Run times are extracted from the BEAST
   ``.out`` screen-log files.
"""
import argparse

import nbformat as nbf
from nbconvert import HTMLExporter
from beast_pype.nb_utils import execute_notebook, make_kernelspec


# --------------------------------------------------------------------------- #
# BEAST run times
# --------------------------------------------------------------------------- #
def _beast_runtimes_cells(beast_out_paths):
    """Build notebook cells reporting per-model BEAST run start/end/duration.

    ``beast_out_paths`` is a mapping of ``{model_name: [out_file_path, ...]}``.

    Run duration is extracted from the ``Total calculation time:`` line of each
    BEAST ``.out`` screen-log file. The run *end* time is taken from the
    ``.out`` file modification time, and the run *start* time is derived as
    ``end - duration``.
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
def gen_beast_run_stats_report(output_report_path,
                               beast_out_paths,
                               kernel_name='beast_pype',
                               as_version=4):
    """Generate a notebook report summarising BEAST run statistics.

    Parameters
    ----------
    output_report_path : str
        Path (including filename) where the report notebook will be saved.
    beast_out_paths : dict of {str: list of str}
        Mapping of BEAST model name to a list of its BEAST ``.out`` file paths.
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
            "# BEAST Run Stats Report\n\n"
            "This report summarises the run statistics of the BEAST models "
            "used in this analysis."),
    ]

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
        description="Generate a BEAST run stats report notebook, then execute "
                    "it and export to HTML.")
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
        '-k', '--kernel_name', default='beast_pype',
        help='Name of the Jupyter kernel to use when executing the notebook. '
             'Default: "beast_pype".')
    args = parser.parse_args()

    if len(args.beast_model) != len(args.beast_out):
        parser.error(
            'Each --beast_model must be paired with a --beast_out (same '
            'number, same order).')

    beast_out_paths = _parse_beast_out(args.beast_model, args.beast_out)

    notebook_path = gen_beast_run_stats_report(
        output_report_path=args.output,
        beast_out_paths=beast_out_paths,
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
