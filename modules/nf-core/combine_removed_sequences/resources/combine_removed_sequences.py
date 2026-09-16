#!/usr/bin/env python3
"""Combine excluded-samples and clock-outlier records into one report.

Merges:
  * clean_metadata's ``*_excluded_samples.csv`` (exclude-list hits and
    invalid-date removals), columns:
    accessionVersion, accession, reason_for_removal, source_for_reason
  * iterative_treetime's ``all_outliers.tsv`` (TreeTime clock outliers),
    columns: iteration, name, given_date, apparent_date, z-score, diagnosis

into a single tidy CSV describing which sequences were removed and why.
"""

import argparse
import pandas as pd


OUT_COLUMNS = [
    "accessionVersion",
    "accession",
    "removal_step",
    "reason_for_removal",
    "source_for_reason",
    "details",
]

# Human-readable descriptions for TreeTime outlier diagnoses.
DIAGNOSIS_REASON = {
    "date_too_late": "Sampling date inconsistent with molecular clock (too late)",
    "date_too_early": "Sampling date inconsistent with molecular clock (too early)",
    "excess_mutations": "Excess mutations relative to molecular clock",
}


def _strip_version(accession_version):
    """Return the accession without a trailing .<version> suffix."""
    if not isinstance(accession_version, str):
        return ""
    return accession_version.rsplit(".", 1)[0]


def load_excluded(path):
    df = pd.read_csv(path, dtype=str).fillna("")
    if df.empty:
        return pd.DataFrame(columns=OUT_COLUMNS)
    out = pd.DataFrame({
        "accessionVersion": df.get("accessionVersion", ""),
        "accession": df.get("accession", ""),
        "removal_step": "clean_metadata",
        "reason_for_removal": df.get("reason_for_removal", ""),
        "source_for_reason": df.get("source_for_reason", ""),
        "details": "",
    })
    return out[OUT_COLUMNS]


def load_outliers(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        return pd.DataFrame(columns=OUT_COLUMNS)
    diagnosis = df.get("diagnosis", "")
    reason = diagnosis.map(
        lambda d: DIAGNOSIS_REASON.get(d, d if d else "Clock outlier")
    )
    details = df.apply(
        lambda r: "; ".join(
            part for part in [
                f"iteration={r['iteration']}" if r.get("iteration", "") else "",
                f"given_date={r['given_date']}" if r.get("given_date", "") else "",
                f"apparent_date={r['apparent_date']}" if r.get("apparent_date", "") else "",
                f"z-score={r['z-score']}" if r.get("z-score", "") else "",
            ] if part
        ),
        axis=1,
    ) if len(df) else ""
    accession_version = df.get("name", "")
    out = pd.DataFrame({
        "accessionVersion": accession_version,
        "accession": accession_version.map(_strip_version),
        "removal_step": "iterative_treetime",
        "reason_for_removal": reason,
        "source_for_reason": "TreeTime clock filter",
        "details": details,
    })
    return out[OUT_COLUMNS]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--excluded", required=True,
                        help="clean_metadata excluded_samples.csv")
    parser.add_argument("--outliers", required=True,
                        help="iterative_treetime all_outliers.tsv")
    parser.add_argument("-o", "--output", required=True,
                        help="output combined CSV path")
    args = parser.parse_args()

    frames = [load_excluded(args.excluded), load_outliers(args.outliers)]
    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["accessionVersion"].astype(str).str.strip() != ""]
    combined = combined.sort_values(
        ["removal_step", "accessionVersion"]
    ).reset_index(drop=True)
    combined.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
