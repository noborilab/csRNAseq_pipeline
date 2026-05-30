#!/usr/bin/env python3
"""
Generate a csRNA-seq pipeline sample sheet from a directory of FASTQ files.

Usage
-----
    python scripts/make_sample_sheet.py /fastq > config/samples.tsv
    python scripts/make_sample_sheet.py /fastq -o config/samples.tsv

The script walks <fastq_root> and discovers .fastq.gz files.  It groups them
using the *top-level* directory name directly under <fastq_root>, which must
follow the naming convention:

    {condition}_{cs|in}{replicate}[_redo]

Examples of valid top-level names:
    col0_12d_seedlings_cs1          → (col0_12d_seedlings, csrna, 1)
    col0_12d_seedlings_in2          → (col0_12d_seedlings, input, 2)
    col0_flg22_mock_30m_cs3_redo    → merged into (col0_flg22_mock_30m, csrna, 3)

R1 / R2 are detected from _R1_ / _R2_ in the file name (standard Illumina).
Files that are neither R1 nor R2 are skipped.

Multiple files for the same (sample, type, replicate, strand) are joined with
commas in the output (multi-lane / multi-run / redo merging).

Directories whose names do not match the expected pattern are skipped with a
warning printed to stderr.  Add those rows manually or use --rename to map them.

--rename
--------
    Pass one or more rename rules to handle non-standard top-level directory
    names.  Each rule is of the form:

        DIRNAME=condition,type,replicate

    where type is 'cs' or 'in'.  Example:

        --rename 'col0_3d_seedlings_cs5__input_=col0_3d_seedlings,in,5'

Output columns
--------------
    sample_name  sample_type  replicate  read_r1  read_r2

    read_r2 is empty for single-end libraries (no R2 files found for a sample).
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# regex: {anything}_{cs|in}{digits}[_redo][end-of-string]
_PATTERN = re.compile(r"^(.+?)_(cs|in)(\d+)(_redo)?$", re.IGNORECASE)


def parse_topdir(name: str) -> tuple[str, str, int, bool] | None:
    """Return (condition, type_abbrev, replicate, is_redo) or None."""
    m = _PATTERN.match(name)
    if not m:
        return None
    condition, stype, rep, redo = m.groups()
    return condition, stype.lower(), int(rep), redo is not None


def classify_strand(filename: str) -> str | None:
    """Return 'R1', 'R2', or None."""
    if re.search(r"[_.]R1[_.]", filename) or filename.endswith("_R1.fastq.gz"):
        return "R1"
    if re.search(r"[_.]R2[_.]", filename) or filename.endswith("_R2.fastq.gz"):
        return "R2"
    return None


def collect_files(root: str, rename_map: dict[str, tuple]) -> dict:
    """
    Walk root; return dict keyed by (condition, sample_type, replicate)
    mapping to {"R1": [paths], "R2": [paths]}.
    """
    groups: dict = defaultdict(lambda: {"R1": [], "R2": []})
    skipped: set[str] = set()

    for entry in sorted(os.scandir(root), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        topdir = entry.name

        if topdir in rename_map:
            condition, stype, rep = rename_map[topdir]
        else:
            parsed = parse_topdir(topdir)
            if parsed is None:
                skipped.add(topdir)
                continue
            condition, stype, rep, _ = parsed  # redo flag already handled by key merge

        key = (condition, stype, rep)

        for dirpath, _, filenames in os.walk(entry.path):
            for fname in sorted(filenames):
                if not fname.endswith(".fastq.gz"):
                    continue
                strand = classify_strand(fname)
                if strand is None:
                    continue
                full = os.path.join(dirpath, fname)
                groups[key][strand].append(full)

    if skipped:
        for d in sorted(skipped):
            print(f"WARNING: skipping unrecognised directory: {d}", file=sys.stderr)

    return groups


def build_rows(groups: dict) -> list[dict]:
    type_map = {"cs": "csrna", "in": "input"}
    rows = []
    for (condition, stype, rep), strands in groups.items():
        r1_files = strands["R1"]
        r2_files = strands["R2"]
        if not r1_files:
            print(
                f"WARNING: no R1 files for ({condition}, {stype}, {rep}), skipping",
                file=sys.stderr,
            )
            continue
        rows.append(
            {
                "sample_name": condition,
                "sample_type": type_map.get(stype, stype),
                "replicate": rep,
                "read_r1": ",".join(r1_files),
                "read_r2": ",".join(r2_files),
            }
        )

    # sort: condition alphabetically, then input before csrna, then replicate
    type_order = {"csrna": 0, "input": 1}
    rows.sort(
        key=lambda r: (r["sample_name"], type_order.get(r["sample_type"], 99), r["replicate"])
    )
    return rows


def parse_rename(specs: list[str]) -> dict[str, tuple]:
    rename_map: dict[str, tuple] = {}
    for spec in specs:
        try:
            dirname, rest = spec.split("=", 1)
            condition, stype, rep_str = rest.split(",")
            rename_map[dirname.strip()] = (condition.strip(), stype.strip().lower(), int(rep_str.strip()))
        except ValueError:
            print(f"ERROR: invalid --rename spec: {spec!r}", file=sys.stderr)
            sys.exit(1)
    return rename_map


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("fastq_root", help="Root directory containing per-sample FASTQ subdirectories")
    parser.add_argument("-o", "--output", default="-", help="Output TSV path (default: stdout)")
    parser.add_argument(
        "--rename",
        metavar="DIRNAME=condition,type,replicate",
        action="append",
        default=[],
        help="Override parsing for a non-standard top-level directory name (repeatable)",
    )
    args = parser.parse_args()

    rename_map = parse_rename(args.rename)
    groups = collect_files(args.fastq_root, rename_map)
    rows = build_rows(groups)

    if not rows:
        print("ERROR: no samples found; check fastq_root and naming convention", file=sys.stderr)
        sys.exit(1)

    header = ["sample_name", "sample_type", "replicate", "read_r1", "read_r2"]
    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        print("\t".join(header), file=out)
        for row in rows:
            print("\t".join(str(row[c]) for c in header), file=out)
    finally:
        if args.output != "-":
            out.close()

    print(f"Wrote {len(rows)} rows.", file=sys.stderr)


if __name__ == "__main__":
    main()
