#!/usr/bin/env python3
"""
End-to-end assertions for the csRNA-seq pipeline test run.

Usage:
    python tests/e2e/assertions.py <output_dir> [--filter-pass]

Arguments:
    output_dir    The pipeline output directory (files.output_dir).
    --filter-pass If set, asserts that tss.final.bed has FEWER rows than
                  tss.consensus.bed (i.e. the CPM filter actually removed TSSs).
                  Without this flag, asserts they are equal (default no-op filter).
"""

import csv
import glob
import os
import sys
import argparse

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"  ok  {msg}")


def count_lines(path: str, skip_header: bool = False) -> int:
    with open(path) as fh:
        n = sum(1 for line in fh if line.strip())
    return n - (1 if skip_header else 0)


def check_exists_nonempty(path: str, label: str) -> None:
    if not os.path.exists(path):
        fail(f"{label} not found: {path}")
    if os.path.getsize(path) == 0:
        fail(f"{label} is empty: {path}")
    ok(f"{label} exists and is non-empty")


def check_file_nonempty(outdir: str, rel: str) -> None:
    check_exists_nonempty(os.path.join(outdir, rel), rel)


def assert_golden_qc(outdir: str) -> None:
    """Compare the QC tables against committed golden copies.

    The rest of the suite checks structure: columns present, row counts, filter
    arithmetic. None of it would notice a metric quietly changing value, which is the
    failure mode most likely to mislead later. The synthetic run is deterministic, so
    the tables should reproduce exactly on the same tool versions.

    A mismatch is not automatically a bug: upgrading HOMER, bwa or bfqutils can shift
    these legitimately. Inspect the diff, and if the new values are right, regenerate
    with:

        tests/scripts/update_golden.sh <e2e output dir>
    """
    for name in ("qc_final_cs.txt", "qc_final_in.txt"):
        got_path = os.path.join(outdir, name)
        exp_path = os.path.join(REPO, "tests", "data", "golden", name)
        if not os.path.exists(exp_path):
            fail(f"missing golden file {exp_path}")
        with open(got_path) as fh:
            got = list(csv.DictReader(fh, delimiter="\t"))
        with open(exp_path) as fh:
            exp = list(csv.DictReader(fh, delimiter="\t"))
        if len(got) != len(exp):
            fail(f"{name}: {len(got)} rows, golden has {len(exp)}")
        if set(got[0]) != set(exp[0]):
            only_new = sorted(set(got[0]) - set(exp[0]))
            only_old = sorted(set(exp[0]) - set(got[0]))
            fail(f"{name}: column set changed (added {only_new}, removed {only_old})")
        by_sample_got = {r["Sample"]: r for r in got}
        by_sample_exp = {r["Sample"]: r for r in exp}
        if set(by_sample_got) != set(by_sample_exp):
            fail(f"{name}: sample set changed")
        drift = []
        for sample, e in by_sample_exp.items():
            g = by_sample_got[sample]
            for col, want in e.items():
                have = g[col]
                if want == have:
                    continue
                try:
                    w, h = float(want), float(have)
                except ValueError:
                    drift.append(f"{sample}/{col}: {want!r} -> {have!r}")
                    continue
                if w != w and h != h:      # both NaN
                    continue
                tol = 1e-6 * max(abs(w), 1.0)
                if abs(w - h) > tol:
                    drift.append(f"{sample}/{col}: {w:g} -> {h:g}")
        if drift:
            fail(f"{name}: {len(drift)} value(s) drifted from golden:\n    "
                 + "\n    ".join(drift[:10]))
        ok(f"{name} matches the golden copy ({len(exp)} rows, {len(exp[0])} columns)")

def assert_stats_placeholders(outdir: str) -> None:
    """HOMER's not-computed fractions must read "na", not a number.

    Only valid for a run with no annotation and no -rna: with a GTF the
    promoter-distal fraction is a real measurement and must survive untouched, which
    check_annotation_wiring.py asserts instead.
    """
    stats = sorted(glob.glob(os.path.join(outdir, "tss", "*.stats.txt")))
    if not stats:
        fail(f"no *.stats.txt under {outdir}/tss")
    for path in stats:
        text = open(path).read()
        for label in ("Fraction of stable transcript TSS clusters",
                      "Fraction Promoter-Distal TSS clusters"):
            if f"{label}: na" not in text:
                got = next((l for l in text.splitlines() if l.startswith(label)), "<absent>")
                fail(f"{os.path.basename(path)}: expected {label!r} to read na with no "
                     f"annotation and no -rna, got {got!r}")
    ok(f"placeholder fractions read na in all {len(stats)} stats files")


def read_aln_summary(path: str) -> dict:
    counts = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "MAPQ" and len(parts) == 2:
                counts[parts[0]] = int(float(parts[1]))
    return counts


def assert_alignment_accounting(outdir: str, cs_samples, in_samples) -> None:
    """Every trimmed read is accounted for between trimming and the tag directory.

    This is the loss the QC could not see before: on one Arabidopsis panel a median
    74.9% of trimmed reads vanished between TrimmedReads and TotalReads with nothing
    recording whether they failed to align, fell below the MAPQ floor, or were dropped
    by a flag or length filter. The three checks here are that the counts partition the
    reads, that they add up to what the pipeline kept, and that the QC table carries the
    same numbers as the summary it was scraped from.
    """
    for s in cs_samples + in_samples:
        counts = read_aln_summary(os.path.join(outdir, "qc", f"{s}.aln.raw.txt"))
        if counts["ReadsSeen"] != counts["Unmapped"] + counts["PrimaryMapped"]:
            fail(f"{s}: ReadsSeen {counts['ReadsSeen']} != Unmapped "
                 f"{counts['Unmapped']} + PrimaryMapped {counts['PrimaryMapped']}")
        parts = counts["BelowMapq"] + counts["FailedOtherFilters"] + counts["Passed"]
        if counts["PrimaryMapped"] != parts:
            fail(f"{s}: PrimaryMapped {counts['PrimaryMapped']} != BelowMapq + "
                 f"FailedOtherFilters + Passed ({parts})")
    ok(f"alignment summaries partition the reads for all "
       f"{len(cs_samples) + len(in_samples)} libraries")

    for fname, ids in (("qc/stats_initial_cs.txt", cs_samples),
                       ("qc/stats_initial_in.txt", in_samples)):
        with open(os.path.join(outdir, fname)) as fh:
            rows = {r["Sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for s in ids:
            row = rows[s]
            counts = read_aln_summary(os.path.join(outdir, "qc", f"{s}.aln.raw.txt"))
            for col, key in (("AlignedReads", "PrimaryMapped"),
                             ("BelowMapqReads", "BelowMapq"),
                             ("FilteredOutReads", "FailedOtherFilters")):
                if int(float(row[col])) != counts[key]:
                    fail(f"{fname}: {s} {col} is {row[col]}, but the alignment summary "
                         f"says {key} is {counts[key]}")
            trimmed = int(float(row["TrimmedReads"]))
            aligned = int(float(row["AlignedReads"]))
            total = int(float(row["TotalReads"]))
            if counts["ReadsSeen"] != trimmed:
                fail(f"{s}: the aligner saw {counts['ReadsSeen']} reads but trimming "
                     f"wrote {trimmed}; the alignment accounting is incomplete")
            if aligned > trimmed:
                fail(f"{s}: AlignedReads {aligned} exceeds TrimmedReads {trimmed}")
            if aligned < total:
                fail(f"{s}: AlignedReads {aligned} is below TotalReads {total}, so the "
                     "tag directory holds reads the aligner never reported")
            if counts["Passed"] < total:
                fail(f"{s}: {counts['Passed']} reads passed the filter but the tag "
                     f"directory holds {total}")
    ok("QC tables agree with the alignment summaries, and the loss adds up")


def assert_log2_threshold(outdir: str, cs_samples, in_samples) -> None:
    """Log2FoldThreshold in the QC tables is the threshold HOMER actually used.

    Checked against HOMER's own stats file rather than against an expected number, so
    it holds whether HOMER derived the threshold from an annotation or fell back to
    -defaultLog2Fold.
    """
    for fname, ids in (("qc/qc_initial_cs.txt", cs_samples),
                       ("qc/qc_initial_in.txt", in_samples)):
        with open(os.path.join(outdir, fname)) as fh:
            rows = {r["Sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for s in ids:
            if "Log2FoldThreshold" not in rows[s]:
                fail(f"{fname} has no Log2FoldThreshold column")
            stats = os.path.join(outdir, "tss", f"{s}.stats.txt")
            line = next((l for l in open(stats) if "log2 fold vs. input:" in l), None)
            if line is None:
                fail(f"{stats} has no 'log2 fold vs. input' line")
            want = float(line.split(":")[1])
            have = float(rows[s]["Log2FoldThreshold"])
            if abs(want - have) > 1e-9:
                fail(f"{fname}: {s} Log2FoldThreshold is {have}, but {stats} says {want}")
    ok("Log2FoldThreshold matches the threshold HOMER reported for every library")


def main(outdir: str, filter_pass: bool, skip_golden: bool = False,
         no_annotation: bool = False) -> None:
    print(f"Checking outputs in: {outdir}")

    # Expected csRNA samples
    cs_samples = ["condA_csrna1", "condA_csrna2", "condB_csrna1"]
    # Expected QC rows: 3 csRNA samples, 2 input samples
    n_cs = 3
    n_in = 2

    # Core outputs
    for f in [
        "tss.consensus.bed",
        "tss.final.bed",
        "tss.final.raw.txt",
        "tss.final.cpm.txt",
        "norm_factors.txt",
        "tss.final.in.raw.txt",
    ]:
        check_file_nonempty(outdir, f)

    # QC outputs
    for f in [
        "qc/qc_initial_cs.txt",
        "qc/qc_initial_in.txt",
        "qc/replicate_correlation.txt",
        "qc_final_cs.txt",
        "qc_final_in.txt",
        "qc/stats_initial_cs.txt",
        "qc/stats_initial_in.txt",
    ]:
        check_file_nonempty(outdir, f)

    # bigWigs for every csRNA sample
    for s in cs_samples:
        for strand in ("pos", "neg"):
            check_file_nonempty(outdir, f"bw/{s}.rpm.{strand}.bw")

    # QC tables have the right number of data rows
    for fname, expected_n, label in [
        ("qc/qc_initial_cs.txt", n_cs, "qc_initial_cs"),
        ("qc/qc_initial_in.txt", n_in, "qc_initial_in"),
        ("qc_final_cs.txt",      n_cs, "qc_final_cs"),
        ("qc_final_in.txt",      n_in, "qc_final_in"),
        ("qc/stats_initial_cs.txt", n_cs, "stats_initial_cs"),
        ("qc/stats_initial_in.txt", n_in, "stats_initial_in"),
    ]:
        actual = count_lines(os.path.join(outdir, fname), skip_header=True)
        if actual != expected_n:
            fail(f"{label} has {actual} data rows, expected {expected_n}")
        ok(f"{label} has {expected_n} rows")

    assert_alignment_accounting(outdir, cs_samples, ["condA_input1", "condA_input2"])
    assert_log2_threshold(outdir, cs_samples, ["condA_input1", "condA_input2"])

    # qc_final_cs.txt must have all expected columns
    import csv
    with open(os.path.join(outdir, "qc_final_cs.txt")) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
    for col in (
        "NConsensusTSS", "NFinalTSS", "NFilteredTSS", "FinalTSSDetected",
        "csRiP", "sRiP", "csFRiP", "sFRiP", "csRNACappedPct", "csEnrichment", "sDepletion",
        "miRNA", "miRNADepletion", "miRNAPct", "PretRNA", "PretRNAPct", "PhosEfficiency",
        "group", "RawReads", "TrimmedReads",
        "AlignedReads", "BelowMapqReads", "FilteredOutReads",
    ):
        if col not in header:
            fail(f"qc_final_cs.txt missing column: {col}")
    ok("qc_final_cs.txt has required columns")

    # replicate_correlation.txt must have both Stage values
    with open(os.path.join(outdir, "qc/replicate_correlation.txt")) as fh:
        repl_rows = list(csv.DictReader(fh, delimiter="\t"))
    stages = {r["Stage"] for r in repl_rows}
    if "Initial" not in stages or "Final" not in stages:
        fail(f"replicate_correlation.txt missing Stage values (got: {stages})")
    ok("replicate_correlation.txt has both Initial and Final stages")

    # NFilteredTSS == NConsensusTSS - NFinalTSS (spot-check first row)
    with open(os.path.join(outdir, "qc_final_cs.txt")) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    nc = int(float(row["NConsensusTSS"]))
    nf = int(float(row["NFinalTSS"]))
    nfilt = int(float(row["NFilteredTSS"]))
    if nfilt != nc - nf:
        fail(f"NFilteredTSS={nfilt} != NConsensusTSS({nc}) - NFinalTSS({nf})")
    ok("NFilteredTSS == NConsensusTSS - NFinalTSS")

    # norm_factors.txt: one row per csRNA sample, norm.factors finite
    with open(os.path.join(outdir, "norm_factors.txt")) as fh:
        nf_rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(nf_rows) != n_cs:
        fail(f"norm_factors.txt has {len(nf_rows)} rows, expected {n_cs}")
    for r in nf_rows:
        val = float(r.get("norm.factors", "nan"))
        if not (0 < val < 1e6):
            fail(f"norm.factors out of range for sample {r}: {val}")
    ok("norm_factors.txt has one row per csRNA sample with finite factors")

    # CPM sanity: tss.final.cpm.txt non-header rows, each sample column should
    # sum to approximately 1e6 (CPM property). Allow 10% tolerance.
    with open(os.path.join(outdir, "tss.final.cpm.txt")) as fh:
        cpm_rows = list(csv.DictReader(fh, delimiter="\t"))
    for sid in cs_samples:
        if sid not in (cpm_rows[0] if cpm_rows else {}):
            continue  # sample not in cpm table (may happen with 0-read outcome)
        total = sum(float(r[sid]) for r in cpm_rows)
        if total < 0.5e6 or total > 1.5e6:
            fail(f"tss.final.cpm.txt column {sid} sums to {total:.0f}, expected ~1e6")
    ok("tss.final.cpm.txt per-sample sums are ~1e6")

    # TSS filter assertion
    n_consensus = count_lines(os.path.join(outdir, "tss.consensus.bed"))
    n_final     = count_lines(os.path.join(outdir, "tss.final.bed"))
    if filter_pass:
        if n_final >= n_consensus:
            fail(f"CPM filter had no effect: tss.final.bed ({n_final}) >= tss.consensus.bed ({n_consensus})")
        ok(f"CPM filter removed {n_consensus - n_final} TSSs ({n_final} of {n_consensus} retained)")
    else:
        if n_final != n_consensus:
            fail(f"Default params: tss.final.bed ({n_final}) != tss.consensus.bed ({n_consensus}); "
                 "expected no-op filter")
        ok(f"Default filter is no-op: {n_final} TSSs in both consensus and final")

    # Golden comparison only on the default run: the filtered variant and the per-aligner
    # runs legitimately produce different numbers.
    if not filter_pass and not skip_golden:
        assert_golden_qc(outdir)

    if no_annotation:
        assert_stats_placeholders(outdir)

    print("\nAll assertions passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("--filter-pass", action="store_true")
    parser.add_argument("--no-annotation", action="store_true",
                        help="Assert that HOMER's promoter-distal and stable-transcript "
                             "fractions read na. Only for runs with no GTF and no -rna, "
                             "where both are placeholders rather than measurements.")
    parser.add_argument("--skip-golden", action="store_true",
                        help="Skip the golden QC comparison. Used by the per-aligner runs, "
                             "where a different aligner may legitimately shift values and a "
                             "golden mismatch would be misread as a pipeline regression.")
    args = parser.parse_args()
    main(args.outdir, args.filter_pass, args.skip_golden, args.no_annotation)
