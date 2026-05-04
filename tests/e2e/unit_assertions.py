#!/usr/bin/env python3
"""
Post-run assertions for individual unit tests.

Usage:
    python tests/e2e/unit_assertions.py <rule> <outdir> [--filtered]

    rule      One of: normalize, collect_consensus, qc_initial, qc_final
    outdir    Directory where the unit mini-Snakefile wrote its outputs
    --filtered  (normalize only) Expect that the CPM filter removed some TSSs
"""

import csv
import os
import sys
import argparse

N_CS_TSSS  = 10   # rows in tss.consensus.bed fixture
N_IN_TSSS  = 5    # rows in merged_in.bed fixture


def fail(msg: str) -> None:
    print(f"  FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"  ok  {msg}")


def count_lines(path: str, skip_header: bool = False) -> int:
    with open(path) as fh:
        n = sum(1 for line in fh if line.strip())
    return n - (1 if skip_header else 0)


def assert_normalize(outdir: str, filtered: bool) -> None:
    bed_path  = os.path.join(outdir, "tss.final.bed")
    raw_path  = os.path.join(outdir, "tss.final.raw.txt")
    norm_path = os.path.join(outdir, "tss.final.cpm.txt")
    nf_path   = os.path.join(outdir, "norm_factors.txt")

    for p in (bed_path, raw_path, norm_path, nf_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    n_final = count_lines(bed_path)

    if filtered:
        if n_final >= N_CS_TSSS:
            fail(f"Expected CPM filter to remove TSSs, but tss.final.bed has {n_final} rows (same as consensus {N_CS_TSSS})")
        ok(f"CPM filter kept {n_final} of {N_CS_TSSS} TSSs")
    else:
        if n_final != N_CS_TSSS:
            fail(f"Default filter should be no-op: tss.final.bed has {n_final} rows, expected {N_CS_TSSS}")
        ok(f"Default no-op filter: {n_final} TSSs retained")

    # norm_factors.txt: 3 csRNA samples, norm.factors > 0 and finite
    with open(nf_path) as fh:
        nf_rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(nf_rows) != 3:
        fail(f"norm_factors.txt has {len(nf_rows)} rows, expected 3 csRNA samples")
    for r in nf_rows:
        val = float(r.get("norm.factors", "nan"))
        if not (0 < val < 1e6):
            fail(f"norm.factors out of range: {val}")
    ok("norm_factors.txt has 3 rows with finite factors")

    # raw and cpm tables have the same number of rows as the final BED
    n_raw  = count_lines(raw_path,  skip_header=True)
    n_cpm  = count_lines(norm_path, skip_header=True)
    if n_raw != n_final:
        fail(f"tss.final.raw.txt has {n_raw} rows but tss.final.bed has {n_final}")
    if n_cpm != n_final:
        fail(f"tss.final.cpm.txt has {n_cpm} rows but tss.final.bed has {n_final}")
    ok("raw and cpm tables have same row count as filtered BED")


def assert_collect_consensus(outdir: str) -> None:
    bed_path = os.path.join(outdir, "tss.consensus.bed")
    if not os.path.exists(bed_path):
        fail(f"Missing output: {bed_path}")

    # All features must be >= 150 bp
    widths = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t")
            widths.append(int(parts[2]) - int(parts[1]))
    if not widths:
        fail("tss.consensus.bed is empty")
    too_narrow = [w for w in widths if w < 150]
    if too_narrow:
        fail(f"{len(too_narrow)} TSSs are narrower than 150 bp: {too_narrow[:5]}")
    ok(f"All {len(widths)} TSSs are >= 150 bp")

    # Only valid chromosomes (testchr1 from the test chrom.sizes)
    chroms = set()
    with open(bed_path) as fh:
        for line in fh:
            if line.strip():
                chroms.add(line.split("\t")[0])
    if chroms - {"testchr1"}:
        fail(f"Unexpected chromosomes in output: {chroms - {'testchr1'}}")
    ok("All TSSs are on testchr1")

    # No overlapping intervals on the same strand
    entries = []
    with open(bed_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                entries.append((parts[0], int(parts[1]), int(parts[2]), parts[5]))
    by_strand = {}
    for chrom, start, end, strand in entries:
        by_strand.setdefault((chrom, strand), []).append((start, end))
    for key, ivs in by_strand.items():
        ivs.sort()
        for i in range(len(ivs) - 1):
            if ivs[i][1] > ivs[i + 1][0]:
                fail(f"Overlapping TSSs on {key}: {ivs[i]} and {ivs[i+1]}")
    ok("No overlapping TSSs on the same strand")


def assert_qc_initial(outdir: str) -> None:
    qc_cs_path = os.path.join(outdir, "qc_cs.txt")
    qc_in_path = os.path.join(outdir, "qc_in.txt")
    cor_path   = os.path.join(outdir, "replicate_correlation.txt")
    for p in (qc_cs_path, qc_in_path, cor_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    # 3 csRNA samples in qc_cs.txt
    n_cs = count_lines(qc_cs_path, skip_header=True)
    if n_cs != 3:
        fail(f"qc_cs.txt has {n_cs} rows, expected 3")
    ok("qc_cs.txt has 3 csRNA rows")

    # 2 input samples in qc_in.txt
    n_in = count_lines(qc_in_path, skip_header=True)
    if n_in != 2:
        fail(f"qc_in.txt has {n_in} rows, expected 2")
    ok("qc_in.txt has 2 input rows")

    # Required columns present in qc_cs.txt
    with open(qc_cs_path) as fh:
        hdr = fh.readline().strip().split("\t")
    for col in ("csFRiP", "csEnrichment", "TSSDetected"):
        if col not in hdr:
            fail(f"qc_cs.txt missing column: {col}")
    ok("qc_cs.txt has expected columns")

    # Verify shared-input pairing: condB_csrna1 should appear in qc_cs.txt
    with open(qc_cs_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = [r["Sample"] for r in rows]
    if "condB_csrna1" not in samples:
        fail("condB_csrna1 not found in qc_cs.txt — shared-input pairing may be broken")
    ok("condB_csrna1 present in qc_cs.txt (shared-input pairing works)")


def assert_qc_final(outdir: str) -> None:
    qc_cs_path = os.path.join(outdir, "qc_final_cs.txt")
    qc_in_path = os.path.join(outdir, "qc_final_in.txt")
    for p in (qc_cs_path, qc_in_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    # Required columns
    with open(qc_cs_path) as fh:
        hdr = fh.readline().strip().split("\t")
    for col in ("FinalFRiP", "FinalEnrichment", "NConsensusTSS", "NFinalTSS", "NFilteredTSS"):
        if col not in hdr:
            fail(f"qc_final_cs.txt missing column: {col}")
    ok("qc_final_cs.txt has required columns")

    # NFilteredTSS == NConsensusTSS - NFinalTSS for each row
    with open(qc_cs_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            nc    = int(float(row["NConsensusTSS"]))
            nf    = int(float(row["NFinalTSS"]))
            nfilt = int(float(row["NFilteredTSS"]))
            if nfilt != nc - nf:
                fail(f"Sample {row['Sample']}: NFilteredTSS({nfilt}) != NConsensusTSS({nc}) - NFinalTSS({nf})")
    ok("NFilteredTSS == NConsensusTSS - NFinalTSS for all rows")

    # 3 csRNA rows, 2 input rows
    n_cs = count_lines(qc_cs_path, skip_header=True)
    n_in = count_lines(qc_in_path, skip_header=True)
    if n_cs != 3:
        fail(f"qc_final_cs.txt has {n_cs} rows, expected 3")
    if n_in != 2:
        fail(f"qc_final_in.txt has {n_in} rows, expected 2")
    ok(f"qc_final_cs.txt ({n_cs} rows) and qc_final_in.txt ({n_in} rows) correct")


RULES = {
    "normalize":          assert_normalize,
    "collect_consensus":  assert_collect_consensus,
    "qc_initial":         assert_qc_initial,
    "qc_final":           assert_qc_final,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rule", choices=list(RULES))
    parser.add_argument("outdir")
    parser.add_argument("--filtered", action="store_true",
                        help="For normalize: expect CPM filter removed TSSs")
    args = parser.parse_args()

    fn = RULES[args.rule]
    if args.rule == "normalize":
        fn(args.outdir, args.filtered)
    else:
        fn(args.outdir)

    print()
