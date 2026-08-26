#!/usr/bin/env python3
"""
Post-run assertions for individual unit tests.

Usage:
    python tests/e2e/unit_assertions.py <rule> <outdir> [--filtered]

    rule      One of: normalize, collect_consensus, qc_initial, qc_final
    outdir    Directory where the unit mini-Snakefile wrote its outputs
    --filtered  (normalize only) Expect that the CPM filter removed some TSSs
    --expect-status   (qc_initial, qc_final) Status every library should have
    --expect-reason   (qc_initial, qc_final) Gate name every StatusReason should name
"""

import csv
import os
import sys
import argparse

N_CS_TSSS  = 10   # rows in tss.consensus.bed fixture
N_IN_TSSS  = 5    # rows in merged_in.bed fixture
TOP_SIZES_MAX = 5 # tests/data/config.yaml default for tss_top_sizes_max

REPO     = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIXTURES = os.path.join(REPO, "tests", "data", "unit_fixtures")


def fail(msg: str) -> None:
    print(f"  FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"  ok  {msg}")


def count_lines(path: str, skip_header: bool = False) -> int:
    with open(path) as fh:
        n = sum(1 for line in fh if line.strip())
    return n - (1 if skip_header else 0)


def assert_normalize(outdir: str, filtered: bool, expect_kept: int = None) -> None:
    bed_path  = os.path.join(outdir, "tss.final.bed")
    raw_path  = os.path.join(outdir, "tss.final.raw.txt")
    norm_path = os.path.join(outdir, "tss.final.cpm.txt")
    nf_path   = os.path.join(outdir, "norm_factors.txt")

    for p in (bed_path, raw_path, norm_path, nf_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    n_final = count_lines(bed_path)

    if expect_kept is not None:
        if n_final != expect_kept:
            fail(f"tss.final.bed has {n_final} rows, expected exactly {expect_kept}")
        ok(f"Filter kept exactly {n_final} of {N_CS_TSSS} TSSs")
    elif filtered:
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


def assert_collect_consensus(outdir: str, expect_clusters: int = None) -> None:
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

    if expect_clusters is not None:
        if len(widths) != expect_clusters:
            fail(f"consensus has {len(widths)} TSSs, expected {expect_clusters}")
        ok(f"Consensus has exactly {expect_clusters} TSSs")

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


def check_status(path: str, basis: str, expect_status: str = None,
                 expect_reason: str = None) -> None:
    """Status, StatusReason and StatusBasis in one of the csRNA QC tables.

    StatusBasis is what says which table a Status came from, which is the difference
    between a FAIL that changed the run (the initial table, which
    collect_consensus_tss reads) and one that is only advice (the final table).
    """
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for col in ("Status", "StatusReason", "StatusBasis"):
        if col not in rows[0]:
            fail(f"{os.path.basename(path)} missing column: {col}")
    wrong_basis = [r["Sample"] for r in rows if r["StatusBasis"] != basis]
    if wrong_basis:
        fail(f"{os.path.basename(path)}: StatusBasis should be {basis!r} for "
             f"{wrong_basis}")
    ok(f"{os.path.basename(path)} has Status, StatusReason and StatusBasis={basis}")

    if expect_status is not None:
        wrong = [(r["Sample"], r["Status"]) for r in rows if r["Status"] != expect_status]
        if wrong:
            fail(f"{os.path.basename(path)}: expected Status {expect_status} for every "
                 f"library, got {wrong}")
        ok(f"every library is {expect_status}")
    if expect_reason is not None:
        wrong = [(r["Sample"], r["StatusReason"]) for r in rows
                 if expect_reason not in r["StatusReason"].split(",")]
        if wrong:
            fail(f"{os.path.basename(path)}: expected {expect_reason!r} among the failed "
                 f"gates for every library, got {wrong}")
        ok(f"every StatusReason names {expect_reason}")
    elif expect_status == "Ok":
        noisy = [(r["Sample"], r["StatusReason"]) for r in rows if r["StatusReason"]]
        if noisy:
            fail(f"{os.path.basename(path)}: passing libraries should have an empty "
                 f"StatusReason, got {noisy}")
        ok("passing libraries have an empty StatusReason")


def assert_qc_initial(outdir: str, expect_status: str = None,
                      expect_reason: str = None) -> None:
    qc_cs_path = os.path.join(outdir, "qc_initial_cs.txt")
    qc_in_path = os.path.join(outdir, "qc_initial_in.txt")
    cor_path   = os.path.join(outdir, "replicate_correlation.txt")
    for p in (qc_cs_path, qc_in_path, cor_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    # 3 csRNA samples in qc_initial_cs.txt
    n_cs = count_lines(qc_cs_path, skip_header=True)
    if n_cs != 3:
        fail(f"qc_initial_cs.txt has {n_cs} rows, expected 3")
    ok("qc_initial_cs.txt has 3 csRNA rows")

    # 2 input samples in qc_initial_in.txt
    n_in = count_lines(qc_in_path, skip_header=True)
    if n_in != 2:
        fail(f"qc_initial_in.txt has {n_in} rows, expected 2")
    ok("qc_initial_in.txt has 2 input rows")

    # Required columns present in qc_initial_cs.txt
    with open(qc_cs_path) as fh:
        hdr = fh.readline().strip().split("\t")
    for col in ("csFRiP", "csEnrichment", "TSSDetected", "group", "RawReads", "TrimmedReads",
                "AlignedReads", "BelowMapqReads", "FilteredOutReads"):
        if col not in hdr:
            fail(f"qc_initial_cs.txt missing column: {col}")
    ok("qc_initial_cs.txt has expected columns")

    # Verify shared-input pairing: condB_csrna1 should appear in qc_initial_cs.txt
    with open(qc_cs_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = [r["Sample"] for r in rows]
    if "condB_csrna1" not in samples:
        fail("condB_csrna1 not found in qc_initial_cs.txt; shared-input pairing may be broken")
    ok("condB_csrna1 present in qc_initial_cs.txt (shared-input pairing works)")

    check_status(qc_cs_path, "initial", expect_status, expect_reason)


def assert_qc_final(outdir: str, expect_status: str = None,
                    expect_reason: str = None) -> None:
    qc_cs_path = os.path.join(outdir, "qc_final_cs.txt")
    qc_in_path = os.path.join(outdir, "qc_final_in.txt")
    cor_path   = os.path.join(outdir, "replicate_correlation.txt")
    for p in (qc_cs_path, qc_in_path, cor_path):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    # Required columns
    with open(qc_cs_path) as fh:
        hdr = fh.readline().strip().split("\t")
    for col in (
        "NConsensusTSS", "NFinalTSS", "NFilteredTSS", "FinalTSSDetected",
        "csRiP", "sRiP", "csFRiP", "sFRiP", "csRNACappedPct", "csEnrichment", "sDepletion",
        "group", "RawReads", "TrimmedReads",
        "AlignedReads", "BelowMapqReads", "FilteredOutReads",
    ):
        if col not in hdr:
            fail(f"qc_final_cs.txt missing column: {col}")
    ok("qc_final_cs.txt has required columns")

    # replicate_correlation.txt must have both Stage values
    with open(cor_path) as fh:
        cor_rows = list(csv.DictReader(fh, delimiter="\t"))
    stages = {r["Stage"] for r in cor_rows}
    if "Initial" not in stages or "Final" not in stages:
        fail(f"replicate_correlation.txt missing Stage values (got: {stages})")
    ok("replicate_correlation.txt has both Initial and Final stages")

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

    check_status(qc_cs_path, "final-advisory", expect_status, expect_reason)


def assert_size_composition(outdir: str, srna_enabled: bool, top_enabled: bool) -> None:
    path = os.path.join(outdir, "tss.consensus.sizes.txt")
    if not os.path.exists(path):
        fail(f"Missing output: {path}")

    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != N_CS_TSSS:
        fail(f"tss.consensus.sizes.txt has {len(rows)} rows, expected {N_CS_TSSS}")

    if not (srna_enabled or top_enabled):
        # Both filters off: cluster names only, so no tag directory was read.
        if list(rows[0]) != ["TSS"]:
            fail(f"With no filter configured, expected a TSS column only, got {list(rows[0])}")
        ok(f"Filters disabled: {len(rows)} cluster names, no per-sample columns")
        return

    # The fixture holds every column the script can emit; a run writes .reads always,
    # .srna only when tss_srna_sizes is set, and .topn only when the top-lengths filter
    # is on. Check the column set, then the values for whatever was written.
    expected_path = os.path.join(FIXTURES, "tss.consensus.sizes.txt")
    with open(expected_path) as fh:
        expected = {r["TSS"]: r for r in csv.DictReader(fh, delimiter="\t")}
    got = {r["TSS"]: r for r in rows}
    if set(got) != set(expected):
        fail(f"cluster names differ from fixture: {set(got) ^ set(expected)}")

    cols = [c for c in rows[0] if c != "TSS"]
    suffixes = {c.rsplit(".", 1)[1] for c in cols}
    want = {"reads"} | ({"srna"} if srna_enabled else set())
    if top_enabled:
        # top1..topN are all precomputed so that tss_top_sizes_n can be retuned freely
        want |= {f"top{n}" for n in range(1, TOP_SIZES_MAX + 1)}
    if suffixes != want:
        fail(f"expected column kinds {sorted(want)}, got {sorted(suffixes)}")
    missing = [c for c in cols if c not in expected[next(iter(expected))]]
    if missing:
        fail(f"columns not present in the fixture: {missing}")
    for tss, exp in expected.items():
        for c in cols:
            if float(got[tss][c]) != float(exp[c]):
                fail(f"{tss} {c}: got {got[tss][c]}, expected {exp[c]}")
    ok(f"{', '.join(sorted(want))} counts match the fixture for all {len(rows)} clusters")


# aln_summary.awk against tests/data/unit_fixtures/aln_summary.sam, hand-counted:
#   14 records: 3 unmapped, 1 secondary, 1 supplementary, so 9 with a primary alignment
#   of those 9: 2 below MAPQ 30, 3 failing the rest (NM 5 > 2, rlen 10 < 15, no NM tag),
#   4 passing (including one 20M5D5M, whose reference length is 30 rather than 25)
ALN_SUMMARY_EXPECTED = {
    "ReadsSeen": 12,
    "Records": 14,
    "Unmapped": 3,
    "Secondary": 1,
    "Supplementary": 1,
    "PrimaryMapped": 9,
    "BelowMapq": 2,
    "FailedOtherFilters": 3,
    "Passed": 4,
    "UnmappedSampled": 2,
}
ALN_SUMMARY_MAPQ = {0: 1, 10: 1, 40: 7}


def assert_aln_summary(outdir: str) -> None:
    import gzip

    summary_path = os.path.join(outdir, "aln.raw.txt")
    passthrough   = os.path.join(outdir, "passthrough.sam")
    unmapped_fq   = os.path.join(outdir, "unmapped.sample.fastq.gz")
    for p in (summary_path, passthrough, unmapped_fq):
        if not os.path.exists(p):
            fail(f"Missing output: {p}")

    # The awk sits in the alignment pipe, so anything it changes about the stream would
    # change the BAM. This is the assertion that matters most.
    with open(os.path.join(FIXTURES, "aln_summary.sam"), "rb") as fh:
        want_stream = fh.read()
    with open(passthrough, "rb") as fh:
        got_stream = fh.read()
    if got_stream != want_stream:
        fail("the SAM that came out of aln_summary.awk is not the SAM that went in")
    ok("pass-through stream is byte-identical to the input")

    counts, mapq_hist = {}, {}
    with open(summary_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "MAPQ":
                mapq_hist[int(parts[1])] = int(parts[2])
            elif len(parts) == 2:
                counts[parts[0]] = parts[1]
    wrong = {k: (v, counts.get(k)) for k, v in ALN_SUMMARY_EXPECTED.items()
             if counts.get(k) != str(v)}
    if wrong:
        fail(f"aln.raw.txt counts wrong (expected, got): {wrong}")
    ok(f"all {len(ALN_SUMMARY_EXPECTED)} counts match the hand-counted fixture")

    # The three columns gather_stats.sh scrapes have to partition the reads, or the QC
    # table would show a loss that does not add up.
    if int(counts["PrimaryMapped"]) != (int(counts["BelowMapq"])
                                        + int(counts["FailedOtherFilters"])
                                        + int(counts["Passed"])):
        fail("PrimaryMapped != BelowMapq + FailedOtherFilters + Passed")
    if int(counts["ReadsSeen"]) != int(counts["Unmapped"]) + int(counts["PrimaryMapped"]):
        fail("ReadsSeen != Unmapped + PrimaryMapped")
    ok("counts partition the reads")

    if mapq_hist != ALN_SUMMARY_MAPQ:
        fail(f"MAPQ histogram is {mapq_hist}, expected {ALN_SUMMARY_MAPQ}")
    ok("MAPQ histogram matches")

    with gzip.open(unmapped_fq, "rt") as fh:
        fq_lines = [l for l in fh if l.strip()]
    if len(fq_lines) != 4 * ALN_SUMMARY_EXPECTED["UnmappedSampled"]:
        fail(f"unmapped sample has {len(fq_lines)} lines, expected "
             f"{4 * ALN_SUMMARY_EXPECTED['UnmappedSampled']}")
    ok(f"unmapped sample holds {ALN_SUMMARY_EXPECTED['UnmappedSampled']} reads")


RULES = {
    "normalize":          assert_normalize,
    "collect_consensus":  assert_collect_consensus,
    "qc_initial":         assert_qc_initial,
    "qc_final":           assert_qc_final,
    "size_composition":   assert_size_composition,
    "aln_summary":        assert_aln_summary,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rule", choices=list(RULES))
    parser.add_argument("outdir")
    parser.add_argument("--filtered", action="store_true",
                        help="For normalize: expect a filter removed TSSs")
    parser.add_argument("--expect-kept", type=int, default=None,
                        help="For normalize: exact number of TSSs expected in tss.final.bed")
    parser.add_argument("--srna-enabled", action="store_true",
                        help="For size_composition: expect .srna columns")
    parser.add_argument("--top-enabled", action="store_true",
                        help="For size_composition: expect top1..topN columns")
    parser.add_argument("--expect-clusters", type=int, default=None,
                        help="For collect_consensus: exact number of consensus TSSs")
    parser.add_argument("--expect-status", choices=["Ok", "FAIL"], default=None,
                        help="For qc_initial / qc_final: Status every library should have")
    parser.add_argument("--expect-reason", default=None,
                        help="For qc_initial / qc_final: a gate name every StatusReason "
                             "should list")
    args = parser.parse_args()

    fn = RULES[args.rule]
    if args.rule == "normalize":
        fn(args.outdir, args.filtered, args.expect_kept)
    elif args.rule == "size_composition":
        fn(args.outdir, args.srna_enabled, args.top_enabled)
    elif args.rule == "collect_consensus":
        fn(args.outdir, args.expect_clusters)
    elif args.rule == "aln_summary":
        fn(args.outdir)
    else:
        fn(args.outdir, args.expect_status, args.expect_reason)

    print()
