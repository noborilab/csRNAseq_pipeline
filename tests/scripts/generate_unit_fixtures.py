#!/usr/bin/env python3
"""
Generate hand-crafted unit test fixtures in tests/data/unit_fixtures/.

Run from the repository root:
    python tests/scripts/generate_unit_fixtures.py

These are small, human-readable TSVs committed alongside the synthetic FASTQ
data.  Unlike the FASTQ fixtures they are NOT derived from the synthetic genome;
they are designed to exercise specific code paths in the R scripts.
"""

import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
UF = os.path.join(REPO_ROOT, "tests", "data", "unit_fixtures")
os.makedirs(UF, exist_ok=True)

CHROM = "testchr1"

# ── Sample IDs ──────────────────────────────────────────────────────────────────
CS_IDS = ["condA_csrna1", "condA_csrna2", "condB_csrna1"]
IN_IDS = ["condA_input1", "condA_input2"]

# paired_in_ids for Snakemake params:
#   condA_csrna1 (rep1)  → condA_input1  (same name + rep)
#   condA_csrna2 (rep2)  → condA_input2  (same name + rep)
#   condB_csrna1 (rep1, input_name=condA) → condA_input1
PAIRED_IN_IDS = ["condA_input1", "condA_input2", "condA_input1"]

# ── Statistics fixtures ─────────────────────────────────────────────────────────
# Column order matches workflow/scripts/gather_stats.sh:
#   Sample  RawReads  TrimmedReads  AlignedReads  BelowMapqReads  FilteredOutReads
#   TotalReads  P20ReadLength  MedianReadLength  P80ReadLength  ModeReadLength
#   OrganelleReads  Freq1A  PosReads  NegReads
# The three alignment columns are chosen so that they account for the whole loss:
# AlignedReads - BelowMapqReads - FilteredOutReads == TotalReads, and
# TrimmedReads >= AlignedReads.
# The csRNA rows carry a wider gap between the median and P80 than the input rows,
# which is the shape difference the columns exist to show.
STATS_CS = [
    ("condA_csrna1", 550, 510, 505, 3, 2, 500, 25, 32, 42, 31, 0, 0.80, 400, 100),
    ("condA_csrna2", 530, 490, 486, 4, 2, 480, 26, 33, 43, 24, 0, 0.75, 360, 120),
    ("condB_csrna1", 570, 530, 526, 4, 2, 520, 24, 30, 40, 24, 0, 0.82, 426,  94),
]
STATS_IN = [
    ("condA_input1", 500, 460, 455, 3, 2, 450, 23, 27, 33, 24, 0, 0.50, 225, 225),
    ("condA_input2", 510, 470, 465, 3, 2, 460, 24, 28, 34, 24, 0, 0.45, 207, 253),
]

def write_stats(path, rows):
    with open(path, "w") as fh:
        fh.write("Sample\tRawReads\tTrimmedReads\tAlignedReads\tBelowMapqReads"
                 "\tFilteredOutReads\tTotalReads\tP20ReadLength\tMedianReadLength"
                 "\tP80ReadLength\tModeReadLength"
                 "\tOrganelleReads\tFreq1A\tPosReads\tNegReads\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

# ── CS TSS BED (consensus set, 10 TSSs all ≥150 bp) ────────────────────────────
# Designed so that with min_cpm=1, min_samples=2: TSS_3,4,7,8 are filtered out.
CS_TSS = [
    (CHROM, 1000,  1150,  "TSS_1",  0, "+"),
    (CHROM, 2000,  2200,  "TSS_2",  0, "+"),
    (CHROM, 3000,  3150,  "TSS_3",  0, "-"),
    (CHROM, 4000,  4200,  "TSS_4",  0, "+"),
    (CHROM, 5000,  5200,  "TSS_5",  0, "+"),
    (CHROM, 6000,  6150,  "TSS_6",  0, "+"),
    (CHROM, 7000,  7200,  "TSS_7",  0, "-"),
    (CHROM, 8000,  8150,  "TSS_8",  0, "+"),
    (CHROM, 9000,  9200,  "TSS_9",  0, "+"),
    (CHROM, 10000, 10150, "TSS_10", 0, "+"),
]

# ── IN TSS BED (input merged TSS set, 5 TSSs) ──────────────────────────────────
IN_TSS = [
    (CHROM, 1050,  1200,  "INTSS_1", 0, "+"),
    (CHROM, 3100,  3250,  "INTSS_2", 0, "-"),
    (CHROM, 5100,  5250,  "INTSS_3", 0, "+"),
    (CHROM, 7100,  7250,  "INTSS_4", 0, "-"),
    (CHROM, 20000, 20150, "INTSS_5", 0, "+"),
]

def write_bed(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

# ── Tag-count matrices ──────────────────────────────────────────────────────────
# quant_cs_cs: 10 CS-TSSs × 3 csRNA samples
# Designed so that cols sum to ~200 each (library size for CPM calculation)
QUANT_CS_CS = {
    "TSS_1":  [100, 90, 80],
    "TSS_2":  [50,  40, 30],
    "TSS_3":  [5,   0,  0],   # only 1 sample; filtered at min_cpm=1, min_samples=2
    "TSS_4":  [0,   0,  5],   # only 1 sample, filtered
    "TSS_5":  [10,  8,  0],   # 2 samples, kept at min_samples=2
    "TSS_6":  [3,   2,  1],
    "TSS_7":  [0,   0,  0],   # zero everywhere, filtered at min_cpm=1
    "TSS_8":  [0,   3,  0],   # only 1 sample, filtered
    "TSS_9":  [1,   1,  1],
    "TSS_10": [30,  25, 20],
}

# quant_cs_in: 10 CS-TSSs × 2 input samples (lower counts, input is background)
QUANT_CS_IN = {
    "TSS_1":  [20, 18],
    "TSS_2":  [10,  9],
    "TSS_3":  [1,   0],
    "TSS_4":  [0,   1],
    "TSS_5":  [2,   2],
    "TSS_6":  [1,   1],
    "TSS_7":  [0,   0],
    "TSS_8":  [0,   1],
    "TSS_9":  [0,   0],
    "TSS_10": [5,   4],
}

# quant_in_cs: 5 IN-TSSs × 3 csRNA samples
QUANT_IN_CS = {
    "INTSS_1": [5,  4,  3],
    "INTSS_2": [2,  1,  2],
    "INTSS_3": [3,  2,  1],
    "INTSS_4": [1,  1,  0],
    "INTSS_5": [0,  1,  1],
}

# quant_in_in: 5 IN-TSSs × 2 input samples (higher than csRNA, input enriched here)
QUANT_IN_IN = {
    "INTSS_1": [30, 25],
    "INTSS_2": [15, 12],
    "INTSS_3": [20, 18],
    "INTSS_4": [10,  8],
    "INTSS_5": [5,   6],
}

# quant_final_in: 10 CS-TSSs × 2 input samples (for qc_final_tss.R, HOMER format)
# same as quant_cs_in but in HOMER annotatePeaks.pl format
QUANT_FINAL_IN = QUANT_CS_IN  # reuse same values


def homer_header(peak_ids_unused, sample_ids, tag_label="Tag Count"):
    """Build a HOMER annotatePeaks.pl-style header row."""
    extra_cols = ["Chr", "Start", "End", "Strand", "Annotation"]
    count_cols = [
        f"{sid} {tag_label} ({i + 1} of {len(sample_ids)})"
        for i, sid in enumerate(sample_ids)
    ]
    return ["PeakID"] + extra_cols + count_cols


def homer_row(peak_id, chrom, start, end, strand, counts):
    extra = [chrom, str(start), str(end), strand, "NA"]
    return [peak_id] + extra + [str(c) for c in counts]


def write_homer_quant(path, tss_list, sample_ids, counts_dict):
    header = homer_header(None, sample_ids)
    rows = []
    for tss in tss_list:
        pid = tss[3]
        chrom, start, end, _, _, strand = tss[0], tss[1], tss[2], tss[3], tss[4], tss[5]
        counts = counts_dict[pid]
        rows.append(homer_row(pid, chrom, start, end, strand, counts))
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


# ── Per-sample TSS BEDs for collect_consensus_tss.R unit test ──────────────────
# condA rep1 and condA rep2 see TSSs 1-6; condB rep1 sees TSSs 1-4 and 6
# (TSS_5 not present in condB so condB acts as "different condition" replicates
# for testing that min_reps=2 keeps what's shared between condA reps)
PER_SAMPLE_TSS = {
    "condA_csrna1": CS_TSS[:6],          # TSS_1..TSS_6
    "condA_csrna2": CS_TSS[:6],          # same set
    # TSS_1..TSS_4, TSS_6, plus TSS_7 which only condB sees. condB is marked FAIL in
    # qc_initial_cs.txt, so excluding failed libraries from the consensus drops TSS_7 and
    # nothing else, which makes the flag observable in the output.
    "condB_csrna1": CS_TSS[:4] + [CS_TSS[5], CS_TSS[6]],
}


# ── Initial QC table, for the consensus QC gate ────────────────────────────────
QC_INITIAL_CS = [
    ("condA_csrna1", "Ok"),
    ("condA_csrna2", "Ok"),
    ("condB_csrna1", "FAIL"),
]


def write_qc_initial_cs(path):
    with open(path, "w") as fh:
        fh.write("Sample\tStatus\n")
        for sample, status in QC_INITIAL_CS:
            fh.write(f"{sample}\t{status}\n")

# ── Per-library HOMER stats files, for the enrichment-threshold gate ───────────
# The three csRNA thresholds are real values from an Arabidopsis panel: 1.686 is what
# HOMER picks for a sound library, while 0.144 and -0.680 are what it picks for libraries
# that have lost their capped signal. A negative threshold accepts clusters carrying less
# signal than their own input, so with the gate falling back to default_log2_fold (1) the
# first library passes and the other two fail.
LOG2_THRESHOLDS = {
    "condA_csrna1": 1.6858595532089,
    "condA_csrna2": 0.144,
    "condB_csrna1": -0.680,
    "condA_input1": 1.5,
    "condA_input2": 1.4,
}


def write_tss_stats(path, sample, threshold):
    """A trimmed findcsRNATSS.pl stats.txt: the shape the parser has to cope with,
    including the tab-indented threshold line and the -rna block that follows it with
    HOMER's no-op sentinel in place."""
    with open(path, "w") as fh:
        fh.write(
            f"cmd = findcsRNATSS.pl tagdir/{sample} -i tagdir/paired -o tss/{sample}"
            " -genome genome.fa -size 100 -gtf annotation.gtf\n\n"
            "Total csRNA reads: 500.0\n"
            "Total input reads: 450.0\n\n"
            "total putative TSS clusters\t12\n"
            "Valid TSS clusters\t10\n\n"
            "Fraction Promoter-Distal TSS clusters: 20.00%\n"
            "Fraction of stable transcript TSS clusters: na\n\n"
            "vs. Input (-i):\n"
            f"\tlog2 fold vs. input: {threshold}\n"
            "\tMaximum CDF difference: 0.81\n"
            "\tTotal TP (tss) regions: 15\n"
            "\tTotal FP (exon) regions: 5\n\n"
            "vs. RNA-seq (-rna):\n"
            "\tlog2 fold vs. rna: -10000000000\n"
            "\tMaximum CDF difference: -10000000000\n"
            "\tTotal TP (tss) regions: 0\n"
            "\tTotal FP (exon) regions: 0\n"
        )


# ── HOMER tag directories for tss_size_composition.R ───────────────────────────
# Reads placed inside the consensus clusters with chosen lengths, so both read-size
# filters have a known right answer. Each entry is (n_reads, length) and a list per TSS
# is a mixture; lengths within a TSS are distinct, as in a collapsed tag file.
#
# With srna sizes 21-25, max srna fraction 0.5, min reads 30:
#   TSS_3  87.5% small-RNA sized in all three libraries → dropped at min_samples 1 or 2
#   TSS_9  75% in condB only                            → dropped at 1, kept at 2
#   TSS_6  100% small-RNA sized but only 20 reads       → kept (below the read floor)
#   TSS_5  25% small-RNA sized                          → kept (below the fraction)
#   → 8 of 10 kept at min_samples 1, 9 at min_samples 2
#
# With top_sizes_n 2 and max top fraction 0.8:
#   TSS_3  two lengths only, so top-2 is 100%           → dropped in all three
#   TSS_7  two lengths only                             → dropped in all three
#   TSS_9  one length in condB after the siRNA mix      → dropped at 1, kept at 2
#   TSS_6  100% but under the read floor                → kept
#   others spread over 5-10 lengths (top-2 of 0.20-0.40) → kept
#   → 7 of 10 kept at min_samples 1, 8 at min_samples 2
SRNA_SIZES = [21, 22, 23, 24, 25]
TOP_N = 2      # the n the filter tests use
TOP_MAX = 5    # tss_top_sizes_max: top1..top5 are precomputed


def spread(n_per, lengths):
    """n_per reads at each of `lengths`, i.e. a heterogeneous cluster."""
    return [(n_per, ln) for ln in lengths]


TAG_MIX_COMMON = {
    "TSS_1":  spread(10, range(30, 40)),                    # 100 reads, top2 0.20
    "TSS_2":  spread(10, range(36, 41)),                    # 50 reads,  top2 0.40
    "TSS_3":  [(5, 35), (35, 24)],                          # 40 reads,  top2 1.00
    "TSS_4":  spread(5, range(28, 36)),                     # 40 reads,  top2 0.25
    "TSS_5":  spread(6, range(30, 35)) + [(10, 24)],        # 40 reads,  top2 0.40
    "TSS_6":  [(20, 24)],                                   # 20 reads,  under the floor
    "TSS_7":  [(20, 45), (20, 46)],                         # 40 reads,  top2 1.00
    "TSS_8":  spread(5, range(46, 54)),                     # 40 reads,  top2 0.25
    "TSS_9":  spread(4, range(32, 42)),                     # 40 reads,  top2 0.20
    "TSS_10": spread(7, range(57, 62)),                     # 35 reads,  top2 0.40
}
TAG_MIX = {
    "condA_csrna1": dict(TAG_MIX_COMMON),
    "condA_csrna2": dict(TAG_MIX_COMMON),
    # condB turns TSS_9 into a siRNA-dominated, single-length-dominated cluster
    "condB_csrna1": dict(TAG_MIX_COMMON, **{"TSS_9": [(10, 35), (30, 22)]}),
}


def write_tagdir(path, mix):
    """Write a HOMER-format <chrom>.tags.tsv: name, chr, pos, strand, count, len.

    Positions are spread across the cluster so no two rows share (pos, strand, len),
    which is what a real tag directory looks like after collapsing.
    """
    os.makedirs(path, exist_ok=True)
    rows = []
    for tss in CS_TSS:
        chrom, bed_start, bed_end, pid, _, strand = tss
        if pid not in mix:
            continue
        # BED is half-open and 0-based; tag positions are 1-based inclusive.
        first, last = bed_start + 1, bed_end
        strand_code = 0 if strand == "+" else 1
        offset = 0
        for n_reads, length in mix[pid]:
            for i in range(n_reads):
                pos = first + (offset % (last - first + 1))
                rows.append((chrom, pos, strand_code, 1, length))
                offset += 1
    rows.sort(key=lambda r: (r[1], r[2], r[4]))
    with open(os.path.join(path, f"{CHROM}.tags.tsv"), "w") as fh:
        for chrom, pos, strand_code, count, length in rows:
            fh.write(f"\t{chrom}\t{pos}\t{strand_code}\t{count}.0\t{length}\n")


def write_size_composition(path, tag_mix):
    """Hand-computed expected output of tss_size_composition.R for TAG_MIX.

    Holds every column the scripts can emit: .reads always, .srna when
    tss_srna_sizes is set, and .top1 .. .topTOP_MAX when the top-lengths filter is
    on. A run with only one filter enabled writes a subset of these columns, and the
    assertions compare whichever columns it produced.
    """
    samples = list(tag_mix)
    header = ["TSS"]
    for s in samples:
        header += [f"{s}.reads", f"{s}.srna"]
        header += [f"{s}.top{n}" for n in range(1, TOP_MAX + 1)]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for tss in CS_TSS:
            pid = tss[3]
            row = [pid]
            for s in samples:
                mix = tag_mix[s].get(pid, [])
                total = sum(n for n, _ in mix)
                srna = sum(n for n, ln in mix if ln in SRNA_SIZES)
                # lengths are distinct within a mix, so each entry is one length's total
                desc = sorted((n for n, _ in mix), reverse=True)
                row += [str(total), str(srna)]
                row += [str(sum(desc[:n])) for n in range(1, TOP_MAX + 1)]
            fh.write("\t".join(row) + "\n")


def main():
    # Stats
    write_stats(os.path.join(UF, "stats_cs.txt"), STATS_CS)
    write_stats(os.path.join(UF, "stats_in.txt"), STATS_IN)
    print("  stats_cs.txt, stats_in.txt")

    # BEDs
    write_bed(os.path.join(UF, "tss.consensus.bed"), CS_TSS)
    write_bed(os.path.join(UF, "merged_cs.bed"), CS_TSS)   # alias for qc_initial
    write_bed(os.path.join(UF, "merged_in.bed"), IN_TSS)
    print("  tss.consensus.bed, merged_cs.bed, merged_in.bed")

    # Per-sample TSS BEDs for collect_consensus_tss
    ps_dir = os.path.join(UF, "per_sample_tss")
    os.makedirs(ps_dir, exist_ok=True)
    for sid, tss_list in PER_SAMPLE_TSS.items():
        write_bed(os.path.join(ps_dir, f"{sid}.tss.bed"), tss_list)
    print(f"  per_sample_tss/*.tss.bed")

    # Tag directories + the hand-computed size composition they should produce
    for sid, mix in TAG_MIX.items():
        write_tagdir(os.path.join(UF, "tagdir", sid), mix)
    print("  tagdir/*/%s.tags.tsv" % CHROM)
    write_size_composition(os.path.join(UF, "tss.consensus.sizes.txt"), TAG_MIX)
    print("  tss.consensus.sizes.txt")
    write_qc_initial_cs(os.path.join(UF, "qc_initial_cs.txt"))
    print("  qc_initial_cs.txt")

    tss_dir = os.path.join(UF, "tss")
    os.makedirs(tss_dir, exist_ok=True)
    for sid, threshold in LOG2_THRESHOLDS.items():
        write_tss_stats(os.path.join(tss_dir, f"{sid}.stats.txt"), sid, threshold)
    print("  tss/*.stats.txt")

    # HOMER quant matrices
    write_homer_quant(os.path.join(UF, "tss.consensus.homer.raw.txt"), CS_TSS, CS_IDS, QUANT_CS_CS)
    write_homer_quant(os.path.join(UF, "quant_cs_in.txt"),              CS_TSS, IN_IDS, QUANT_CS_IN)
    write_homer_quant(os.path.join(UF, "quant_in_cs.txt"),              IN_TSS, CS_IDS, QUANT_IN_CS)
    write_homer_quant(os.path.join(UF, "quant_in_in.txt"),              IN_TSS, IN_IDS, QUANT_IN_IN)
    write_homer_quant(os.path.join(UF, "tss.final.in.raw.txt"),         CS_TSS, IN_IDS, QUANT_FINAL_IN)
    print("  tss.consensus.homer.raw.txt, quant_cs_in.txt, quant_in_cs.txt, quant_in_in.txt")
    print("  tss.final.in.raw.txt")

    print(f"Done. Unit fixtures written to {UF}")


if __name__ == "__main__":
    main()
