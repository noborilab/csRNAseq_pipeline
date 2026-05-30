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
# Header: Sample  TotalReads  OrganelleReads  Freq1A  PosReads  NegReads
STATS_CS = [
    ("condA_csrna1", 500, 0, 0.80, 400, 100),
    ("condA_csrna2", 480, 0, 0.75, 360, 120),
    ("condB_csrna1", 520, 0, 0.82, 426,  94),
]
STATS_IN = [
    ("condA_input1", 450, 0, 0.50, 225, 225),
    ("condA_input2", 460, 0, 0.45, 207, 253),
]

def write_stats(path, rows):
    with open(path, "w") as fh:
        fh.write("Sample\tTotalReads\tOrganelleReads\tFreq1A\tPosReads\tNegReads\n")
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
    "condB_csrna1": CS_TSS[:4] + [CS_TSS[5]],  # TSS_1..TSS_4, TSS_6
}


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
