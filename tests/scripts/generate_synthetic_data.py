#!/usr/bin/env python3
"""
Generate synthetic test fixtures for the csRNA-seq pipeline test suite.

Run from the repository root:
    python tests/scripts/generate_synthetic_data.py

All outputs land in tests/data/ and should be committed.  The script is
deterministic (fixed seed) so re-running produces byte-identical FASTQs,
BEDs, and FASTA.
"""

import gzip
import os
import random
import string

SEED = 42
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_DIR = os.path.join(REPO_ROOT, "tests", "data")
READS_DIR = os.path.join(OUT_DIR, "reads")

# ── Genome parameters ──────────────────────────────────────────────────────────
CHROM = "testchr1"
GENOME_LEN = 80_000
READ_LEN = 70
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_PROB = 0.30   # fraction of reads that have adapter appended (before truncation)
QUAL = "I" * READ_LEN  # all Q=40

# ── TSS positions (0-based start on + strand) ─────────────────────────────────
# 15 evenly spaced TSS clusters well inside the chromosome
TSS_POSITIONS = [3000 + i * 5000 for i in range(15)]  # 3k, 8k, 13k, … 73k

# ── Read-count parameters ─────────────────────────────────────────────────────
# For a SINGLE-END csRNA sample: 500 reads, 80% from TSS regions
# For a SINGLE-END input  sample: 500 reads, 30% from TSS regions
N_READS = 500
CS_TSS_FRAC = 0.80
IN_TSS_FRAC = 0.30
TSS_SPREAD = 10   # reads start within ±SPREAD bp of TSS

SAMPLES = {
    # name                : (sample_type, replicate, random_seed_offset)
    "condA_csrna1": ("csrna", 1, 0),
    "condA_csrna2": ("csrna", 2, 1),
    "condB_csrna1": ("csrna", 1, 2),
    "condA_input1": ("input", 1, 3),
    "condA_input2": ("input", 2, 4),
}

MINUS_FRAC = 0.30  # fraction of reads reverse-complemented → minus-strand alignments

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_COMP)


def make_genome(rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=GENOME_LEN))


def make_read(genome: str, start: int, length: int, add_adapter: bool, minus: bool) -> str:
    end = start + length
    if end > len(genome):
        end = len(genome)
        start = max(0, end - length)
    seq = genome[start:end]
    if minus:
        seq = reverse_complement(seq)
    if add_adapter:
        seq = seq + ADAPTER
    return seq[:length]


def make_fastq_reads(genome: str, sample_type: str, n_reads: int, seed: int) -> list[str]:
    rng = random.Random(SEED * 100 + seed)
    tss_frac = CS_TSS_FRAC if sample_type == "csrna" else IN_TSS_FRAC
    records = []
    for i in range(n_reads):
        use_tss = rng.random() < tss_frac
        if use_tss:
            tss_pos = rng.choice(TSS_POSITIONS)
            offset = rng.randint(-TSS_SPREAD, TSS_SPREAD)
            start = max(0, tss_pos + offset)
        else:
            start = rng.randint(0, GENOME_LEN - READ_LEN - 1)
        minus = rng.random() < MINUS_FRAC
        add_adapter = rng.random() < ADAPTER_PROB
        seq = make_read(genome, start, READ_LEN, add_adapter, minus)
        records.append(f"@synth_{sample_type}_{i}\n{seq}\n+\n{QUAL}\n")
    return records


def write_fastq_gz(path: str, records: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.writelines(records)


def write_genome(genome: str) -> None:
    fa_path = os.path.join(OUT_DIR, "genome.fa")
    sizes_path = os.path.join(OUT_DIR, "chrom.sizes")
    with open(fa_path, "w") as fh:
        fh.write(f">{CHROM}\n")
        # 60 bp per line
        for i in range(0, len(genome), 60):
            fh.write(genome[i : i + 60] + "\n")
    with open(sizes_path, "w") as fh:
        fh.write(f"{CHROM}\t{GENOME_LEN}\n")
    print(f"  genome.fa ({GENOME_LEN // 1000} kb), chrom.sizes")


def write_beds() -> None:
    # mirnas.bed: 1 locus between TSS_0 (3000) and TSS_1 (8000), not overlapping either
    mirna_path = os.path.join(OUT_DIR, "mirnas.bed")
    with open(mirna_path, "w") as fh:
        fh.write(f"{CHROM}\t5200\t5600\tMIR001\t0\t+\n")

    # trnas.bed: 1 locus between TSS_1 and TSS_2, not overlapping TSSs
    trna_path = os.path.join(OUT_DIR, "trnas.bed")
    with open(trna_path, "w") as fh:
        fh.write(f"{CHROM}\t10200\t10600\ttRNA001\t0\t-\n")

    print("  mirnas.bed, trnas.bed")


def write_samples_tsv() -> None:
    tsv = os.path.join(OUT_DIR, "samples.tsv")
    lines = [
        "sample_name\tsample_type\treplicate\tread_r1\tinput_name",
        "condA\tcsrna\t1\ttests/data/reads/condA_csrna1.r1.fq.gz\t",
        "condA\tcsrna\t2\ttests/data/reads/condA_csrna2.r1.fq.gz\t",
        "condA\tinput\t1\ttests/data/reads/condA_input1.r1.fq.gz\t",
        "condA\tinput\t2\ttests/data/reads/condA_input2.r1.fq.gz\t",
        # condB csrna rep 1 shares condA's input (exercises input_name feature)
        "condB\tcsrna\t1\ttests/data/reads/condB_csrna1.r1.fq.gz\tcondA",
    ]
    with open(tsv, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  samples.tsv")


def write_config() -> None:
    cfg = os.path.join(OUT_DIR, "config.yaml")
    content = """\
### Test configuration for the csRNA-seq pipeline
# output_dir and intermediate_dir are overridden at runtime by run_tests.sh

sample_table: "tests/data/samples.tsv"

chrom_sizes: "tests/data/chrom.sizes"

genome_fasta: "tests/data/genome.fa"

program:
  # genome_index is overridden at runtime by run_tests.sh (points to tmpdir)
  genome_index: "PLACEHOLDER_SET_BY_RUNNER"
  bfqutils:
    sequencing_adapter: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    use_pigz: False
  alignment_threads: 1
  alignment_program: "bwa-aln"
  keep_trimmed_fastq: False
  keep_bam: False
  bwa_aln:
    keep_sai: False
  # genomeSAindexNbases must be ≤ log2(GenomeLen)/2-1; for the 80 kb test
  # genome that is floor(log2(80000)/2-1) = 7 (production default is 14).
  star:
    read_files_command: "gzip -dc"
    align_ends_type: "EndToEnd"
    align_intron_max: 1
    out_filter_multimap_nmax: 1000
    genome_sa_index_nbases: 7
  bowtie2:
    extra_params: "--very-sensitive"
  hisat2:
    extra_params: "--very-sensitive"
  homer:
    genome: "tests/data/genome.fa"
    tagdir:
      extra_args: "-checkGC"
    tss:
      size: 100
      ntagthreshold: 3
      extra_args: ""

filtering:
  alignment_mapq: 0
  include_flags: 0
  exclude_flags: 2308
  max_read_length: 70
  min_alignment_length: 15
  max_mismatch: 999
  tss_min_reps: 1
  tss_min_cpm: 0
  tss_min_samples: 1
  mask_srnas_from_bw: False

files:
  output_dir: "PLACEHOLDER_SET_BY_RUNNER"
  intermediate_dir: "PLACEHOLDER_SET_BY_RUNNER"

qc:
  organelle_chroms: []
  mirnas: "tests/data/mirnas.bed"
  trnas: "tests/data/trnas.bed"
  min_cs_frip: 0.0
  min_pct_nuclear: 0.0
"""
    with open(cfg, "w") as fh:
        fh.write(content)
    print("  config.yaml")


def main() -> None:
    rng = random.Random(SEED)
    print("Generating synthetic genome ...")
    genome = make_genome(rng)

    os.makedirs(READS_DIR, exist_ok=True)
    write_genome(genome)
    write_beds()
    write_samples_tsv()
    write_config()

    print("Generating synthetic reads ...")
    for sample_id, (sample_type, replicate, seed_offset) in SAMPLES.items():
        reads = make_fastq_reads(genome, sample_type, N_READS, seed_offset)
        path = os.path.join(READS_DIR, f"{sample_id}.r1.fq.gz")
        write_fastq_gz(path, reads)
        size = os.path.getsize(path)
        print(f"  {sample_id}.r1.fq.gz  ({size // 1024} KB)")

    print("Done.  All fixtures committed under tests/data/.")


if __name__ == "__main__":
    main()
