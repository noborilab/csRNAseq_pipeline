#!/usr/bin/env python3
"""
Generate a complete test config YAML with runtime paths substituted.

Usage:
    python tests/scripts/make_test_config.py <output_dir> <intermediate_dir> <genome_index> \
        [--min-cpm N] [--min-samples N] [--alignment-program PROG]

Prints the path to a temp YAML file (which the caller is responsible for removing).
"""

import sys
import os
import tempfile
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("intermediate_dir")
    p.add_argument("genome_index")
    p.add_argument("--min-cpm",     type=float, default=None)
    p.add_argument("--min-samples", type=int,   default=None)
    p.add_argument("--srna-sizes", default=None,
                   help="Comma-separated read lengths for the size-composition "
                        "filter, e.g. 21,22,23,24,25")
    p.add_argument("--srna-min-samples", type=int, default=None)
    p.add_argument("--max-top-sizes-fraction", type=float, default=None,
                   help="Enable the top-lengths filter at this fraction (1 disables)")
    p.add_argument("--top-sizes-n", type=int, default=None)
    p.add_argument("--alignment-program",
                   choices=["bwa-aln", "bwa-mem", "STAR", "bowtie2", "hisat2"],
                   default=None)
    args = p.parse_args()

    with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
        cfg = yaml.safe_load(fh)

    cfg["files"]["output_dir"]      = args.output_dir
    cfg["files"]["intermediate_dir"] = args.intermediate_dir
    cfg["program"]["genome_index"]  = args.genome_index

    # Make all fixture paths absolute so shadow-mode rules can resolve them
    # regardless of the CWD the rule runs from.
    def abs_fixture(rel):
        if rel and not os.path.isabs(rel):
            return os.path.join(REPO, rel)
        return rel

    cfg["sample_table"]  = abs_fixture(cfg.get("sample_table"))
    cfg["chrom_sizes"]   = abs_fixture(cfg.get("chrom_sizes"))
    cfg["genome_fasta"]  = abs_fixture(cfg.get("genome_fasta"))
    cfg["program"]["homer"]["genome"] = abs_fixture(cfg["program"]["homer"].get("genome"))
    if cfg.get("qc", {}).get("mirnas"):
        cfg["qc"]["mirnas"] = abs_fixture(cfg["qc"]["mirnas"])
    if cfg.get("qc", {}).get("trnas"):
        cfg["qc"]["trnas"]  = abs_fixture(cfg["qc"]["trnas"])

    if args.min_cpm is not None:
        cfg["filtering"]["tss_min_cpm"] = args.min_cpm
    if args.min_samples is not None:
        cfg["filtering"]["tss_min_samples"] = args.min_samples
    if args.srna_sizes is not None:
        cfg["filtering"]["tss_srna_sizes"] = [
            int(x) for x in args.srna_sizes.split(",") if x.strip()
        ]
    if args.srna_min_samples is not None:
        cfg["filtering"]["tss_srna_min_samples"] = args.srna_min_samples
    if args.max_top_sizes_fraction is not None:
        cfg["filtering"]["tss_max_top_sizes_fraction"] = args.max_top_sizes_fraction
    if args.top_sizes_n is not None:
        cfg["filtering"]["tss_top_sizes_n"] = args.top_sizes_n
    if args.alignment_program is not None:
        cfg["program"]["alignment_program"] = args.alignment_program
        # STAR's SA pre-indexing string must satisfy genomeSAindexNbases ≤ log2(L)/2-1.
        # For the 80 kb synthetic genome: log2(80000)/2 - 1 ≈ 7.1 → use 7.
        if args.alignment_program == "STAR":
            cfg["program"].setdefault("star", {})["genome_sa_index_nbases"] = 7

    # Fix relative fixture paths to be relative to repo root (already correct
    # since run_tests.sh cds to repo root before running snakemake)
    fd, path = tempfile.mkstemp(suffix=".yaml", prefix="csRNAseq_test_cfg_")
    with os.fdopen(fd, "w") as fh:
        yaml.dump(cfg, fh, default_flow_style=False)
    print(path)

if __name__ == "__main__":
    main()
