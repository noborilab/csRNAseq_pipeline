"""
Schema validation tests for samples.schema.yaml and config.schema.yaml.

Run from the repository root:
    python -m pytest tests/schema/test_schemas.py -v

Requires: snakemake (already a pipeline dependency).
"""

import os
import sys
import pytest
import pandas as pd
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

SAMPLES_SCHEMA = os.path.join(REPO, "workflow", "schemas", "samples.schema.yaml")
CONFIG_SCHEMA  = os.path.join(REPO, "workflow", "schemas", "config.schema.yaml")


def validate_samples(tsv_path):
    from snakemake.utils import validate
    df = pd.read_table(tsv_path, dtype=str)
    # Mirror the Snakefile's type handling: replicate must be int for the schema
    if "replicate" in df.columns:
        df["replicate"] = pd.to_numeric(df["replicate"].fillna("1"), errors="coerce").astype("Int64")
    df = df.fillna("")
    validate(df, SAMPLES_SCHEMA)


def validate_config(yaml_path):
    from snakemake.utils import validate
    with open(yaml_path) as fh:
        cfg = yaml.safe_load(fh)
    validate(cfg, CONFIG_SCHEMA)


# ── samples.schema.yaml ────────────────────────────────────────────────────────

class TestSamplesSchema:
    def test_valid_sheet_passes(self):
        """A well-formed sample sheet validates without errors."""
        validate_samples(os.path.join(REPO, "tests", "schema", "samples_valid.tsv"))

    def test_shared_input_passes(self):
        """A sheet using input_name (shared input) validates."""
        validate_samples(os.path.join(REPO, "tests", "schema", "samples_shared_input.tsv"))

    def test_example_sheet_passes(self):
        """The committed config/samples_example.tsv validates."""
        validate_samples(os.path.join(REPO, "config", "samples_example.tsv"))

    def test_synthetic_test_sheet_passes(self):
        """The synthetic test sheet (used for e2e) validates."""
        validate_samples(os.path.join(REPO, "tests", "data", "samples.tsv"))

    def test_missing_required_column_fails(self):
        """A sheet missing read_r1 must fail validation."""
        with pytest.raises(Exception):
            validate_samples(os.path.join(REPO, "tests", "schema", "samples_missing_col.tsv"))

    def test_bad_sample_type_fails(self):
        """A sheet with sample_type='CSRNA' (wrong case) must fail validation."""
        with pytest.raises(Exception):
            validate_samples(os.path.join(REPO, "tests", "schema", "samples_bad_type.tsv"))


# ── config.schema.yaml ─────────────────────────────────────────────────────────

class TestConfigSchema:
    def test_test_config_passes(self):
        """tests/data/config.yaml validates against config.schema.yaml."""
        validate_config(os.path.join(REPO, "tests", "data", "config.yaml"))

    def test_example_config_passes(self):
        """config/config.yaml validates against config.schema.yaml."""
        validate_config(os.path.join(REPO, "config", "config.yaml"))

    def test_new_filter_keys_accepted(self):
        """tss_min_cpm and tss_min_samples are accepted in config."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_min_cpm"] = 1.5
        cfg["filtering"]["tss_min_samples"] = 3
        validate(cfg, CONFIG_SCHEMA)

    def test_srna_size_keys_accepted(self):
        """The read-size composition filter keys are accepted in config."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_srna_sizes"] = [21, 22, 23, 24, 25]
        cfg["filtering"]["tss_max_srna_fraction"] = 0.5
        cfg["filtering"]["tss_srna_min_reads"] = 30
        cfg["filtering"]["tss_srna_min_samples"] = 2
        validate(cfg, CONFIG_SCHEMA)

    def test_srna_sizes_must_be_integers(self):
        """tss_srna_sizes given as strings must fail validation."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_srna_sizes"] = ["21", "22"]
        with pytest.raises(Exception):
            validate(cfg, CONFIG_SCHEMA)

    def test_max_srna_fraction_is_bounded(self):
        """tss_max_srna_fraction above 1 must fail validation."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_max_srna_fraction"] = 1.5
        with pytest.raises(Exception):
            validate(cfg, CONFIG_SCHEMA)

    def test_top_sizes_keys_accepted(self):
        """The top-lengths filter keys are accepted in config."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_top_sizes_n"] = 3
        cfg["filtering"]["tss_max_top_sizes_fraction"] = 0.8
        validate(cfg, CONFIG_SCHEMA)

    def test_top_sizes_n_must_be_at_least_one(self):
        """tss_top_sizes_n of 0 must fail validation."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_top_sizes_n"] = 0
        with pytest.raises(Exception):
            validate(cfg, CONFIG_SCHEMA)

    def test_tss_annotation_and_program_keys_accepted(self):
        """The optional GTF, program choice and HOMER knobs validate."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["program"]["homer"]["tss"]["gtf"] = "tests/data/genome.fa"
        cfg["program"]["homer"]["tss"]["program"] = "findcsRNATSR.pl"
        cfg["program"]["homer"]["tss"]["pseudo_count"] = 3
        cfg["program"]["homer"]["tss"]["default_log2_fold"] = 2
        cfg["program"]["homer"]["tss"]["local_fold"] = 4
        validate(cfg, CONFIG_SCHEMA)

    def test_tss_program_rejects_unknown(self):
        """Only the two known HOMER callers are accepted."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["program"]["homer"]["tss"]["program"] = "findPeaks"
        with pytest.raises(Exception):
            validate(cfg, CONFIG_SCHEMA)

    def test_consensus_qc_gate_and_top_max_accepted(self):
        """exclude_failed_from_consensus and tss_top_sizes_max validate."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["exclude_failed_from_consensus"] = True
        cfg["filtering"]["tss_top_sizes_max"] = 8
        validate(cfg, CONFIG_SCHEMA)

    def test_max_top_sizes_fraction_is_bounded(self):
        """tss_max_top_sizes_fraction above 1 must fail validation."""
        from snakemake.utils import validate
        with open(os.path.join(REPO, "tests", "data", "config.yaml")) as fh:
            cfg = yaml.safe_load(fh)
        cfg["filtering"]["tss_max_top_sizes_fraction"] = 1.2
        with pytest.raises(Exception):
            validate(cfg, CONFIG_SCHEMA)
