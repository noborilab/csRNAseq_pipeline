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
