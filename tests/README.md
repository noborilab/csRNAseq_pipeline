# Test Suite

Tests live in `tests/` and are organised into four layers, from fastest to slowest:

| Layer | Command | Tools needed | Coverage |
|-------|---------|--------------|---------|
| `schema` | `./tests/run_tests.sh schema` | Python + snakemake | samples.schema.yaml, config.schema.yaml |
| `validate` | `./tests/run_tests.sh validate` | snakemake | DAG resolves; all rules present |
| `unit` | `./tests/run_tests.sh unit` | R + edgeR + rtracklayer | R scripts in isolation |
| `e2e` | `./tests/run_tests.sh e2e` | all pipeline tools | Full pipeline on synthetic data |
| `all` | `./tests/run_tests.sh all` | all | Everything, cheapest first |

Run from the repository root inside your Singularity container (or with HOMER, bfqutils, bwa, samtools, bedtools, and R on PATH).

## Synthetic fixtures

All committed test data lives in `tests/data/`:

| File | Description |
|------|-------------|
| `genome.fa` | 80 kb synthetic chromosome `testchr1` |
| `chrom.sizes` | One-line chrom sizes file |
| `samples.tsv` | 5-sample sheet with shared-input via `input_name` |
| `config.yaml` | Pipeline config pointing at synthetic fixtures |
| `mirnas.bed`, `trnas.bed` | One placeholder locus each |
| `reads/*.r1.fq.gz` | 500 synthetic reads per sample (~5–10 KB each) |
| `unit_fixtures/` | Hand-crafted TSVs for R-script unit tests |
| `unit_fixtures/aln_summary.sam` | 14 SAM records, one for each branch of `aln_summary.awk` (unmapped, secondary, supplementary, below MAPQ, too short, too many mismatches, no NM tag), with the hand-counted expectations in `unit_assertions.py` |
| `unit_fixtures/tagdir/` | Miniature HOMER tag directories whose read lengths are chosen so the read-size composition filter has a known right answer, alongside the expected `tss.consensus.sizes.txt` |
| `unit_fixtures/qc_initial_cs.txt` | Minimal QC table for the consensus QC gate; `condB_csrna1` is FAIL and is the only library seeing `TSS_7`, so excluding failed libraries is observable |
| `golden/` | Committed `qc_final_*.txt` from a default e2e run, compared value-by-value so a metric changing silently fails the suite. Refresh with `tests/scripts/update_golden.sh <e2e outdir>` once you have checked the diff is intended |
| `unit_fixtures/tss/` | Trimmed per-library `stats.txt` files carrying the enrichment threshold HOMER chose, for the `qc/min_log2_fold` gate. 1.686 passes a floor of 1; 0.144 and -0.680, both real values from an Arabidopsis panel, do not |

### Regenerating fixtures

Regenerate FASTQs and genome (only needed when changing the synthetic data design):

```bash
python tests/scripts/generate_synthetic_data.py
```

Regenerate unit-test TSV fixtures:

```bash
python tests/scripts/generate_unit_fixtures.py
```

Commit the resulting files. (Both scripts are deterministic, with a fixed seed, so regenerating without changing the design should give you identical output.)

### Refreshing the golden tables

Adding or removing a QC column fails the golden comparison by design, and `run_tests.sh`
deletes its own output directory on the way out, so the refresh needs a run you keep:

```bash
tmp=$(mktemp -d)
cfg=$(python tests/scripts/make_test_config.py "$tmp" "$tmp" "$tmp/index/genome.fa")
snakemake --cores 4 -s workflow/Snakefile --configfile "$cfg"
tests/scripts/update_golden.sh "$tmp"    # prints the diff before overwriting
```

Read the diff before committing it. A changed column set is expected after a deliberate
change; a changed *value* in a column you did not touch is a regression.

## Unit tests

Each file in `tests/unit/*.smk` is a minimal Snakefile with just a single rule in it.
Its inputs point at `tests/data/unit_fixtures/`, and its outputs go to a temporary
directory the runner supplies. The R script under test is left completely unmodified,
so it receives exactly the same `snakemake@*` object that it would in a production run.

To run a single unit manually:
```bash
tmp=$(mktemp -d)
snakemake -s tests/unit/normalize.smk \
  --configfile tests/data/config.yaml \
  --config test_outdir=$tmp \
  --cores 1
python tests/e2e/unit_assertions.py normalize $tmp
```

## End-to-end test

The e2e test runs the full pipeline twice:

1. **Default params** (`tss_min_cpm=0`): asserts `tss.final.bed` == `tss.consensus.bed` row-count.
2. **Strict filter** (`tss_min_cpm=10, tss_min_samples=2`): asserts `tss.final.bed` has fewer rows.

Outputs go to temporary directories and are cleaned up automatically.

## Adding a new test

- **New schema case**: add a TSV to `tests/schema/` and a test method in `test_schemas.py`.
- **New R-script unit test**: copy one of `tests/unit/*.smk`, wire in fixtures, add assertions to `unit_assertions.py`.
- **New e2e assertion**: add a check in `tests/e2e/assertions.py`.
