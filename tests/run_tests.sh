#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# csRNA-seq pipeline test runner
#
# Usage:
#   ./tests/run_tests.sh [schema|validate|unit|e2e|all] [-k]
#
# Subcommands:
#   schema    Schema validation tests (snakemake + pandas required; fast)
#   validate  Snakemake --lint and -n dry-run (fast)
#   unit      R-script unit tests via mini Snakefiles (requires R + edgeR + rtracklayer)
#   e2e       Full end-to-end pipeline run (requires all pipeline tools)
#   all       Run all four in order, cheapest first (default)
#
# Options:
#   -k        Keep going on failure (default: stop at first failure)
#
# Run from the repository root:
#   ./tests/run_tests.sh all
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

KEEP_GOING=0
CMD="${1:-all}"
shift || true
for arg in "$@"; do
    [[ "$arg" == "-k" ]] && KEEP_GOING=1
done

PASS=0
FAIL=0
FAILED_TESTS=()

run_step() {
    local name="$1"
    shift
    echo ""
    echo "══════════════════════════════════════════════"
    echo "  $name"
    echo "══════════════════════════════════════════════"
    if "$@"; then
        echo "  → PASSED"
        PASS=$((PASS + 1))
    else
        echo "  → FAILED" >&2
        FAIL=$((FAIL + 1))
        FAILED_TESTS+=("$name")
        if [[ $KEEP_GOING -eq 0 ]]; then
            exit 1
        fi
    fi
}

# Helper: generate a temp config YAML with runtime paths substituted.
# Prints the path to the temp file; caller removes it.
make_cfg() {
    python tests/scripts/make_test_config.py "$@"
}

# ── schema ────────────────────────────────────────────────────────────────────
run_schema() {
    python -m pytest tests/schema/test_schemas.py -v --tb=short
}

# ── validate ──────────────────────────────────────────────────────────────────
run_validate() {
    local tmp cfg
    tmp=$(mktemp -d)
    cfg=$(make_cfg "$tmp" "$tmp" "$tmp/index/genome.fa")
    trap "rm -rf $tmp $cfg" EXIT

    echo "Running snakemake --lint (style warnings are informational, not fatal) ..."
    snakemake --lint \
        --cores 3 \
        -s workflow/Snakefile \
        --configfile "$cfg" \
        --quiet 2>&1 || true

    echo ""
    echo "Running snakemake -n (dry-run) ..."
    local dry_out
    dry_out=$(snakemake -n \
        --cores 3 \
        -s workflow/Snakefile \
        --configfile "$cfg" \
        2>&1)
    echo "$dry_out"

    # Verify all expected rules appear in the dry-run output
    local expected_rules=(
        trim align make_tagdir find_tss_initial make_raw_bedgraph
        gather_stats merge_initial_tss
        quantify_initial_cs_tss quantify_initial_in_tss qc_initial_tss
        collect_consensus_tss quantify_final_tss normalize_tss_quantification
        quantify_final_in_tss qc_final_tss generate_normalized_bw
    )
    local missing=()
    for rule in "${expected_rules[@]}"; do
        if ! echo "$dry_out" | grep -qE "^$rule[[:space:]]"; then
            missing+=("$rule")
        fi
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "ERROR: these rules were not in dry-run output: ${missing[*]}" >&2
        return 1
    fi
    echo ""
    echo "All ${#expected_rules[@]} expected rules found in dry-run."
}

# ── unit ──────────────────────────────────────────────────────────────────────
run_unit_smk() {
    local rule_name="$1"
    local smk="tests/unit/${rule_name}.smk"

    echo "--- $rule_name (default params) ---"
    local tmp cfg
    tmp=$(mktemp -d)
    cfg=$(make_cfg "$tmp" "$tmp" "$tmp/index/genome.fa")
    trap "rm -rf $tmp $cfg" EXIT

    snakemake -s "$smk" \
        --configfile "$cfg" \
        --config "test_outdir=$tmp" \
        --cores 1 \
        --quiet 2>&1
    python tests/e2e/unit_assertions.py "$rule_name" "$tmp"

    if [[ "$rule_name" == "normalize" ]]; then
        local tmp2 cfg2
        tmp2=$(mktemp -d)
        cfg2=$(make_cfg "$tmp2" "$tmp2" "$tmp2/index/genome.fa" --min-cpm 1 --min-samples 2)
        trap "rm -rf $tmp2 $cfg2" EXIT
        echo ""
        echo "--- normalize (CPM filter: min_cpm=1, min_samples=2) ---"
        snakemake -s "$smk" \
            --configfile "$cfg2" \
            --config "test_outdir=$tmp2" \
            --cores 1 \
            --quiet 2>&1
        python tests/e2e/unit_assertions.py normalize "$tmp2" --filtered
    fi
}

run_unit() {
    for rule in normalize collect_consensus qc_initial qc_final; do
        run_step "unit/$rule" run_unit_smk "$rule"
    done
}

# ── e2e ───────────────────────────────────────────────────────────────────────
run_e2e() {
    echo "--- e2e: default params (no-op CPM filter) ---"
    local tmp cfg
    tmp=$(mktemp -d)
    cfg=$(make_cfg "$tmp" "$tmp" "$tmp/index/genome.fa")
    trap "rm -rf $tmp $cfg" EXIT

    snakemake \
        --cores 4 \
        -s workflow/Snakefile \
        --configfile "$cfg" \
        2>&1
    python tests/e2e/assertions.py "$tmp"

    echo ""
    echo "--- e2e: CPM filter (min_cpm=10000, min_samples=3) ---"
    local tmp2 cfg2
    tmp2=$(mktemp -d)
    cfg2=$(make_cfg "$tmp2" "$tmp2" "$tmp2/index/genome.fa" --min-cpm 10000 --min-samples 3)
    trap "rm -rf $tmp2 $cfg2" EXIT

    snakemake \
        --cores 4 \
        -s workflow/Snakefile \
        --configfile "$cfg2" \
        2>&1
    python tests/e2e/assertions.py "$tmp2" --filter-pass
}

# ── dispatch ──────────────────────────────────────────────────────────────────
case "$CMD" in
    schema)
        run_step "schema" run_schema
        ;;
    validate)
        run_step "validate" run_validate
        ;;
    unit)
        run_unit
        ;;
    e2e)
        run_step "e2e" run_e2e
        ;;
    all)
        run_step "schema"   run_schema
        run_step "validate" run_validate
        run_unit
        run_step "e2e"      run_e2e
        ;;
    *)
        echo "Usage: $0 [schema|validate|unit|e2e|all] [-k]" >&2
        exit 1
        ;;
esac

echo ""
echo "══════════════════════════════════════════════"
echo "  Results: $PASS passed, $FAIL failed"
if [[ ${#FAILED_TESTS[@]} -gt 0 ]]; then
    echo "  Failed: ${FAILED_TESTS[*]}"
fi
echo "══════════════════════════════════════════════"
[[ $FAIL -eq 0 ]]
