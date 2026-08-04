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
#   star      End-to-end with STAR alignment
#   bowtie2   End-to-end with bowtie2 alignment
#   bwa-mem   End-to-end with bwa-mem alignment
#   hisat2    End-to-end with hisat2 alignment
#   all       Run all in order, cheapest first (default)
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
        collect_consensus_tss quantify_final_tss
        tss_size_composition_library tss_size_composition
        normalize_tss_quantification run_info
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
    # Every step records into rc rather than aborting, so a failure in one
    # invocation is still reported when later invocations follow it.
    local rc=0

    echo "--- $rule_name (default params) ---"
    local tmp cfg
    tmp=$(mktemp -d)
    cfg=$(make_cfg "$tmp" "$tmp" "$tmp/index/genome.fa")
    trap "rm -rf $tmp $cfg" EXIT

    snakemake -s "$smk" \
        --configfile "$cfg" \
        --config "test_outdir=$tmp" \
        --cores 1 \
        --quiet 2>&1 || rc=1
    python tests/e2e/unit_assertions.py "$rule_name" "$tmp" || rc=1

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
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py normalize "$tmp2" --filtered || rc=1

        # Size-composition filter. The fixture is built so TSS_3 is small-RNA
        # dominated in all three libraries and TSS_9 in condB only, so any-library
        # voting drops both and requiring two libraries drops only TSS_3.
        local tmp3 cfg3
        tmp3=$(mktemp -d)
        cfg3=$(make_cfg "$tmp3" "$tmp3" "$tmp3/index/genome.fa" --srna-sizes 21,22,23,24,25)
        trap "rm -rf $tmp3 $cfg3" EXIT
        echo ""
        echo "--- normalize (size filter: srna sizes 21-25, min_samples=1) ---"
        snakemake -s "$smk" \
            --configfile "$cfg3" \
            --config "test_outdir=$tmp3" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py normalize "$tmp3" --expect-kept 8 || rc=1

        local tmp4 cfg4
        tmp4=$(mktemp -d)
        cfg4=$(make_cfg "$tmp4" "$tmp4" "$tmp4/index/genome.fa" \
            --srna-sizes 21,22,23,24,25 --srna-min-samples 2)
        trap "rm -rf $tmp4 $cfg4" EXIT
        echo ""
        echo "--- normalize (size filter: srna sizes 21-25, min_samples=2) ---"
        snakemake -s "$smk" \
            --configfile "$cfg4" \
            --config "test_outdir=$tmp4" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py normalize "$tmp4" --expect-kept 9 || rc=1

        # Top-lengths filter. TSS_3 and TSS_7 have only two lengths each in every
        # library, and condB's TSS_9 collapses to one, so any-library voting drops
        # three clusters and requiring two libraries drops two.
        local tmp6 cfg6
        tmp6=$(mktemp -d)
        cfg6=$(make_cfg "$tmp6" "$tmp6" "$tmp6/index/genome.fa" --max-top-sizes-fraction 0.8)
        trap "rm -rf $tmp6 $cfg6" EXIT
        echo ""
        echo "--- normalize (top-lengths filter: n=2, max 0.8, min_samples=1) ---"
        snakemake -s "$smk" \
            --configfile "$cfg6" \
            --config "test_outdir=$tmp6" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py normalize "$tmp6" --expect-kept 7 || rc=1

        local tmp7 cfg7
        tmp7=$(mktemp -d)
        cfg7=$(make_cfg "$tmp7" "$tmp7" "$tmp7/index/genome.fa" \
            --max-top-sizes-fraction 0.8 --srna-min-samples 2)
        trap "rm -rf $tmp7 $cfg7" EXIT
        echo ""
        echo "--- normalize (top-lengths filter: n=2, max 0.8, min_samples=2) ---"
        snakemake -s "$smk" \
            --configfile "$cfg7" \
            --config "test_outdir=$tmp7" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py normalize "$tmp7" --expect-kept 8 || rc=1
    fi

    if [[ "$rule_name" == "collect_consensus" ]]; then
        # condB_csrna1 is FAIL in the QC fixture and is the only library seeing TSS_7, so
        # the union has 7 clusters by default and 6 with failed libraries excluded.
        python tests/e2e/unit_assertions.py collect_consensus "$tmp" --expect-clusters 7 || rc=1
        local tmp9 cfg9
        tmp9=$(mktemp -d)
        cfg9=$(make_cfg "$tmp9" "$tmp9" "$tmp9/index/genome.fa" --exclude-failed)
        trap "rm -rf $tmp9 $cfg9" EXIT
        echo ""
        echo "--- collect_consensus (exclude_failed_from_consensus) ---"
        snakemake -s "$smk" \
            --configfile "$cfg9" \
            --config "test_outdir=$tmp9" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py collect_consensus "$tmp9" --expect-clusters 6 || rc=1
    fi

    if [[ "$rule_name" == "size_composition" ]]; then
        local tmp5 cfg5
        tmp5=$(mktemp -d)
        cfg5=$(make_cfg "$tmp5" "$tmp5" "$tmp5/index/genome.fa" --srna-sizes 21,22,23,24,25)
        trap "rm -rf $tmp5 $cfg5" EXIT
        echo ""
        echo "--- size_composition (srna sizes 21-25) ---"
        snakemake -s "$smk" \
            --configfile "$cfg5" \
            --config "test_outdir=$tmp5" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py size_composition "$tmp5" --srna-enabled || rc=1

        local tmp8 cfg8
        tmp8=$(mktemp -d)
        cfg8=$(make_cfg "$tmp8" "$tmp8" "$tmp8/index/genome.fa" --max-top-sizes-fraction 0.8)
        trap "rm -rf $tmp8 $cfg8" EXIT
        echo ""
        echo "--- size_composition (top-lengths only, n=2) ---"
        snakemake -s "$smk" \
            --configfile "$cfg8" \
            --config "test_outdir=$tmp8" \
            --cores 1 \
            --quiet 2>&1 || rc=1
        python tests/e2e/unit_assertions.py size_composition "$tmp8" --top-enabled || rc=1
    fi

    return $rc
}

run_unit() {
    for rule in normalize size_composition collect_consensus qc_initial qc_final; do
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

    # The optional GTF. Asserted on HOMER's own report rather than on the flag being
    # present, since -gtf is easy to wire up so that it looks fine and does nothing.
    echo ""
    echo "--- e2e: optional -gtf wiring ---"
    local tmp3 cfg3
    tmp3=$(mktemp -d)
    cfg3=$(make_cfg "$tmp3" "$tmp3" "$tmp3/index/genome.fa" --gtf tests/data/annotation.gtf)
    trap "rm -rf $tmp3 $cfg3" EXIT

    snakemake \
        --cores 4 \
        -s workflow/Snakefile \
        --configfile "$cfg3" \
        2>&1
    python tests/e2e/assertions.py "$tmp3" --skip-golden
    python tests/e2e/check_annotation_wiring.py "$tmp3"
}

# ── per-aligner e2e ───────────────────────────────────────────────────────────
# run_aligner_e2e PROG
# Runs a full pipeline end-to-end with a non-default alignment program and
# asserts all outputs are present (same checks as the default e2e run).
run_aligner_e2e() {
    local aln="$1"
    echo "--- e2e: $aln ---"
    local tmp cfg index_path
    tmp=$(mktemp -d)

    # STAR needs a directory path; other aligners use a file prefix.
    if [[ "$aln" == "STAR" ]]; then
        index_path="$tmp/index_star"
    else
        index_path="$tmp/index/genome"
    fi

    cfg=$(make_cfg "$tmp" "$tmp" "$index_path" --alignment-program "$aln")
    trap "rm -rf $tmp $cfg" EXIT

    snakemake \
        --cores 4 \
        -s workflow/Snakefile \
        --configfile "$cfg" \
        2>&1
    python tests/e2e/assertions.py "$tmp" --skip-golden
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
    star)
        run_step "e2e/STAR" run_aligner_e2e "STAR"
        ;;
    bowtie2)
        run_step "e2e/bowtie2" run_aligner_e2e "bowtie2"
        ;;
    bwa-mem)
        run_step "e2e/bwa-mem" run_aligner_e2e "bwa-mem"
        ;;
    hisat2)
        run_step "e2e/hisat2" run_aligner_e2e "hisat2"
        ;;
    all)
        run_step "schema"      run_schema
        run_step "validate"    run_validate
        run_unit
        run_step "e2e"         run_e2e
        run_step "e2e/STAR"    run_aligner_e2e "STAR"
        run_step "e2e/bowtie2" run_aligner_e2e "bowtie2"
        run_step "e2e/bwa-mem" run_aligner_e2e "bwa-mem"
        run_step "e2e/hisat2"  run_aligner_e2e "hisat2"
        ;;
    *)
        echo "Usage: $0 [schema|validate|unit|e2e|star|bowtie2|bwa-mem|hisat2|all] [-k]" >&2
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
