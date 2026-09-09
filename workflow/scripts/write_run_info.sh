#!/usr/bin/env bash
# Record what produced this results directory: pipeline version and commit, the resolved
# config and sample table, and the version of every tool the run depends on.
#
# Written because results outlive memory. Reopening an old results folder otherwise leaves
# no way to tell which annotation, which HOMER build or which pipeline commit made it, and
# the answer can change the interpretation.
set -uo pipefail

# Never fail the run over provenance collection: a missing tool is recorded as such.
ver() {
    local name="$1"; shift
    if ! command -v "$1" > /dev/null 2>&1; then
        printf '  %-12s (not on PATH)\n' "$name"
        return 0
    fi
    local out
    out=$("$@" 2>&1 | grep -aoEm1 '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
    printf '  %-12s %s\n' "$name" "${out:-unknown}"
}

{
    printf 'csRNAseq_pipeline run\n'
    printf 'finished        %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'pipeline        %s\n' "${PIPELINE_VERSION:-unknown}"
    if git -C "$REPO_DIR" rev-parse --git-dir > /dev/null 2>&1; then
        printf 'commit          %s%s\n' \
            "$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null)" \
            "$(git -C "$REPO_DIR" diff --quiet HEAD 2>/dev/null || echo ' (working tree dirty)')"
    else
        printf 'commit          (not a git checkout)\n'
    fi
    printf 'host            %s\n' "$(hostname)"
    printf 'snakemake       %s\n' "${SNAKEMAKE_VERSION:-unknown}"

    printf '\nsample table    %s\n' "$SAMPLE_TABLE"
    if [[ -f "$SAMPLE_TABLE" ]]; then
        printf 'sample table sha256  %s\n' \
            "$( (shasum -a 256 "$SAMPLE_TABLE" 2>/dev/null || sha256sum "$SAMPLE_TABLE" 2>/dev/null) | awk '{print $1}')"
    fi

    printf '\ntool versions\n'
    ver snakemake snakemake --version
    ver R R --version
    ver samtools samtools --version
    ver bedtools bedtools --version
    ver bwa bwa
    ver STAR STAR --version
    ver bowtie2 bowtie2 --version
    ver hisat2 hisat2 --version
    ver bfqutils bfqutils --version
    ver bam2td bam2td --version
    # HOMER help text contains numeric thresholds that are not software versions.
    # Its config.txt and executable fingerprints below identify the installed build.
    printf '  %-12s %s\n' HOMER "$(command -v findcsRNATSS.pl > /dev/null 2>&1 && dirname "$(command -v findcsRNATSS.pl)" || echo '(not on PATH)')"

    printf '\ncontent fingerprints and installed packages\n'
    python "$REPO_DIR/workflow/scripts/provenance.py" "$CONFIG_DUMP" "$REPO_DIR" || exit 1

    printf '\nresolved config\n'
    sed 's/^/  /' "$CONFIG_DUMP"
} > "$OUT"
