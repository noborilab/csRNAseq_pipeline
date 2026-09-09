suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(edgeR))))
suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

read_homer_quant <- function(path) {
    df  <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(path))))
    ids <- df[[1]]
    idx <- grep(" Tag Count", colnames(df))
    m   <- as.matrix(df[, idx, drop = FALSE])
    colnames(m) <- basename(gsub(" Tag Count .+", "", colnames(m)))
    rownames(m) <- ids
    m
}

tss <- import(snakemake@input[["bed"]])   # tss.consensus.bed

quant <- as.data.frame(suppressWarnings(suppressMessages(readr::read_tsv(trimws(snakemake@input[["quant"]])))))
tag_cols <- grep(" Tag Count", colnames(quant))
quant_m <- as.matrix(quant[, tag_cols, drop = FALSE])
colnames(quant_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_m)))
rownames(quant_m) <- quant[[1]]

# CPM filter on library-size CPM (before TMM, per edgeR best practice)
min_cpm     <- snakemake@params[["min_cpm"]]
min_samples <- snakemake@params[["min_samples"]]
libsizes <- colSums(quant_m)
cpm_lib  <- t(t(quant_m) / libsizes) * 1e6
keep     <- rowSums(cpm_lib >= min_cpm) >= min_samples
quant_m  <- quant_m[keep, , drop = FALSE]
tss      <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
cat("CPM filter: kept ", sum(keep), " of ", length(keep), " consensus TSSs\n", sep = "")

# Read-size composition filters.  Both read tss.consensus.sizes.txt from the
# tss_size_composition rule, so neither needs an extra HOMER run, and both work the same
# way: a cluster is dropped when some measure of its reads exceeds a fraction, in at least
# min_samples libraries that have at least min_reads reads there.
#
# Known-size and top-length concentration screens are empirical indicators of
# processed-RNA contamination. Size selection, trimming and sampling can also narrow
# genuine initiation products; neither screen determines RNA origin or cap state.
srna_sizes       <- trimws(as.character(snakemake@params[["srna_sizes"]]))
srna_sizes       <- if (nzchar(srna_sizes)) strsplit(srna_sizes, "[[:space:]]+")[[1]] else character(0)
max_srna_frac    <- snakemake@params[["max_srna_fraction"]]
srna_min_reads   <- snakemake@params[["srna_min_reads"]]
srna_min_samples <- snakemake@params[["srna_min_samples"]]
top_n            <- snakemake@params[["top_sizes_n"]]
max_top_frac     <- snakemake@params[["max_top_sizes_fraction"]]

# TSS names to keep: those where fewer than min_samples libraries have more than max_frac
# of their reads in the `suffix` measure.  Only libraries with enough reads in a cluster
# get a vote, so a couple of stray reads cannot condemn it.
size_filter <- function(comp, suffix, max_frac, min_reads, min_samples, label, why) {
    read_cols <- grep("\\.reads$", colnames(comp), value = TRUE)
    num_cols  <- sub("\\.reads$", suffix, read_cols)
    if (!length(read_cols) || !all(num_cols %in% colnames(comp))) {
        stop("tss.consensus.sizes.txt lacks the ", suffix, " columns this filter needs (",
             why, " is set). Re-run the size composition rules: delete ",
             "tss.consensus.sizes.txt and results/tss_sizes. If the missing column is a ",
             "top-n one, check that filtering.tss_top_sizes_n is no greater than ",
             "filtering.tss_top_sizes_max.")
    }
    ids <- intersect(rownames(quant_m), rownames(comp))
    r <- as.matrix(comp[ids, read_cols, drop = FALSE])
    v <- as.matrix(comp[ids, num_cols, drop = FALSE])
    dominated <- (r >= min_reads) & (v / pmax(r, 1)) > max_frac
    keep <- rowSums(dominated) < min_samples
    cat(label, ": kept ", sum(keep), " of ", length(keep), " TSSs\n", sep = "")
    ids[keep]
}

if (length(srna_sizes) > 0 || max_top_frac < 1) {
    comp <- as.data.frame(suppressWarnings(suppressMessages(
        readr::read_tsv(trimws(snakemake@input[["sizes"]]), progress = FALSE))))
    rownames(comp) <- comp[["TSS"]]

    if (length(srna_sizes) > 0) {
        keep_ids <- size_filter(comp, ".srna", max_srna_frac, srna_min_reads,
            srna_min_samples, sprintf(
                "Small-RNA size filter (>%g%% of reads at %s nt in >=%d librar%s with >=%d reads)",
                100 * max_srna_frac, paste(srna_sizes, collapse = "/"), srna_min_samples,
                if (srna_min_samples == 1) "y" else "ies", srna_min_reads),
            "tss_srna_sizes")
        quant_m <- quant_m[keep_ids, , drop = FALSE]
        tss     <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    }

    if (max_top_frac < 1) {
        # top1..topN are all precomputed by tss_size_composition, so tss_top_sizes_n can be
        # retuned without re-reading the tag directories; pick the requested one here.
        keep_ids <- size_filter(comp, paste0(".top", top_n), max_top_frac, srna_min_reads,
            srna_min_samples, sprintf(
                "Top-%d-lengths filter (>%g%% of reads in the %d commonest lengths in >=%d librar%s with >=%d reads)",
                top_n, 100 * max_top_frac, top_n, srna_min_samples,
                if (srna_min_samples == 1) "y" else "ies", srna_min_reads),
            "tss_max_top_sizes_fraction")
        quant_m <- quant_m[keep_ids, , drop = FALSE]
        tss     <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    }
}

# Optional enrichment filter: both counts must describe exactly the same consensus
# regions. Normalize against all retained non-organelle tags, not counts in the
# selected TSS set, so changing that set does not change the depth denominator.
min_cs_in_ratio   <- snakemake@params[["min_cs_in_ratio"]]
min_ratio_samples <- snakemake@params[["min_ratio_samples"]]

if (min_cs_in_ratio > 0) {
    input_table <- as.data.frame(readr::read_tsv(snakemake@input[["quant_in"]],
                                               show_col_types = FALSE))
    quant_in_m <- read_homer_quant(snakemake@input[["quant_in"]])
    if (anyDuplicated(quant[[1]]) || anyDuplicated(input_table[[1]]) ||
        !setequal(quant[[1]], input_table[[1]])) {
        stop("csRNA and input quantification must contain the same unique consensus TSS IDs")
    }
    input_order <- match(quant[[1]], input_table[[1]])
    # HOMER tables share a coordinate convention; compare their location columns
    # directly rather than assuming that sequential TSS IDs identify the same loci.
    if (!identical(unname(as.matrix(quant[, 2:5])),
                   unname(as.matrix(input_table[input_order, 2:5])))) {
        stop("csRNA and input consensus TSS coordinates/strands disagree")
    }
    cs_ids <- strsplit(snakemake@params[["cs_ids"]], " ", fixed = TRUE)[[1]]
    in_ids <- strsplit(snakemake@params[["paired_in_ids"]], " ", fixed = TRUE)[[1]]
    if (length(cs_ids) != length(in_ids) || !all(cs_ids %in% colnames(quant_m)) ||
        !all(in_ids %in% colnames(quant_in_m))) stop("invalid csRNA/input sample pairing")

    nuclear_depth <- function(path, ids) {
        stats <- read.delim(path, check.names = FALSE)
        if (anyDuplicated(stats$Sample) || !all(ids %in% stats$Sample))
            stop("missing or duplicate library depth records")
        depth <- (stats$TotalReads - stats$OrganelleReads)[match(ids, stats$Sample)]
        if (length(depth) != length(ids) || any(!is.finite(depth) | depth <= 0))
            stop("ratio filter requires positive finite non-organelle library depths")
        depth
    }
    cs_depth <- nuclear_depth(snakemake@input[["stats_cs"]], cs_ids)
    in_depth <- nuclear_depth(snakemake@input[["stats_in"]], in_ids)
    # One input read is the explicit denominator floor. This avoids infinite
    # ratios at zero input; it is a screening heuristic, not a significance test.
    ratio_m <- matrix(0, nrow(quant_m), length(cs_ids))
    for (i in seq_along(cs_ids)) {
        ratio_m[, i] <- (quant_m[, cs_ids[i]] / cs_depth[i]) /
            (pmax(quant_in_m[rownames(quant_m), in_ids[i]], 1) / in_depth[i])
    }
    keep_ratio <- rowSums(ratio_m >= min_cs_in_ratio) >= min_ratio_samples
    quant_m <- quant_m[keep_ratio, , drop = FALSE]
    tss <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    cat("Depth-normalized ratio filter: kept ", sum(keep_ratio), " of ",
        length(keep_ratio), " TSSs\n", sep = "")
}

if (!nrow(quant_m)) {
    stop("no TSSs left after filtering. Loosen filtering.tss_min_cpm / tss_min_samples, ",
         "tss_min_cs_in_ratio, or the size-composition filter (tss_srna_sizes / ",
         "tss_max_srna_fraction).")
}

# Row order. annotatePeaks.pl emits the clusters in a different order on every run and
# nothing above imposes one, so two runs on the same data produced count tables holding
# the same numbers with their rows shuffled, while tss.final.bed agreed with neither
# because it keeps the consensus BED's order all the way through. Ordering the counts by
# that same BED, which is coordinate-sorted and reproducible, settles all three at once:
# the tables come out byte-identical between runs, they read in genomic order beside the
# BED, and MinReplCorrPearson stops drifting in the last bit of a double, which is what
# cor() did when it added the same pairs up in a new order each run.
ord <- match(mcols(tss)[["name"]], rownames(quant_m))
if (anyNA(ord) || length(ord) != nrow(quant_m)) {
    stop("the consensus BED and the count matrix disagree after filtering: ",
         length(ord), " clusters in the BED against ", nrow(quant_m), " rows of counts. ",
         "That is a bug in the filtering above rather than anything in the config.")
}
quant_m <- quant_m[ord, , drop = FALSE]

export.bed(tss, snakemake@output[["bed"]])

readr::write_tsv(as.data.frame(cbind(TSS = rownames(quant_m), quant_m)),
    snakemake@output[["raw"]])

y <- DGEList(quant_m, group = gsub('_csrna[0-9]+$', '', colnames(quant_m)))
y <- calcNormFactors(y, method = 'TMMwsp')

nf <- y[['samples']]
nf[['sample']] <- rownames(nf)
nf[['final.factors']] <- nf[['norm.factors']] * nf[['lib.size']]
nf[['mult.per.million']] <- (nf[['final.factors']] / 1000000)^-1
readr::write_tsv(nf, snakemake@output[["normfactors"]])

y_cpm <- cpm(y)
readr::write_tsv(as.data.frame(cbind(TSS = rownames(y_cpm), y_cpm)),
    snakemake@output[["norm"]])
