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
quant_m <- as.matrix(quant[, tag_cols])
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
#   1. Small-RNA sizes (filtering.tss_srna_sizes, in plants the 21-25 nt siRNAs).  Those
#      species are uncapped but survive the enzymatic depletion well enough to be called
#      as TSS clusters, and enrichment over the input cannot reject them where the input
#      has too little coverage to measure a background.
#   2. Top-n lengths (filtering.tss_max_top_sizes_fraction), whatever those lengths are.
#      Genuine initiation is heterogeneous: a promoter yields reads across tens of
#      lengths, so its two commonest lengths hold only a fifth or so of it.  A discretely
#      processed RNA puts nearly everything into one or two lengths.  This catches
#      contaminants that fall outside any fixed size list, e.g. the 27-28 nt
#      5'-polyphosphate species found over some transposons.
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
             why, " is set). Re-run the tss_size_composition rule ",
             "(delete tss.consensus.sizes.txt).")
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
        keep_ids <- size_filter(comp, ".topn", max_top_frac, srna_min_reads,
            srna_min_samples, sprintf(
                "Top-%d-lengths filter (>%g%% of reads in the %d commonest lengths in >=%d librar%s with >=%d reads)",
                top_n, 100 * max_top_frac, top_n, srna_min_samples,
                if (srna_min_samples == 1) "y" else "ies", srna_min_reads),
            "tss_max_top_sizes_fraction")
        quant_m <- quant_m[keep_ids, , drop = FALSE]
        tss     <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    }
}

# csRNA/sRNA ratio filter: keep TSSs where (csRNA reads / paired input reads)
# >= tss_min_cs_in_ratio in >= tss_min_ratio_samples samples.
# Uses input reads at the initial csRNA TSS set (all_cs.tss_merged_quant_in.txt),
# subsetted to the consensus TSS IDs that survived the CPM filter.
min_cs_in_ratio  <- snakemake@params[["min_cs_in_ratio"]]
min_ratio_samples <- snakemake@params[["min_ratio_samples"]]

if (min_cs_in_ratio > 0) {
    quant_in_m <- read_homer_quant(snakemake@input[["quant_in"]])
    common_ids <- intersect(rownames(quant_m), rownames(quant_in_m))
    quant_in_m <- quant_in_m[common_ids, , drop = FALSE]
    quant_ratio_m <- quant_m[common_ids, , drop = FALSE]

    cs_ids <- strsplit(snakemake@params[["cs_ids"]], " ", fixed = TRUE)[[1]]
    in_ids <- strsplit(snakemake@params[["paired_in_ids"]], " ", fixed = TRUE)[[1]]

    ratio_m <- sapply(seq_along(cs_ids), function(i) {
        quant_ratio_m[, cs_ids[i]] / pmax(quant_in_m[, in_ids[i]], 1)
    })
    if (is.null(dim(ratio_m))) dim(ratio_m) <- c(length(ratio_m), 1)
    rownames(ratio_m) <- common_ids

    keep_ratio <- rowSums(ratio_m >= min_cs_in_ratio) >= min_ratio_samples
    quant_m <- quant_m[common_ids, , drop = FALSE][keep_ratio, , drop = FALSE]
    tss     <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    cat("Ratio filter: kept ", sum(keep_ratio), " of ", length(keep_ratio), " TSSs\n", sep = "")
}

if (!nrow(quant_m)) {
    stop("no TSSs left after filtering. Loosen filtering.tss_min_cpm / tss_min_samples, ",
         "tss_min_cs_in_ratio, or the size-composition filter (tss_srna_sizes / ",
         "tss_max_srna_fraction).")
}

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
