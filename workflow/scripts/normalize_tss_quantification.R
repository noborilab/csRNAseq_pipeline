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

# Read-size composition filter: drop TSSs whose reads are dominated by the small
# RNA size classes listed in filtering.tss_srna_sizes (in plants, the 21-25 nt
# siRNAs).  Those species are uncapped, but they survive the enzymatic depletion
# well enough to be called as TSS clusters, and enrichment over the input cannot
# reject them where the input has too little coverage to measure a background.
# Read length separates them from genuine initiation, which is not confined to one
# small size class.  Counts come from tss.consensus.sizes.txt (tss_size_composition
# rule), so no extra HOMER run is needed.
srna_sizes       <- trimws(as.character(snakemake@params[["srna_sizes"]]))
srna_sizes       <- if (nzchar(srna_sizes)) strsplit(srna_sizes, "[[:space:]]+")[[1]] else character(0)
max_srna_frac    <- snakemake@params[["max_srna_fraction"]]
srna_min_reads   <- snakemake@params[["srna_min_reads"]]
srna_min_samples <- snakemake@params[["srna_min_samples"]]

if (length(srna_sizes) > 0) {
    comp <- as.data.frame(suppressWarnings(suppressMessages(
        readr::read_tsv(trimws(snakemake@input[["sizes"]]), progress = FALSE))))
    rownames(comp) <- comp[["TSS"]]
    read_cols <- grep("\\.reads$", colnames(comp), value = TRUE)
    srna_cols <- sub("\\.reads$", ".srna", read_cols)
    if (!length(read_cols) || !all(srna_cols %in% colnames(comp))) {
        stop("tss.consensus.sizes.txt has no per-sample columns, but tss_srna_sizes is set. ",
             "Re-run the tss_size_composition rule (delete tss.consensus.sizes.txt).")
    }
    ids <- intersect(rownames(quant_m), rownames(comp))
    r <- as.matrix(comp[ids, read_cols, drop = FALSE])
    v <- as.matrix(comp[ids, srna_cols, drop = FALSE])
    # Only libraries with enough reads in a cluster get a vote on it, so a couple of
    # stray small-RNA reads cannot condemn a cluster.
    testable <- r >= srna_min_reads
    dominated <- testable & (v / pmax(r, 1)) > max_srna_frac
    n_over <- rowSums(dominated)
    keep_size <- n_over < srna_min_samples
    quant_m <- quant_m[ids, , drop = FALSE][keep_size, , drop = FALSE]
    tss     <- tss[mcols(tss)[["name"]] %in% rownames(quant_m)]
    cat("Size-composition filter (>", 100 * max_srna_frac, "% of reads at ",
        paste(srna_sizes, collapse = "/"), " nt in >=", srna_min_samples,
        " librar", if (srna_min_samples == 1) "y" else "ies", " with >=",
        srna_min_reads, " reads): kept ", sum(keep_size), " of ", length(keep_size),
        " TSSs\n", sep = "")
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
