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
