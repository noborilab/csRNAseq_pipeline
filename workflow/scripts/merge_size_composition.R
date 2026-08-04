suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

# Join the per-library read-size composition tables into tss.consensus.sizes.txt.
#
# Each input holds one row per consensus TSS with columns reads / srna / top1..topN
# (whichever the configured filters need). Here they are prefixed with the library name so
# that normalize_tss_quantification.R can key them to the columns of
# tss.consensus.homer.raw.txt, which HOMER names after the tag directory.

tss <- import(snakemake@input[["bed"]])
tss_names <- mcols(tss)[["name"]]

files <- as.character(snakemake@input[["sizes"]])
# results/tss_sizes/<library>.sizes.txt
samples <- sub("\\.sizes\\.txt$", "", basename(files))

out <- data.frame(TSS = tss_names, stringsAsFactors = FALSE)
for (i in seq_along(files)) {
    d <- as.data.frame(suppressWarnings(suppressMessages(
        readr::read_tsv(trimws(files[i]), progress = FALSE))))
    if (!identical(as.character(d[["TSS"]]), as.character(tss_names))) {
        stop("cluster order in ", files[i], " does not match tss.consensus.bed; ",
             "delete results/tss_sizes and re-run the tss_size_composition rules")
    }
    cols <- setdiff(colnames(d), "TSS")
    for (cl in cols) out[[paste0(samples[i], ".", cl)]] <- d[[cl]]
}

readr::write_tsv(out, snakemake@output[[1]])
cat("Collected size composition for ", length(files), " librar",
    if (length(files) == 1) "y" else "ies", " over ", length(tss_names),
    " clusters", if (ncol(out) == 1) " (no size filter configured)" else "", "\n", sep = "")
