suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

# Read-size composition of every consensus TSS cluster.
#
# HOMER tag directories store the aligned length of each read next to its 5'
# position, so the size distribution of the reads initiating in a cluster can be
# measured directly.  For each csRNA library this counts the reads whose 5' end
# falls in each cluster, and how many of those fall in the configured small-RNA
# size classes (filtering.tss_srna_sizes).  normalize_tss_quantification.R turns
# those two numbers into the size-composition filter.
#
# Why the measurement is worth making: abundant uncapped small RNAs (in plants the
# 21-25 nt siRNAs, above all the 24 nt Pol IV class over transposons) survive the
# enzymatic depletion well enough to be called as TSS clusters, and enrichment over
# the input library cannot reject them when the input is shallow.  Their read size
# gives them away, because a genuine capped RNA population is not concentrated in a
# single small size class.
#
# Clusters never overlap within a strand (collect_consensus_tss.R resolves every
# overlap before writing tss.consensus.bed), so a read's 5' position identifies its
# cluster unambiguously and findInterval() is exact.

TAG_COLS <- c("name", "chr", "pos", "strand", "count", "len")
TAG_TYPES <- "cciidi"

split_param <- function(x) {
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(character(0))
    strsplit(x, "[[:space:]]+", perl = TRUE)[[1]]
}

# Accumulate weights w into vec at the (1-based) positions given by row.
accumulate <- function(vec, row, w) {
    agg <- rowsum(w, row)
    vec[as.integer(rownames(agg))] <- vec[as.integer(rownames(agg))] + agg[, 1L]
    vec
}

srna_sizes <- suppressWarnings(as.integer(split_param(snakemake@params[["srna_sizes"]])))
srna_sizes <- sort(unique(srna_sizes[!is.na(srna_sizes)]))

tss <- import(snakemake@input[["bed"]])
tss_names <- mcols(tss)[["name"]]

# Filter disabled: write the cluster list only and leave the tag directories alone,
# since reading them is the whole cost of this rule.
if (!length(srna_sizes)) {
    cat("tss_srna_sizes is empty: size-composition filter disabled, no tag directories read\n")
    readr::write_tsv(data.frame(TSS = tss_names), snakemake@output[[1]])
    quit(save = "no", status = 0)
}

tagdirs <- as.character(snakemake@input[["td"]])
# HOMER names its quantification columns after the tag directory, so the same
# basename keys this table to the columns of tss.consensus.homer.raw.txt.
samples <- basename(tagdirs)

cat("Small-RNA size classes: ", paste(srna_sizes, collapse = ", "), " nt\n", sep = "")
cat("Measuring ", length(tss_names), " clusters across ", length(samples),
    " csRNA libraries\n", sep = "")

tss_chr <- as.character(seqnames(tss))
tss_str <- as.character(strand(tss))
chroms  <- unique(tss_chr)
strands <- c("+", "-")

# Per chromosome and strand: cluster bounds in ascending start order, with the row
# each cluster occupies in the output matrices.
idx <- list()
for (ch in chroms) {
    for (st in strands) {
        i <- which(tss_chr == ch & tss_str == st)
        if (!length(i)) next
        i <- i[order(start(tss)[i])]
        idx[[paste0(ch, "\t", st)]] <- list(
            row = i, start = start(tss)[i], end = end(tss)[i])
    }
}

reads <- matrix(0, nrow = length(tss_names), ncol = length(samples),
    dimnames = list(tss_names, samples))
srna <- reads

for (s in seq_along(tagdirs)) {
    for (ch in chroms) {
        f <- file.path(tagdirs[s], paste0(ch, ".tags.tsv"))
        if (!file.exists(f)) next
        first <- readLines(f, n = 1L)                    # a chromosome with no reads
        if (!length(first) || !nzchar(first)) next       # can leave an empty tag file
        ncol_f <- length(strsplit(first, "\t", fixed = TRUE)[[1]])
        if (ncol_f != length(TAG_COLS)) {
            stop("unexpected tag file layout in ", f, ": ", ncol_f, " columns, expected ",
                 length(TAG_COLS), " (name, chr, pos, strand, count, len)")
        }
        tags <- suppressWarnings(readr::read_tsv(f, col_names = TAG_COLS,
            col_types = TAG_TYPES, progress = FALSE))
        is_srna <- tags$len %in% srna_sizes
        for (st in seq_along(strands)) {
            k <- idx[[paste0(ch, "\t", strands[st])]]
            if (is.null(k)) next
            sel <- which(tags$strand == (st - 1L))   # HOMER strand is 0 (+) / 1 (-)
            if (!length(sel)) next
            j <- findInterval(tags$pos[sel], k$start)
            hit <- j > 0L
            hit[hit] <- tags$pos[sel][hit] <= k$end[j[hit]]
            if (!any(hit)) next
            row <- k$row[j[hit]]
            sel <- sel[hit]
            reads[, s] <- accumulate(reads[, s], row, tags$count[sel])
            si <- is_srna[sel]
            if (any(si)) srna[, s] <- accumulate(srna[, s], row[si], tags$count[sel][si])
        }
    }
    cat("  ", samples[s], ": ", format(sum(reads[, s]), big.mark = ","),
        " reads in clusters, ", round(100 * sum(srna[, s]) / max(sum(reads[, s]), 1), 1),
        "% in the small-RNA sizes\n", sep = "")
}

out <- data.frame(TSS = tss_names, stringsAsFactors = FALSE)
for (s in seq_along(samples)) {
    out[[paste0(samples[s], ".reads")]] <- reads[, s]
    out[[paste0(samples[s], ".srna")]]  <- srna[, s]
}
readr::write_tsv(out, snakemake@output[[1]])
