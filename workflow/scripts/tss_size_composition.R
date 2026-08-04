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
# enzymatic depletion well enough to be called as TSS clusters, and enrichment over the
# input library stops rejecting them once the csRNA library degrades, because a library
# that has lost its capped signal is proportionally richer in siRNA than its own input.
# Read size gives them away either way, because a genuine capped RNA population is not
# concentrated in a single small size class.
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

# Second, independent measure: how much of a cluster sits in its top few read lengths,
# whatever those lengths are. Real initiation is heterogeneous (a promoter yields reads
# across tens of lengths), while a discretely processed RNA puts nearly everything into
# one or two. That catches contaminants the fixed size list cannot, e.g. the 27-28 nt
# 5'-polyphosphate species over some transposons.
top_n <- suppressWarnings(as.integer(snakemake@params[["top_sizes_n"]]))
if (is.na(top_n) || top_n < 1L) top_n <- 2L
max_top_frac <- suppressWarnings(as.numeric(snakemake@params[["max_top_sizes_fraction"]]))
if (is.na(max_top_frac)) max_top_frac <- 1
want_top <- max_top_frac < 1

tss <- import(snakemake@input[["bed"]])
tss_names <- mcols(tss)[["name"]]
n_tss <- length(tss_names)

# Both filters off: write the cluster list only and leave the tag directories alone,
# since reading them is the whole cost of this rule.
if (!length(srna_sizes) && !want_top) {
    cat("no size filter configured (tss_srna_sizes empty, tss_max_top_sizes_fraction 1):",
        "no tag directories read\n")
    readr::write_tsv(data.frame(TSS = tss_names), snakemake@output[[1]])
    quit(save = "no", status = 0)
}

tagdirs <- as.character(snakemake@input[["td"]])
# HOMER names its quantification columns after the tag directory, so the same
# basename keys this table to the columns of tss.consensus.homer.raw.txt.
samples <- basename(tagdirs)

if (length(srna_sizes)) {
    cat("Small-RNA size classes: ", paste(srna_sizes, collapse = ", "), " nt\n", sep = "")
}
if (want_top) cat("Also measuring the share held by the top ", top_n, " read lengths\n", sep = "")
cat("Measuring ", n_tss, " clusters across ", length(samples),
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

reads <- matrix(0, nrow = n_tss, ncol = length(samples),
    dimnames = list(tss_names, samples))
srna <- reads
topn <- reads

# Reads in a cluster's top `top_n` lengths, from per-(cluster, length) totals keyed as
# row + (len - 1) * n_tss. Ranking within a cluster is a sort plus a within-group index,
# so no per-cluster loop is needed. The key is built in double rather than integer
# arithmetic: cluster count times read length exceeds the integer range for a large
# genome sequenced with long reads (a million clusters and multi-kilobase reads), where
# integer multiplication would silently return NA.
top_share <- function(key, val) {
    agg <- rowsum(val, key)
    k <- as.numeric(rownames(agg)); v <- agg[, 1L]
    row <- ((k - 1) %% n_tss) + 1
    o <- order(row, -v)
    row <- row[o]; v <- v[o]
    rank_in_row <- sequence(rle(row)$lengths)
    keep <- rank_in_row <= top_n
    out <- numeric(n_tss)
    agg2 <- rowsum(v[keep], row[keep])
    out[as.numeric(rownames(agg2))] <- agg2[, 1L]
    out
}

for (s in seq_along(tagdirs)) {
    len_keys <- if (want_top) vector("list", length(chroms) * length(strands)) else NULL
    len_vals <- len_keys
    nchunk <- 0L
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
            if (want_top) {
                # Collapse to per-(cluster, length) totals now, so what is carried to the
                # end of the sample stays small.
                agg <- rowsum(tags$count[sel],
                    row + (as.numeric(tags$len[sel]) - 1) * n_tss)
                nchunk <- nchunk + 1L
                len_keys[[nchunk]] <- as.numeric(rownames(agg))
                len_vals[[nchunk]] <- agg[, 1L]
            }
        }
    }
    if (want_top && nchunk) {
        topn[, s] <- top_share(unlist(len_keys[seq_len(nchunk)]),
            unlist(len_vals[seq_len(nchunk)]))
        len_keys <- len_vals <- NULL
    }
    cat("  ", samples[s], ": ", format(sum(reads[, s]), big.mark = ","),
        " reads in clusters", sep = "")
    if (length(srna_sizes)) cat(", ",
        round(100 * sum(srna[, s]) / max(sum(reads[, s]), 1), 1), "% in the small-RNA sizes",
        sep = "")
    cat("\n")
}

out <- data.frame(TSS = tss_names, stringsAsFactors = FALSE)
for (s in seq_along(samples)) {
    out[[paste0(samples[s], ".reads")]] <- reads[, s]
    if (length(srna_sizes)) out[[paste0(samples[s], ".srna")]] <- srna[, s]
    if (want_top)           out[[paste0(samples[s], ".topn")]] <- topn[, s]
}
readr::write_tsv(out, snakemake@output[[1]])
