suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

# Read-size composition of every consensus TSS cluster, for ONE csRNA library.
# merge_size_composition.R joins the per-library files into tss.consensus.sizes.txt.
#
# HOMER tag directories store the aligned length of each read next to its 5' position, so
# the size distribution of the reads initiating in a cluster can be measured directly. Two
# things are counted here:
#
#   reads        total reads whose 5' end falls in the cluster
#   srna         of those, how many fall in filtering.tss_srna_sizes
#   top1..topN   of those, how many fall in the cluster's 1..N commonest read lengths
#
# Why the second one: abundant uncapped small RNAs (in plants the 21-25 nt siRNAs, above
# all the 24 nt Pol IV class over transposons) survive the enzymatic depletion well enough
# to be called as TSS clusters, and enrichment over the input library stops rejecting them
# once the csRNA library degrades, because a library that has lost its capped signal is
# proportionally richer in siRNA than its own input. Read size gives them away either way.
# The top-n form needs no list of suspect lengths, so it also catches species outside the
# usual small RNA sizes.
#
# top1..topN are all precomputed so that retuning filtering.tss_top_sizes_n costs nothing:
# only changing filtering.tss_srna_sizes or tss_top_sizes_max forces a re-read of the tag
# directories, which is the whole cost of this rule. Keeping the full per-cluster length
# histogram instead would serve any n exactly, but it runs to tens of millions of rows on a
# plant-sized dataset (about 2.8M per library over 113k clusters) against a handful of
# columns here, so the cumulative sums are kept instead.
#
# Clusters never overlap within a strand (collect_consensus_tss.R resolves every overlap
# before writing tss.consensus.bed), so a read's 5' position identifies its cluster
# unambiguously and findInterval() is exact.

TAG_COLS <- c("name", "chr", "pos", "strand", "count", "len")
TAG_TYPES <- "cciidi"

split_param <- function(x) {
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(character(0))
    strsplit(x, "[[:space:]]+", perl = TRUE)[[1]]
}

srna_sizes <- suppressWarnings(as.integer(split_param(snakemake@params[["srna_sizes"]])))
srna_sizes <- sort(unique(srna_sizes[!is.na(srna_sizes)]))

top_max <- suppressWarnings(as.integer(snakemake@params[["top_sizes_max"]]))
if (is.na(top_max) || top_max < 1L) top_max <- 5L
max_top_frac <- suppressWarnings(as.numeric(snakemake@params[["max_top_sizes_fraction"]]))
if (is.na(max_top_frac)) max_top_frac <- 1
want_top <- max_top_frac < 1

tss <- import(snakemake@input[["bed"]])
tss_names <- mcols(tss)[["name"]]
n_tss <- length(tss_names)

# Both filters off: write the cluster list only and leave the tag directory alone, since
# reading it is the whole cost of this rule.
if (!length(srna_sizes) && !want_top) {
    cat("no size filter configured: no tag directory read\n")
    readr::write_tsv(data.frame(TSS = tss_names), snakemake@output[[1]])
    quit(save = "no", status = 0)
}

tagdir <- as.character(snakemake@input[["td"]])
if (length(srna_sizes)) {
    cat("Small-RNA size classes: ", paste(srna_sizes, collapse = ", "), " nt\n", sep = "")
}
if (want_top) cat("Precomputing the top 1..", top_max, " read-length shares\n", sep = "")

tss_chr <- as.character(seqnames(tss))
tss_str <- as.character(strand(tss))
chroms  <- unique(tss_chr)
strands <- c("+", "-")

# Per chromosome and strand: cluster bounds in ascending start order, with the row each
# cluster occupies in the output vectors.
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

# Accumulate weights w into vec at the (1-based) positions given by row.
accumulate <- function(vec, row, w) {
    agg <- rowsum(w, row)
    at <- as.numeric(rownames(agg))
    vec[at] <- vec[at] + agg[, 1L]
    vec
}

reads <- numeric(n_tss)
srna  <- numeric(n_tss)
len_keys <- vector("list", length(chroms) * length(strands))
len_vals <- len_keys
nchunk <- 0L

for (ch in chroms) {
    f <- file.path(tagdir, paste0(ch, ".tags.tsv"))
    if (!file.exists(f)) next
    first <- readLines(f, n = 1L)                    # a chromosome with no reads can
    if (!length(first) || !nzchar(first)) next       # leave an empty tag file
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
        sel <- which(tags$strand == (st - 1L))       # HOMER strand is 0 (+) / 1 (-)
        if (!length(sel)) next
        j <- findInterval(tags$pos[sel], k$start)
        hit <- j > 0L
        hit[hit] <- tags$pos[sel][hit] <= k$end[j[hit]]
        if (!any(hit)) next
        row <- k$row[j[hit]]
        sel <- sel[hit]
        reads <- accumulate(reads, row, tags$count[sel])
        si <- is_srna[sel]
        if (any(si)) srna <- accumulate(srna, row[si], tags$count[sel][si])
        if (want_top) {
            # Collapse to per-(cluster, length) totals now, so what is carried to the end
            # stays small. The key is built in double arithmetic because cluster count
            # times read length overflows the integer range on a large genome with long
            # reads, where integer multiplication would silently give NA.
            agg <- rowsum(tags$count[sel],
                row + (as.numeric(tags$len[sel]) - 1) * n_tss)
            nchunk <- nchunk + 1L
            len_keys[[nchunk]] <- as.numeric(rownames(agg))
            len_vals[[nchunk]] <- agg[, 1L]
        }
    }
}

out <- data.frame(TSS = tss_names, reads = reads, stringsAsFactors = FALSE)
if (length(srna_sizes)) out[["srna"]] <- srna

if (want_top) {
    tops <- matrix(0, nrow = n_tss, ncol = top_max,
        dimnames = list(NULL, paste0("top", seq_len(top_max))))
    if (nchunk) {
        agg <- rowsum(unlist(len_vals[seq_len(nchunk)]),
                      unlist(len_keys[seq_len(nchunk)]))
        k <- as.numeric(rownames(agg)); v <- agg[, 1L]
        row <- ((k - 1) %% n_tss) + 1
        o <- order(row, -v)
        row <- row[o]; v <- v[o]
        # Rank of each length within its cluster, then a cumulative sum per rank cut-off.
        rank_in_row <- sequence(rle(row)$lengths)
        for (n in seq_len(top_max)) {
            keep <- rank_in_row <= n
            a <- rowsum(v[keep], row[keep])
            tops[as.numeric(rownames(a)), n] <- a[, 1L]
        }
    }
    out <- cbind(out, as.data.frame(tops))
}

readr::write_tsv(out, snakemake@output[[1]])
cat("  ", basename(tagdir), ": ", format(sum(reads), big.mark = ","),
    " reads in clusters\n", sep = "")
