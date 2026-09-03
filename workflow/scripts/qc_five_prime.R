# 5'-end precision at annotated small-RNA gene classes.
#
# csRNA-seq should place its reads on capped 5' ends. Two annotation classes make
# that testable inside a single library, with no matched input and no
# normalisation, because each is scored as a ratio against its own neighbourhood:
#
#   snRNA  Sm-class spliceosomal snRNAs are Pol II transcripts and are capped, so
#          a working library puts nearly all of their signal on the annotated 5'
#          end. Measured across an enzyme panel this ran 94% with alkaline
#          phosphatase present and 72% without.
#   tRNA   RNase P leaves a 5'-monophosphate exactly at the mature 5' end. That
#          end is ligation-competent the moment the phosphatase step fails, so
#          tRNA precision moves the opposite way, 48% to 62% in the same panel.
#
# snRNA5pPct is the one to gate on. Measured across the enzyme panel it splits
# cleanly: every library with the phosphatase scored 92.4 or better and every
# library without it 88.9 or worse, a gap of 4.5 SD. tRNA5pPct does not separate
# on its own (30.3 to 65.9 with the phosphatase present) because it only moves in
# the most extreme arm, so it is reported as a mechanistic cross-check rather
# than a gate. Do NOT combine them into a ratio or a difference: the tRNA noise
# swamps the snRNA signal and destroys the clean split (tested, both fail).
#
# Reads are taken from the raw bedGraphs, whose coordinates are the 5' base of
# each read, so no re-counting is needed.

suppressWarnings(suppressMessages(suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
})))

WINDOW <- 500L   # neighbourhood that counts as "at this gene"
ONEND  <- 5L     # distance that counts as "on the 5' end"

read_bg <- function(path, strand) {
    if (!nchar(path) || !file.exists(path)) return(GRanges())
    # HOMER writes a `track ...` header line that rtracklayer parses but warns
    # about while coercing; the ranges themselves come through correctly.
    bg <- suppressWarnings(import(path, format = "bedGraph"))
    if (!length(bg)) return(GRanges())
    strand(bg) <- strand
    bg
}

f_sn <- trimws(snakemake@params[["snrnas"]])
f_tr <- trimws(snakemake@params[["trnas"]])
sample_id <- snakemake@wildcards[["sample_cs"]]

out <- data.frame(Sample = sample_id, snRNA5pPct = NA_real_, tRNA5pPct = NA_real_)

if (nchar(f_sn) || nchar(f_tr)) {
    reads <- c(read_bg(trimws(snakemake@input[["bg_p"]]), "+"),
               read_bg(trimws(snakemake@input[["bg_n"]]), "-"))

    # A feature's 5' base, respecting strand; the ratio is then computed against
    # reads on the same strand, so antisense signal cannot inflate it.
    precision <- function(path) {
        if (!nchar(path) || !file.exists(path) || !length(reads)) return(NA_real_)
        feat <- import(path)
        if (!length(feat)) return(NA_real_)
        # a BED without strand would make "5' end" meaningless
        if (all(as.character(strand(feat)) == "*")) {
            warning("no strand in ", path, "; 5'-end precision not computed")
            return(NA_real_)
        }
        five <- resize(feat, 1L, fix = "start")
        d <- distanceToNearest(reads, five, ignore.strand = FALSE)
        if (!length(d)) return(NA_real_)
        dist <- mcols(d)[["distance"]]
        w <- score(reads)[queryHits(d)]
        near <- dist <= WINDOW
        if (!any(near) || sum(w[near]) == 0) return(NA_real_)
        100 * sum(w[near & dist <= ONEND]) / sum(w[near])
    }

    out[["snRNA5pPct"]] <- precision(f_sn)
    out[["tRNA5pPct"]]  <- precision(f_tr)
}

readr::write_tsv(out, snakemake@output[[1]])
