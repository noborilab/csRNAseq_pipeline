# Regression tests for the production consensus placement, including chromosome
# boundaries and chains that pairwise shifting cannot reliably resolve.
suppressPackageStartupMessages(library(rtracklayer))
expressions <- parse("workflow/scripts/collect_consensus_tss.R")
definition <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("resolve_tss_overlaps")), as.list(expressions))
stopifnot(length(definition) == 1L)
eval(definition[[1]])
place <- function(x, lengths) {
    invisible(capture.output(y <- resolve_tss_overlaps(x, lengths)))
    y
}
check <- function(before, after, lengths) {
    stopifnot(length(before) == length(after),
        identical(sort(width(before)), sort(width(after))),
        all(start(after) >= 1),
        all(end(after) <= lengths[as.character(seqnames(after))]))
    hits <- findOverlaps(after, ignore.strand = FALSE)
    stopifnot(all(queryHits(hits) == subjectHits(hits)),
              identical(after, place(after, lengths)))
}
must_fail <- function(expr, pattern) {
    error <- tryCatch({force(expr); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

# A 1-bp overlap at the chromosome end: the old loop shifts the right interval
# out of bounds, then moves it straight back, repeating until its 100-pass limit.
x <- GRanges("chr1", IRanges(c(200, 349), width = 150), strand = "+")
y <- place(x, c(chr1 = 498))
stopifnot(identical(start(y), c(199L, 349L)))
check(x, y, c(chr1 = 498))

# Both chromosome boundaries and a dense chain must be handled jointly.
for (starts in list(c(-30, 20, 80), c(10, 40, 100, 130), c(501, 600, 720))) {
    x <- GRanges("chr1", IRanges(starts, width = 150), strand = "+")
    check(x, place(x, c(chr1 = 800)), c(chr1 = 800))
}
x <- GRanges("chr1", IRanges(seq(1, by = 20, length.out = 2001), width = 150), strand = "-")
check(x, place(x, c(chr1 = 400000)), c(chr1 = 400000))

# Different chromosomes/strands are independent, including overlapping opposite
# strands. A lone out-of-bounds range is slid back without changing its width.
x <- GRanges(c("chr1", "chr2", "chr1", "chr2"),
    IRanges(c(200, -10, 349, 100), width = c(150, 40, 150, 50)),
    strand = c("+", "+", "-", "-"), label = letters[1:4])
y <- place(x, c(chr1 = 498, chr2 = 300))
check(x, y, c(chr1 = 498, chr2 = 300))
stopifnot(identical(width(x), width(y[match(x$label, y$label)])),
          start(y[y$label == "a"]) == 200L, start(y[y$label == "c"]) == 349L,
          start(y[y$label == "b"]) == 1L)
empty <- GRanges()
stopifnot(identical(place(empty, c(chr1 = 100)), empty))
must_fail(place(GRanges("chr1", IRanges(c(1, 100), width = 150), strand = "+"),
                c(chr1 = 299)), "combined width 300 bp exceeds chromosome length 299")
must_fail(place(GRanges("chr1", IRanges(1, width = 150), strand = "+"),
                c(chr1 = 149)), "Cannot place")
must_fail(place(GRanges("chr1", IRanges(1, width = 150), strand = "*"),
                c(chr1 = 1000)), "explicit")

# Compare small variable-width problems with an exhaustive integer optimum.
# This checks the placement objective independently of the regression algorithm.
set.seed(1109)
for (trial in 1:30) {
    w <- sample(1:3, 3, replace = TRUE)
    len <- sum(w) + sample(0:3, 1)
    desired <- sort(sample(-2:(len + 2), 3))
    x <- GRanges("chr1", IRanges(desired, width = w), strand = "+")
    y <- place(x, c(chr1 = len))
    check(x, y, c(chr1 = len))
    candidates <- as.matrix(expand.grid(lapply(w, function(v) seq_len(len - v + 1))))
    feasible <- candidates[, 2] >= candidates[, 1] + w[1] &
                candidates[, 3] >= candidates[, 2] + w[2]
    costs <- rowSums(sweep(candidates[feasible, , drop = FALSE], 2, desired, "-")^2)
    stopifnot(sum((start(y) - desired)^2) == min(costs))
}

# Exercise the complete production script: original 1-bp calls expand to the
# stalled boundary pair above. The consensus must finish and export valid BED.
tmp <- tempfile("overlap-regression-"); dir.create(tmp)
writeLines(c("chr1\t274\t275\ta\t0\t+", "chr1\t423\t424\tb\t0\t+"),
           file.path(tmp, "sample_csrna1.tss.bed"))
writeLines("chr1\t498", file.path(tmp, "chrom.sizes"))
writeLines(c("Sample\tStatus", "sample_csrna1\tOk"), file.path(tmp, "qc.tsv"))
setClass("OverlapTest", slots = c(input = "list", output = "list", params = "list"))
env <- new.env(parent = globalenv())
env$snakemake <- new("OverlapTest", input = list(qc_cs = file.path(tmp, "qc.tsv")),
    output = list(file.path(tmp, "consensus.bed")), params = list(
        ids = "sample_csrna1", sample_names = "sample", sample_types = "csrna",
        replicates = "1", min_reps = 1L, exclude_failed = FALSE, bed_dir = tmp,
        chrom_sizes = file.path(tmp, "chrom.sizes"), mirnas = "", trnas = ""))
sys.source("workflow/scripts/collect_consensus_tss.R", envir = env)
result <- import(file.path(tmp, "consensus.bed"))
stopifnot(identical(start(result), c(199L, 349L)), all(width(result) == 150L),
          identical(result$name, c("TSS_1", "TSS_2")))
unlink(tmp, recursive = TRUE)
cat("Overlap resolution regressions passed\n")
