# Exercise the production scripts without requiring a workflow engine.
suppressPackageStartupMessages(library(rtracklayer))
setClass("FilterTest", slots = c(input = "list", output = "list", params = "list"))
tmp <- tempfile("filter-regressions-")
dir.create(tmp)
fixtures <- "tests/data/unit_fixtures"
run_script <- function(name, input, output, params) {
    env <- new.env(parent = globalenv())
    env$snakemake <- new("FilterTest", input = input, output = output, params = params)
    sys.source(paste0("workflow/scripts/", name, ".R"), envir = env)
}
must_fail <- function(expr, pattern) {
    error <- tryCatch({force(expr); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
cons_params <- list(ids = "condA_csrna1 condA_csrna2 condB_csrna1",
    sample_names = "condA condA condB", sample_types = "csrna csrna csrna",
    replicates = "1 2 1", min_reps = 2L, exclude_failed = FALSE,
    bed_dir = paste0(fixtures, "/per_sample_tss"), chrom_sizes = "tests/data/chrom.sizes",
    mirnas = "", trnas = "")
cons_input <- list(qc_cs = paste0(fixtures, "/qc_initial_cs.txt"))
cons_output <- list(file.path(tmp, "consensus.bed"))
run_consensus <- function() run_script("collect_consensus_tss", cons_input, cons_output, cons_params)
run_consensus()
bed <- import(cons_output[[1]])
# The condB-only locus at 7000 must not enter a two-replicate consensus.
stopifnot(length(bed) == 6L, !any(start(bed) == 7001L))
qc <- read.delim(cons_input$qc_cs)
qc$Status[qc$Sample == "condA_csrna2"] <- "FAIL"
cons_input$qc_cs <- file.path(tmp, "qc.txt")
write.table(qc, cons_input$qc_cs, sep = "\t", row.names = FALSE, quote = FALSE)
cons_params$exclude_failed <- TRUE
must_fail(run_consensus(), "no TSSs meet tss_min_reps")
cons_params$min_reps <- 1L
run_consensus()
stopifnot(length(import(cons_output[[1]])) == 6L)

cs <- read.delim(paste0(fixtures, "/tss.consensus.homer.raw.txt"), check.names = FALSE)
inp <- read.delim(paste0(fixtures, "/tss.consensus.in.homer.raw.txt"), check.names = FALSE)
cscols <- grep(" Tag Count", names(cs))
incols <- grep(" Tag Count", names(inp))
cs[, cscols] <- 10
inp[, incols] <- 10
cs[1, cscols] <- 40; inp[1, incols] <- 2   # normalized ratio 2: retained
cs[2, cscols] <- 20; inp[2, incols] <- 4   # normalized ratio 0.5: removed
cs[3, cscols] <- 5; inp[3, incols] <- 0    # floor: ratio 0.5 (1 after rescaling), removed
depth_cs <- read.delim(paste0(fixtures, "/stats_cs.txt"))
depth_in <- read.delim(paste0(fixtures, "/stats_in.txt"))
depth_cs$TotalReads <- 1100; depth_cs$OrganelleReads <- 100
depth_in$TotalReads <- 120; depth_in$OrganelleReads <- 20
inputs <- list(quant = file.path(tmp, "cs.txt"), quant_in = file.path(tmp, "in.txt"),
    stats_cs = file.path(tmp, "depth_cs.txt"), stats_in = file.path(tmp, "depth_in.txt"),
    bed = paste0(fixtures, "/tss.consensus.bed"),
    sizes = paste0(fixtures, "/tss.consensus.sizes.txt"))
outputs <- list(bed = file.path(tmp, "final.bed"), raw = file.path(tmp, "raw.txt"),
    norm = file.path(tmp, "cpm.txt"), normfactors = file.path(tmp, "factors.txt"))
params <- list(min_cpm = 0, min_samples = 1, srna_sizes = "", max_srna_fraction = 0.5,
    srna_min_reads = 100, srna_min_samples = 1, top_sizes_n = 2,
    max_top_sizes_fraction = 1, min_cs_in_ratio = 1.5, min_ratio_samples = 2,
    cs_ids = "condA_csrna1 condA_csrna2 condB_csrna1",
    paired_in_ids = "condA_input1 condA_input2 condA_input1")
save_tables <- function() {
    for (entry in list(list(cs, inputs$quant), list(inp, inputs$quant_in),
                      list(depth_cs, inputs$stats_cs), list(depth_in, inputs$stats_in)))
        write.table(entry[[1]], entry[[2]], sep = "\t", row.names = FALSE, quote = FALSE)
}
run_normalize <- function() run_script("normalize_tss_quantification", inputs, outputs, params)
save_tables(); run_normalize()
stopifnot(identical(read.delim(outputs$raw)$TSS, "TSS_1"))
# Shuffled input rows are safe; doubling input depth and its nonzero counts must
# preserve the enriched locus. The one-read floor means exact scale invariance
# is not expected at zero input. Shared input pairing is exercised by condB.
inp <- inp[nrow(inp):1, ]
inp[, incols] <- inp[, incols] * 2
depth_in$TotalReads <- 220
save_tables(); run_normalize()
stopifnot(identical(read.delim(outputs$raw)$TSS, "TSS_1"))
inp[1, 3] <- inp[1, 3] + 1
save_tables(); must_fail(run_normalize(), "coordinates/strands disagree")
inp[1, 3] <- inp[1, 3] - 1
good_inp <- inp
inp <- inp[-1, ]
save_tables(); must_fail(run_normalize(), "same unique consensus TSS IDs")
inp <- good_inp
# A missing/non-positive denominator is not silently treated as an enriched locus.
depth_in$TotalReads <- depth_in$OrganelleReads
save_tables(); must_fail(run_normalize(), "positive finite non-organelle library depths")
depth_in$TotalReads <- 220
# A single csRNA library must preserve its column name/dimensions.
cs <- cs[, c(seq_len(cscols[1] - 1L), cscols[1]), drop = FALSE]
params$cs_ids <- "condA_csrna1"
params$paired_in_ids <- "condA_input1"
params$min_ratio_samples <- 1
save_tables(); run_normalize()
stopifnot(identical(read.delim(outputs$raw)$TSS, "TSS_1"))
# Disabled ratio filter does not read stale/missing ratio inputs.
params$min_cs_in_ratio <- 0
run_normalize()
stopifnot(nrow(read.delim(outputs$raw)) == 10L)
unlink(tmp, recursive = TRUE)
cat("Filter regression tests passed\n")
