suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

cs_ids        <- strsplit(snakemake@params[["cs_ids"]],        " ", fixed = TRUE)[[1]]
paired_in_ids <- strsplit(snakemake@params[["paired_in_ids"]], " ", fixed = TRUE)[[1]]
cs_to_in      <- setNames(paired_in_ids, cs_ids)

# Gates come from params, not snakemake@config, so that editing one is a change Snakemake
# can see; and they are the final gates (qc/min_cs_frip_final and min_pct_nuclear_final,
# each falling back to its initial counterpart), resolved in the Snakefile rather than
# here so the resolved number is what triggers a rerun.
MIN_CS_FRIP     <- as.numeric(snakemake@params[["min_cs_frip"]])
MIN_PCT_NUCLEAR <- as.numeric(snakemake@params[["min_pct_nuclear"]])
# Not a final-specific gate: the threshold HOMER chose is a property of the library rather
# than of the TSS set being measured, so both tables apply the same floor and agree on why
# a library failed.
MIN_LOG2_FOLD   <- as.numeric(snakemake@params[["min_log2_fold"]])
TSS_DIR         <- trimws(snakemake@params[["tss_dir"]])

# The enrichment threshold HOMER chose for this library, out of its own stats file. With
# an annotation it picks the threshold from the data per library and there is no flag to
# floor it, so reading it back is the only way to gate on it. The line is always present
# when -i was given, which this pipeline always does: with nothing to derive a threshold
# from, HOMER falls back to -defaultLog2Fold and writes that value here, so the number is
# the threshold that was actually in force either way. NA when the file is missing or
# truncated, which means TSS calling for that library did not finish, and counts as a
# failure below rather than as a pass.
read_log2_threshold <- function(sample, dir) {
    path <- file.path(dir, paste0(sample, ".stats.txt"))
    if (!file.exists(path)) return(NA_real_)
    hit <- grep("log2 fold vs. input:", readLines(path, warn = FALSE),
                fixed = TRUE, value = TRUE)
    if (!length(hit)) return(NA_real_)
    suppressWarnings(as.numeric(
        sub(".*log2 fold vs\\. input:[[:space:]]*", "", hit[1])
    ))
}

stats_cs <- read.table(trimws(snakemake@input[["stats_cs"]]), header = TRUE, stringsAsFactors = FALSE)
stats_in <- read.table(trimws(snakemake@input[["stats_in"]]), header = TRUE, stringsAsFactors = FALSE)

m <- cs_to_in[stats_cs[["Sample"]]]
if (anyNA(m)) {
  stop("ERROR: could not find paired input for csRNA sample(s): ",
       paste(stats_cs[["Sample"]][is.na(m)], collapse = ", "))
}
stats_in_cs <- stats_in[match(m, stats_in[["Sample"]]), ]

stats_cs[["NuclearReads"]]    <- with(stats_cs, TotalReads - OrganelleReads)
stats_in[["NuclearReads"]]    <- with(stats_in, TotalReads - OrganelleReads)
stats_in_cs[["NuclearReads"]] <- with(stats_in_cs, TotalReads - OrganelleReads)

# Helper: extract tag-count matrix and strip HOMER column-name decoration.
# Rownames are set to the first column (TSS ID) so downstream filtering by name works.
read_homer_quant <- function(path) {
    df  <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(path))))
    ids <- df[[1]]
    idx <- grep(" Tag Count", colnames(df))
    m   <- as.matrix(df[, idx, drop = FALSE])
    colnames(m) <- basename(gsub(" Tag Count .+", "", colnames(m)))
    rownames(m) <- ids
    m
}

# csRNA quant on final TSS set (clean format: TSS-id + integer count columns)
quant_cs_raw <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(snakemake@input[["quant_cs"]]))))
quant_cs_m   <- as.matrix(quant_cs_raw[, -1, drop = FALSE])
storage.mode(quant_cs_m) <- "double"
rownames(quant_cs_m) <- quant_cs_raw[[1]]

# Input quant on final TSS set (HOMER format)
quant_in_m <- read_homer_quant(snakemake@input[["quant_in"]])

# csRNA and input quant on input TSS set (for sRiP / sDepletion and miRNA/tRNA metrics)
quant_inTSS_cs_m <- read_homer_quant(snakemake@input[["quant_in_cs"]])
quant_inTSS_in_m <- read_homer_quant(snakemake@input[["quant_in_in"]])

# ── Metrics re-derived from the final csRNA TSS set ─────────────────────────
#
# Note on naming: in THIS file the csRiP / csFRiP / csEnrichment columns are computed
# against tss.final.bed, whereas the identically named columns in qc_initial_*.txt are
# computed against the initial merged TSS set.  The names are kept so that the same
# column means "reads in the TSS set this file is about" in both, but that does mean the
# two files must not be compared column-by-column as though they measured the same
# regions.  Earlier versions also wrote FinalRiP / FinalFRiP / FinalEnrichment as literal
# copies of these three; those aliases are gone.

stats_cs[["csRiP"]]  <- colSums(quant_cs_m)[stats_cs[["Sample"]]]
stats_in[["csRiP"]]  <- colSums(quant_in_m)[stats_in[["Sample"]]]
stats_in_cs[["csRiP"]] <- colSums(quant_in_m)[stats_in_cs[["Sample"]]]

stats_cs[["csFRiP"]]    <- stats_cs[["csRiP"]]    / stats_cs[["NuclearReads"]]
stats_in[["csFRiP"]]    <- stats_in[["csRiP"]]    / stats_in[["NuclearReads"]]
stats_in_cs[["csFRiP"]] <- stats_in_cs[["csRiP"]] / stats_in_cs[["NuclearReads"]]

stats_cs[["csRNACappedPct"]] <- 100 * stats_cs[["csFRiP"]]
stats_cs[["csEnrichment"]]   <- stats_cs[["csFRiP"]] / stats_in_cs[["csFRiP"]]

# ── Metrics re-derived from the input TSS set ────────────────────────────────
# (input TSSs are unaffected by the CPM filter, but we re-derive here from
# raw quantification rather than copying from the initial QC output)

stats_cs[["sRiP"]]  <- colSums(quant_inTSS_cs_m)[stats_cs[["Sample"]]]
stats_in[["sRiP"]]  <- colSums(quant_inTSS_in_m)[stats_in[["Sample"]]]
stats_in_cs[["sRiP"]] <- colSums(quant_inTSS_in_m)[stats_in_cs[["Sample"]]]

stats_cs[["sFRiP"]]    <- stats_cs[["sRiP"]]    / stats_cs[["NuclearReads"]]
stats_in[["sFRiP"]]    <- stats_in[["sRiP"]]    / stats_in[["NuclearReads"]]
stats_in_cs[["sFRiP"]] <- stats_in_cs[["sRiP"]] / stats_in_cs[["NuclearReads"]]

stats_cs[["sDepletion"]] <- 1 / (stats_cs[["sFRiP"]] / stats_in_cs[["sFRiP"]])

# ── FinalTSSDetected ─────────────────────────────────────────────────────────

stats_cs[["FinalTSSDetected"]] <- colMeans(quant_cs_m > 0)[stats_cs[["Sample"]]]
stats_in[["FinalTSSDetected"]] <- colMeans(quant_in_m > 0)[stats_in[["Sample"]]]

# ── TSS filter summary ───────────────────────────────────────────────────────

tss_final     <- import(trimws(snakemake@input[["tss_final"]]))
tss_consensus <- import(trimws(snakemake@input[["tss_consensus"]]))
n_final     <- length(tss_final)
n_consensus <- length(tss_consensus)
stats_cs[["NConsensusTSS"]] <- n_consensus
stats_cs[["NFinalTSS"]]     <- n_final
stats_cs[["NFilteredTSS"]]  <- n_consensus - n_final
stats_in[["NConsensusTSS"]] <- n_consensus
stats_in[["NFinalTSS"]]     <- n_final
stats_in[["NFilteredTSS"]]  <- n_consensus - n_final

# ── Replicate correlation on final TSS counts ────────────────────────────────

stats_cs[["group"]] <- gsub("_csrna[0-9]+$", "", stats_cs[["Sample"]])
stats_cs[["MinReplCorrSpearman"]] <- NA_real_
stats_cs[["MinReplCorrPearson"]]  <- NA_real_
repl_cor_long <- list()
for (k in unique(stats_cs[["group"]])) {
    samples_k <- stats_cs[["Sample"]][stats_cs[["group"]] == k]
    if (length(samples_k) < 2) next
    m_k    <- as.matrix(quant_cs_m[, samples_k, drop = FALSE])
    cor_sp <- cor(m_k, method = "spearman")
    cor_pe <- cor(log1p(m_k), method = "pearson")
    diag(cor_sp) <- NA
    diag(cor_pe) <- NA
    stats_cs[stats_cs[["group"]] == k, "MinReplCorrSpearman"] <-
        apply(cor_sp, 1, min, na.rm = TRUE)[samples_k]
    stats_cs[stats_cs[["group"]] == k, "MinReplCorrPearson"] <-
        apply(cor_pe, 1, min, na.rm = TRUE)[samples_k]
    pairs <- which(upper.tri(cor_sp), arr.ind = TRUE)
    if (nrow(pairs)) {
        repl_cor_long[[k]] <- data.frame(
            Group      = k,
            SampleA    = rownames(cor_sp)[pairs[, 1]],
            SampleB    = colnames(cor_sp)[pairs[, 2]],
            SpearmanCor = cor_sp[pairs],
            PearsonCor  = cor_pe[pairs],
            stringsAsFactors = FALSE
        )
    }
}
if (length(repl_cor_long)) {
    repl_cor_final <- do.call(rbind, repl_cor_long)
} else {
    repl_cor_final <- data.frame(
        Group = character(), SampleA = character(), SampleB = character(),
        SpearmanCor = numeric(), PearsonCor = numeric(),
        stringsAsFactors = FALSE
    )
}
repl_cor_final[["Stage"]] <- "Final"

repl_cor_initial <- suppressWarnings(suppressMessages(
    readr::read_tsv(trimws(snakemake@input[["repl_cor_initial"]]))
))
repl_cor_combined <- rbind(repl_cor_initial, repl_cor_final)
readr::write_tsv(repl_cor_combined, snakemake@output[["repl_cor"]])

# ── Optional: miRNA and pre-tRNA contamination metrics ──────────────────────

tss_in    <- import(trimws(snakemake@input[["tss_in"]]))
f_mirna   <- trimws(snakemake@params[["mirnas"]])
f_trna    <- trimws(snakemake@params[["trnas"]])
mirna     <- if (nchar(f_mirna)) import(f_mirna) else NULL
trna      <- if (nchar(f_trna))  import(f_trna)  else NULL

if (!is.null(trna)) {
    in_pretrna <- mcols(tss_in)[["name"]][overlapsAny(tss_in, trna)]
    pretrna_rows <- rownames(quant_inTSS_cs_m) %in% in_pretrna
    stats_cs[["PretRNA"]]     <- colSums(quant_inTSS_cs_m[pretrna_rows, , drop = FALSE])[stats_cs[["Sample"]]]
    stats_in[["PretRNA"]]     <- colSums(quant_inTSS_in_m[pretrna_rows, , drop = FALSE])[stats_in[["Sample"]]]
    in_cs_pretrna             <- colSums(quant_inTSS_in_m[pretrna_rows, , drop = FALSE])[stats_in_cs[["Sample"]]]
    stats_cs[["PretRNAPct"]]  <- with(stats_cs, 100 * (PretRNA / (PretRNA + csRiP)))
    stats_in[["PretRNAPct"]]  <- with(stats_in, 100 * (PretRNA / (PretRNA + csRiP)))
    in_cs_pretrnaPct          <- 100 * (in_cs_pretrna / (in_cs_pretrna + stats_in_cs[["csRiP"]]))
    stats_cs[["PhosEfficiency"]] <- 1 / (stats_cs[["PretRNAPct"]] / in_cs_pretrnaPct)
}

if (!is.null(mirna)) {
    in_mirna <- mcols(tss_in)[["name"]][overlapsAny(tss_in, mirna)]
    mirna_rows <- rownames(quant_inTSS_cs_m) %in% in_mirna
    stats_cs[["miRNA"]]        <- colSums(quant_inTSS_cs_m[mirna_rows, , drop = FALSE])[stats_cs[["Sample"]]]
    stats_in[["miRNA"]]        <- colSums(quant_inTSS_in_m[mirna_rows, , drop = FALSE])[stats_in[["Sample"]]]
    in_cs_mirna                <- colSums(quant_inTSS_in_m[mirna_rows, , drop = FALSE])[stats_in_cs[["Sample"]]]
    stats_cs[["miRNADepletion"]] <- 1 / (
        (stats_cs[["miRNA"]] / stats_cs[["NuclearReads"]]) /
        (in_cs_mirna          / stats_in_cs[["NuclearReads"]])
    )
    # Same basis as PretRNAPct: contaminant reads as a share of contaminant plus
    # reads-in-TSS, so the two contamination percentages are comparable.
    stats_cs[["miRNAPct"]] <- with(stats_cs, 100 * (miRNA / (miRNA + csRiP)))
    stats_in[["miRNAPct"]] <- with(stats_in, 100 * (miRNA / (miRNA + csRiP)))
}

# ── Remaining per-library metrics ────────────────────────────────────────────

stats_cs[["StrandBalance"]] <- with(stats_cs, PosReads / (PosReads + NegReads))
stats_in[["StrandBalance"]] <- with(stats_in, PosReads / (PosReads + NegReads))
stats_cs[["PctNuclear"]] <- 100 * (stats_cs[["NuclearReads"]] / stats_cs[["TotalReads"]])
stats_in[["PctNuclear"]] <- 100 * (stats_in[["NuclearReads"]] / stats_in[["TotalReads"]])
stats_cs[["Log2FoldThreshold"]] <- vapply(stats_cs[["Sample"]], read_log2_threshold,
                                          numeric(1), dir = TSS_DIR, USE.NAMES = FALSE)
stats_in[["Log2FoldThreshold"]] <- vapply(stats_in[["Sample"]], read_log2_threshold,
                                          numeric(1), dir = TSS_DIR, USE.NAMES = FALSE)

# Same gate-to-Status logic as qc_initial_tss.R, against the final numbers. StatusReason
# names the gates that tripped; an NA gate counts as a failure, since a gate that cannot
# be evaluated is not evidence that the library is sound.
qc_status <- function(gates) {
    fails <- do.call(cbind, lapply(gates, function(g) is.na(g) | !g))
    colnames(fails) <- names(gates)
    data.frame(
        Status = ifelse(rowSums(fails) > 0, "FAIL", "Ok"),
        StatusReason = apply(fails, 1, function(r) paste(names(gates)[r], collapse = ",")),
        stringsAsFactors = FALSE
    )
}

status <- qc_status(list(
    csFRiP            = stats_cs[["csFRiP"]]     > MIN_CS_FRIP,
    PctNuclear        = stats_cs[["PctNuclear"]] > MIN_PCT_NUCLEAR,
    Log2FoldThreshold = stats_cs[["Log2FoldThreshold"]] >= MIN_LOG2_FOLD
))
stats_cs[["Status"]] <- status[["Status"]]
stats_cs[["StatusReason"]] <- status[["StatusReason"]]
# Advisory, and labelled as such: collect_consensus_tss reads the initial Status, so a
# FAIL here has already had no effect on the TSS set this table describes.
stats_cs[["StatusBasis"]] <- "final-advisory"

readr::write_tsv(stats_cs, snakemake@output[["qc_cs"]])
readr::write_tsv(stats_in, snakemake@output[["qc_in"]])

cat("Final TSS set: ", n_final, " of ", n_consensus,
    " consensus TSSs retained (", n_consensus - n_final, " filtered)\n", sep = "")
cat("Printing summary stats:\n")
print(cbind(stats_cs[, c("Sample", "Status", "StatusReason", "PctNuclear", "csFRiP",
                        "Log2FoldThreshold")],
            csFRiP_in = stats_in_cs[["csFRiP"]]))
