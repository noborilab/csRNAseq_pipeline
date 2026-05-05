suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

cs_ids        <- strsplit(snakemake@params[["cs_ids"]],        " ", fixed = TRUE)[[1]]
paired_in_ids <- strsplit(snakemake@params[["paired_in_ids"]], " ", fixed = TRUE)[[1]]
cs_to_in      <- setNames(paired_in_ids, cs_ids)

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

# csRNA quant on final TSS set (clean format: TSS + integer count columns)
quant_cs <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(snakemake@input[["quant_cs"]]))))
quant_cs_m <- as.matrix(quant_cs[, -1, drop = FALSE])
storage.mode(quant_cs_m) <- "double"
rownames(quant_cs_m) <- quant_cs[[1]]

# Input quant on final TSS set (HOMER annotatePeaks.pl format)
quant_in_raw <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(snakemake@input[["quant_in"]]))))
tag_cols_in <- grep(" Tag Count", colnames(quant_in_raw))
quant_in_m <- as.matrix(quant_in_raw[, tag_cols_in, drop = FALSE])
colnames(quant_in_m) <- basename(gsub(" Tag Count .+", "", colnames(quant_in_m)))

# FRiP on final TSS set
stats_cs[["FinalRiP"]]    <- colSums(quant_cs_m)[stats_cs[["Sample"]]]
stats_in[["FinalRiP"]]    <- colSums(quant_in_m)[stats_in[["Sample"]]]
stats_in_cs[["FinalRiP"]] <- colSums(quant_in_m)[stats_in_cs[["Sample"]]]

stats_cs[["FinalFRiP"]]    <- stats_cs[["FinalRiP"]]    / stats_cs[["NuclearReads"]]
stats_in[["FinalFRiP"]]    <- stats_in[["FinalRiP"]]    / stats_in[["NuclearReads"]]
stats_in_cs[["FinalFRiP"]] <- stats_in_cs[["FinalRiP"]] / stats_in_cs[["NuclearReads"]]

stats_cs[["FinalEnrichment"]] <- stats_cs[["FinalFRiP"]] / stats_in_cs[["FinalFRiP"]]

# TSSDetected on final set (fraction of TSSs with >= 1 tag)
stats_cs[["FinalTSSDetected"]] <- colMeans(quant_cs_m > 0)[stats_cs[["Sample"]]]
stats_in[["FinalTSSDetected"]] <- colMeans(quant_in_m > 0)[stats_in[["Sample"]]]

# TSS filter summary
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

# Replicate correlation on final TSS counts
stats_cs[["key"]] <- gsub("_csrna[0-9]+$", "", stats_cs[["Sample"]])
stats_cs[["MinReplCorrSpearman"]] <- NA_real_
stats_cs[["MinReplCorrPearson"]]  <- NA_real_
repl_cor_long <- list()
for (k in unique(stats_cs[["key"]])) {
    samples_k <- stats_cs[["Sample"]][stats_cs[["key"]] == k]
    if (length(samples_k) < 2) next
    m_k    <- as.matrix(quant_cs_m[, samples_k, drop = FALSE])
    cor_sp <- cor(m_k, method = "spearman")
    cor_pe <- cor(log1p(m_k), method = "pearson")
    diag(cor_sp) <- NA
    diag(cor_pe) <- NA
    stats_cs[stats_cs[["key"]] == k, "MinReplCorrSpearman"] <-
        apply(cor_sp, 1, min, na.rm = TRUE)[samples_k]
    stats_cs[stats_cs[["key"]] == k, "MinReplCorrPearson"] <-
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

# Combine with initial-stage correlations and write
repl_cor_initial <- suppressWarnings(suppressMessages(
    readr::read_tsv(trimws(snakemake@input[["repl_cor_initial"]]))
))
repl_cor_combined <- rbind(repl_cor_initial, repl_cor_final)
readr::write_tsv(repl_cor_combined, snakemake@output[["repl_cor"]])

stats_cs[["StrandBalance"]] <- with(stats_cs, PosReads / (PosReads + NegReads))
stats_in[["StrandBalance"]] <- with(stats_in, PosReads / (PosReads + NegReads))
stats_cs[["PctNuclear"]] <- 100 * (stats_cs[["NuclearReads"]] / stats_cs[["TotalReads"]])
stats_in[["PctNuclear"]] <- 100 * (stats_in[["NuclearReads"]] / stats_in[["TotalReads"]])

# Join library-quality columns from the initial QC output (computed against the
# merged initial TSS set; characterise library quality independently of the
# final TSS filter).
cs_extra <- c("csRiP", "sRiP", "csFRiP", "sFRiP",
               "csRNACappedPct", "csEnrichment", "sDepletion",
               "PretRNA", "PretRNAPct", "PhosEfficiency",
               "miRNA", "miRNADepletion")
in_extra <- c("csRiP", "sRiP", "csFRiP", "sFRiP",
               "PretRNA", "PretRNAPct", "miRNA")

qc_cs_init <- read.table(trimws(snakemake@input[["qc_cs_initial"]]),
                          header = TRUE, stringsAsFactors = FALSE)
qc_in_init <- read.table(trimws(snakemake@input[["qc_in_initial"]]),
                          header = TRUE, stringsAsFactors = FALSE)

join_cols <- function(target, source, cols) {
    cols_present <- intersect(cols, colnames(source))
    src <- source[match(target[["Sample"]], source[["Sample"]]), cols_present, drop = FALSE]
    cbind(target, src)
}

stats_cs <- join_cols(stats_cs, qc_cs_init, cs_extra)
stats_in <- join_cols(stats_in, qc_in_init, in_extra)

stats_cs[["Status"]] <- ifelse(
    stats_cs[["FinalFRiP"]] > snakemake@config[["qc"]][["min_cs_frip"]] &
    stats_cs[["PctNuclear"]] > snakemake@config[["qc"]][["min_pct_nuclear"]],
    "Ok", "FAIL"
)

readr::write_tsv(stats_cs, snakemake@output[["qc_cs"]])
readr::write_tsv(stats_in, snakemake@output[["qc_in"]])

cat("Final TSS set: ", n_final, " of ", n_consensus,
    " consensus TSSs retained (", n_consensus - n_final, " filtered)\n", sep = "")
cat("Printing summary stats:\n")
print(cbind(stats_cs[, c("Sample", "Status", "PctNuclear", "FinalFRiP")],
            FinalFRiP_in = stats_in_cs[["FinalFRiP"]]))
