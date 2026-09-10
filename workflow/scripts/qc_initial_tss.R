suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

QUANT_csTSS_csRNA <- snakemake@input[["quant_cs_cs"]]
QUANT_csTSS_input <- snakemake@input[["quant_cs_in"]]
QUANT_inTSS_csRNA <- snakemake@input[["quant_in_cs"]]
QUANT_inTSS_input <- snakemake@input[["quant_in_in"]]
STATS_CS <- snakemake@input[["stats_cs"]]
STATS_IN <- snakemake@input[["stats_in"]]
TSS_CS <- snakemake@input[["tss_cs"]]
TSS_IN <- snakemake@input[["tss_in"]]
MIRNAS <- snakemake@params[["mirnas"]]
TRNAS <- snakemake@params[["trnas"]]
# The QC gates arrive as params rather than through snakemake@config, so that editing one
# in the config is a change Snakemake can act on. A config value read directly inside the
# script is invisible to every rerun trigger, so the old Status would stay on disk with
# nothing in the log to say the new gate had been ignored.
MIN_CS_FRIP <- as.numeric(snakemake@params[["min_cs_frip"]])
MIN_PCT_NUCLEAR <- as.numeric(snakemake@params[["min_pct_nuclear"]])
MIN_LOG2_FOLD <- as.numeric(snakemake@params[["min_log2_fold"]])
TSS_DIR <- trimws(snakemake@params[["tss_dir"]])
OUT_CS <- snakemake@output[["qc_cs"]]
OUT_IN <- snakemake@output[["qc_in"]]
OUT_REPLCOR <- snakemake@output[["repl_cor"]]

tss_cs <- import(trimws(TSS_CS))
tss_in <- import(trimws(TSS_IN))

# Preserve empty TSV fields from older stats files as missing values.
stats_cs <- read.table(trimws(STATS_CS), header=TRUE, sep="\t", na.strings=c("", "NA", "na"), stringsAsFactors=FALSE)
stats_in <- read.table(trimws(STATS_IN), header=TRUE, sep="\t", na.strings=c("", "NA", "na"), stringsAsFactors=FALSE)

# Build the cs->in pairing from the explicit id lists passed by Snakemake.
# This handles csRNA samples that share an input with a different sample_name
# (via the input_name column in the sample sheet).
cs_ids <- strsplit(snakemake@params[["cs_ids"]], " ", fixed = TRUE)[[1]]
paired_in_ids <- strsplit(snakemake@params[["paired_in_ids"]], " ", fixed = TRUE)[[1]]
cs_to_in <- setNames(paired_in_ids, cs_ids)

m <- cs_to_in[stats_cs[['Sample']]]
if (anyNA(m)) {
  stop("ERROR: could not find paired input for csRNA sample(s): ",
       paste(stats_cs[['Sample']][is.na(m)], collapse = ", "))
}
stats_in_cs <- stats_in[match(m, stats_in[['Sample']]), ]

# group: csRNA sample_name group (for replicate correlation grouping).
# Replicates of the same condition share sample_name, so strip the type+rep suffix.
stats_cs[['group']] <- gsub("_csrna[0-9]+$", "", stats_cs[['Sample']])

quant_csTSS_csRNA <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(QUANT_csTSS_csRNA))))
quant_csTSS_input <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(QUANT_csTSS_input))))
quant_inTSS_csRNA <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(QUANT_inTSS_csRNA))))
quant_inTSS_input <- suppressWarnings(suppressMessages(readr::read_tsv(trimws(QUANT_inTSS_input))))

f_mirna <- trimws(MIRNAS)
f_trna <- trimws(TRNAS)

mirna <- if (nchar(f_mirna)) import(f_mirna) else NULL
trna <- if (nchar(f_trna)) import(f_trna) else NULL

stats_cs[['NuclearReads']] <- with(stats_cs, TotalReads - OrganelleReads)
stats_in[['NuclearReads']] <- with(stats_in, TotalReads - OrganelleReads)
stats_cs[['PctNuclear']] <- 100 * (stats_cs[['NuclearReads']] / stats_cs[['TotalReads']])
stats_in[['PctNuclear']] <- 100 * (stats_in[['NuclearReads']] / stats_in[['TotalReads']])
stats_in_cs[['NuclearReads']] <- with(stats_in_cs, TotalReads - OrganelleReads)

# Detect tag-count columns by HOMER's column name pattern rather than hard-coded index
tag_cols_csTSS_csRNA <- grep(" Tag Count", colnames(quant_csTSS_csRNA))
tag_cols_csTSS_input <- grep(" Tag Count", colnames(quant_csTSS_input))
tag_cols_inTSS_csRNA <- grep(" Tag Count", colnames(quant_inTSS_csRNA))
tag_cols_inTSS_input <- grep(" Tag Count", colnames(quant_inTSS_input))

quant_csTSS_csRNA_m <- quant_csTSS_csRNA[, tag_cols_csTSS_csRNA]
colnames(quant_csTSS_csRNA_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_csTSS_csRNA_m)))

quant_csTSS_input_m <- quant_csTSS_input[, tag_cols_csTSS_input]
colnames(quant_csTSS_input_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_csTSS_input_m)))

quant_inTSS_csRNA_m <- quant_inTSS_csRNA[, tag_cols_inTSS_csRNA]
colnames(quant_inTSS_csRNA_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_inTSS_csRNA_m)))

quant_inTSS_input_m <- quant_inTSS_input[, tag_cols_inTSS_input]
colnames(quant_inTSS_input_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_inTSS_input_m)))

stats_cs[['csRiP']] <- colSums(quant_csTSS_csRNA_m)[stats_cs[['Sample']]]
stats_cs[['sRiP']] <- colSums(quant_inTSS_csRNA_m)[stats_cs[['Sample']]]

stats_in[['csRiP']] <- colSums(quant_csTSS_input_m)[stats_in[['Sample']]]
stats_in[['sRiP']] <- colSums(quant_inTSS_input_m)[stats_in[['Sample']]]
stats_in_cs[['csRiP']] <- colSums(quant_csTSS_input_m)[stats_in_cs[['Sample']]]
stats_in_cs[['sRiP']] <- colSums(quant_inTSS_input_m)[stats_in_cs[['Sample']]]

stats_cs[['csFRiP']] <- stats_cs[['csRiP']] / stats_cs[['NuclearReads']]
stats_cs[['sFRiP']] <- stats_cs[['sRiP']] / stats_cs[['NuclearReads']]
stats_in[['csFRiP']] <- stats_in[['csRiP']] / stats_in[['NuclearReads']]
stats_in[['sFRiP']] <- stats_in[['sRiP']] / stats_in[['NuclearReads']]
stats_in_cs[['csFRiP']] <- stats_in_cs[['csRiP']] / stats_in_cs[['NuclearReads']]
stats_in_cs[['sFRiP']] <- stats_in_cs[['sRiP']] / stats_in_cs[['NuclearReads']]

stats_cs[['csRNACappedPct']] <- 100 * stats_cs[['csFRiP']]
stats_cs[['InputCappedPct']] <- 100 * stats_in_cs[['csFRiP']]
stats_cs[['csEnrichment']] <- stats_cs[['csFRiP']] / stats_in_cs[['csFRiP']]
stats_cs[['sDepletion']] <- 1 / (stats_cs[['sFRiP']] / stats_in_cs[['sFRiP']])

# Strand balance: fraction of reads on the + strand. ~0.5 in csRNA libraries;
# can deviate from 0.5 in input libraries because small-RNA biology is genuinely
# strand-skewed at a few highly expressed loci.
stats_cs[['StrandBalance']] <- with(stats_cs, PosReads / (PosReads + NegReads))
stats_in[['StrandBalance']] <- with(stats_in, PosReads / (PosReads + NegReads))

# TSS saturation: fraction of merged-set TSSs detected (>=1 tag) in this library.
# Low values flag undersequenced libraries.
stats_cs[['TSSDetected']] <- colMeans(as.matrix(quant_csTSS_csRNA_m) > 0)[stats_cs[['Sample']]]
stats_in[['TSSDetected']] <- colMeans(as.matrix(quant_inTSS_input_m) > 0)[stats_in[['Sample']]]

# Replicate correlation: per-sample minimum Spearman and Pearson (log1p) correlation
# with another csRNA replicate of the same sample_name. NA when only one replicate exists.
stats_cs[['MinReplCorrSpearman']] <- NA_real_
stats_cs[['MinReplCorrPearson']]  <- NA_real_
repl_cor_long <- list()
for (k in unique(stats_cs[['group']])) {
    samples_k <- stats_cs[['Sample']][stats_cs[['group']] == k]
    if (length(samples_k) < 2) next
    m_k <- as.matrix(quant_csTSS_csRNA_m[, samples_k])
    cor_sp <- cor(m_k, method = 'spearman')
    cor_pe <- cor(log1p(m_k), method = 'pearson')
    diag(cor_sp) <- NA
    diag(cor_pe) <- NA
    stats_cs[stats_cs[['group']] == k, 'MinReplCorrSpearman'] <-
        apply(cor_sp, 1, min, na.rm = TRUE)[samples_k]
    stats_cs[stats_cs[['group']] == k, 'MinReplCorrPearson'] <-
        apply(cor_pe, 1, min, na.rm = TRUE)[samples_k]
    pairs <- which(upper.tri(cor_sp), arr.ind = TRUE)
    if (nrow(pairs)) {
        repl_cor_long[[k]] <- data.frame(
            Group = k,
            SampleA = rownames(cor_sp)[pairs[, 1]],
            SampleB = colnames(cor_sp)[pairs[, 2]],
            SpearmanCor = cor_sp[pairs],
            PearsonCor  = cor_pe[pairs],
            stringsAsFactors = FALSE
        )
    }
}
if (length(repl_cor_long)) {
    repl_cor_df <- do.call(rbind, repl_cor_long)
} else {
    repl_cor_df <- data.frame(
        Group = character(),
        SampleA = character(),
        SampleB = character(),
        SpearmanCor = numeric(),
        PearsonCor  = numeric(),
        stringsAsFactors = FALSE
    )
}
repl_cor_df[['Stage']] <- 'Initial'
readr::write_tsv(repl_cor_df, OUT_REPLCOR)

if (!is.null(trna)) {
    in_pretrna <- mcols(tss_in)[['name']][overlapsAny(tss_in, trna)]
    stats_cs[['PretRNA']] <- colSums(quant_inTSS_csRNA_m[quant_inTSS_csRNA[[1]] %in% in_pretrna, ])[stats_cs[['Sample']]]
    stats_in[['PretRNA']] <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_pretrna, ])[stats_in[['Sample']]]
    stats_in_cs_pretrna <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_pretrna, ])[stats_in_cs[['Sample']]]
    stats_cs[['PretRNAPct']] <- with(stats_cs, 100 * (PretRNA / (PretRNA + csRiP)))
    stats_in[['PretRNAPct']] <- with(stats_in, 100 * (PretRNA / (PretRNA + csRiP)))
    stats_in_cs_pretrnaPct <- 100 * (stats_in_cs_pretrna / (stats_in_cs_pretrna + stats_in_cs[['csRiP']]))
    stats_cs[['PhosEfficiency']] <- 1 / (stats_cs[['PretRNAPct']] / stats_in_cs_pretrnaPct)
}

if (!is.null(mirna)) {
    in_mirna <- mcols(tss_in)[['name']][overlapsAny(tss_in, mirna)]
    stats_cs[['miRNA']] <- colSums(quant_inTSS_csRNA_m[quant_inTSS_csRNA[[1]] %in% in_mirna, ])[stats_cs[['Sample']]]
    stats_in[['miRNA']] <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_mirna, ])[stats_in[['Sample']]]
    stats_in_cs_miRNA <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_mirna, ])[stats_in_cs[['Sample']]]
    stats_cs[['miRNADepletion']] <- 1 / ((stats_cs[['miRNA']] / stats_cs[['NuclearReads']]) / (stats_in_cs_miRNA / stats_in_cs[['NuclearReads']]))
    # Same basis as PretRNAPct: contaminant reads as a share of contaminant plus
    # reads-in-TSS, so the two contamination percentages are comparable.
    stats_cs[['miRNAPct']] <- with(stats_cs, 100 * (miRNA / (miRNA + csRiP)))
    stats_in[['miRNAPct']] <- with(stats_in, 100 * (miRNA / (miRNA + csRiP)))
}

# The enrichment threshold HOMER chose for this library, out of its own stats file. With
# an annotation it picks the threshold from the data per library and there is no flag to
# floor it, so reading it back is the only way to gate on it. The number is the threshold
# that was in force either way: with nothing to derive one from, HOMER falls back to
# -defaultLog2Fold and writes that value on the same line.
#
# NA when the line is absent, which happens for real. A library with no valid clusters
# takes findcsRNATSS through a division by zero on its way to the promoter-distal
# fraction, so the report stops before the threshold is ever written; the input libraries
# in the test panel do exactly this. NA counts as a failure below rather than as a pass,
# since a library that called nothing has not earned an Ok.
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

stats_cs[['Log2FoldThreshold']] <- vapply(stats_cs[['Sample']], read_log2_threshold,
                                          numeric(1), dir = TSS_DIR, USE.NAMES = FALSE)
stats_in[['Log2FoldThreshold']] <- vapply(stats_in[['Sample']], read_log2_threshold,
                                          numeric(1), dir = TSS_DIR, USE.NAMES = FALSE)

# Status is a bare Ok/FAIL, which says nothing about which gate tripped. StatusReason
# names them, so a library that fails is legible without re-deriving every gate by hand.
# A gate that cannot be evaluated counts as a failure: an NA is not evidence that the
# library is sound.
qc_status <- function(gates) {
    fails <- do.call(cbind, lapply(gates, function(g) is.na(g) | !g))
    colnames(fails) <- names(gates)
    data.frame(
        Status = ifelse(rowSums(fails) > 0, 'FAIL', 'Ok'),
        StatusReason = apply(fails, 1, function(r) paste(names(gates)[r], collapse = ',')),
        stringsAsFactors = FALSE
    )
}

status <- qc_status(list(
    csFRiP            = stats_cs[['csFRiP']]     > MIN_CS_FRIP,
    PctNuclear        = stats_cs[['PctNuclear']] > MIN_PCT_NUCLEAR,
    # >= rather than >, so a library using exactly the configured floor passes. That is
    # the common case: with no annotation every library lands on default_log2_fold, which
    # is also what this gate falls back to.
    Log2FoldThreshold = stats_cs[['Log2FoldThreshold']] >= MIN_LOG2_FOLD
))
# This is the Status that acts: collect_consensus_tss reads this table, not the final
# one, so a FAIL here can change the TSS set when exclude_failed_from_consensus is on.
stats_cs[['Status']] <- status[['Status']]
stats_cs[['StatusReason']] <- status[['StatusReason']]

readr::write_tsv(stats_cs, OUT_CS)
readr::write_tsv(stats_in, OUT_IN)

cat('Printing summary stats:\n')
print(cbind(stats_cs[, c('Sample', 'Status', 'StatusReason', 'PctNuclear', 'csFRiP', 'Log2FoldThreshold')], csFRiP_in = stats_in_cs[['csFRiP']]))
