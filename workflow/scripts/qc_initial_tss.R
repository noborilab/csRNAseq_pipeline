suppressMessages(suppressPackageStartupMessages(library(rtracklayer)))

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
OUT_CS <- snakemake@output[["qc_cs"]]
OUT_IN <- snakemake@output[["qc_in"]]

tss_cs <- import(trimws(TSS_CS))
tss_in <- import(trimws(TSS_IN))

stats_cs <- read.table(trimws(STATS_CS), header=TRUE, stringsAsFactors=FALSE)
stats_in <- read.table(trimws(STATS_IN), header=TRUE, stringsAsFactors=FALSE)

quant_csTSS_csRNA <- suppressMessages(readr::read_tsv(trimws(QUANT_csTSS_csRNA)))
quant_csTSS_input <- suppressMessages(readr::read_tsv(trimws(QUANT_inTSS_csRNA)))
quant_inTSS_csRNA <- suppressMessages(readr::read_tsv(trimws(QUANT_csTSS_input)))
quant_inTSS_input <- suppressMessages(readr::read_tsv(trimws(QUANT_inTSS_input)))

f_mirna <- trimws(MIRNAS)
f_trna <- trimws(TRNAS)

mirna <- if (nchar(f_mirna)) import(f_mirna) else NULL
trna <- if (nchar(f_trna)) import(f_trna) else NULL

stats_cs[['NuclearReads']] <- with(stats_cs, TotalReads - (PtReads + MtReads))
stats_in[['NuclearReads']] <- with(stats_in, TotalReads - (PtReads + MtReads))
stats_cs[['PctNuclear']] <- 100 * (stats_cs[['NuclearReads']] / stats_cs[['TotalReads']])
stats_in[['PctNuclear']] <- 100 * (stats_in[['NuclearReads']] / stats_in[['TotalReads']])

quant_csTSS_csRNA_m <- quant_csTSS_csRNA[, 20:ncol(quant_csTSS_csRNA)]
colnames(quant_csTSS_csRNA_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_csTSS_csRNA_m)))

quant_csTSS_input_m <- quant_csTSS_input[, 20:ncol(quant_csTSS_input)]
colnames(quant_csTSS_input_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_csTSS_input_m)))

quant_inTSS_csRNA_m <- quant_inTSS_csRNA[, 20:ncol(quant_inTSS_csRNA)]
colnames(quant_inTSS_csRNA_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_inTSS_csRNA_m)))

quant_inTSS_input_m <- quant_inTSS_input[, 20:ncol(quant_inTSS_input)]
colnames(quant_inTSS_input_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_inTSS_input_m)))

stats_cs[['csRiP']] <- colSums(quant_csTSS_csRNA_m)[stats_cs[['Sample']]]
stats_cs[['sRiP']] <- colSums(quant_inTSS_csRNA_m)[stats_cs[['Sample']]]

stats_in[['csRiP']] <- colSums(quant_csTSS_input_m)[stats_in[['Sample']]]
stats_in[['sRiP']] <- colSums(quant_inTSS_input_m)[stats_in[['Sample']]]

stats_cs[['csFRiP']] <- stats_cs[['csRiP']] / stats_cs[['NuclearReads']]
stats_cs[['sFRiP']] <- stats_cs[['sRiP']] / stats_cs[['NuclearReads']]
stats_in[['csFRiP']] <- stats_in[['csRiP']] / stats_in[['NuclearReads']]
stats_in[['sFRiP']] <- stats_in[['sRiP']] / stats_in[['NuclearReads']]

stats_cs[['csRNACappedPct']] <- 100 * stats_cs[['csFRiP']]
stats_cs[['InputCappedPct']] <- 100 * stats_in[['csFRiP']]
stats_cs[['csEnrichment']] <- stats_cs[['csFRiP']] / stats_in[['csFRiP']]
stats_cs[['sDepletion']] <- 1 / (stats_cs[['sFRiP']] / stats_in[['sFRiP']])

if (!is.null(trna)) {
    in_pretrna <- mcols(tss_in)[['name']][overlapsAny(tss_in, trna)]
    stats_cs[['PretRNA']] <- colSums(quant_inTSS_csRNA_m[quant_inTSS_csRNA[[1]] %in% in_pretrna, ])[stats_cs[['Sample']]]
    stats_in[['PretRNA']] <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_pretrna, ])[stats_in[['Sample']]]
    stats_cs[['PretRNAPct']] <- with(stats_cs, 100 * (PretRNA / (PretRNA + csRiP)))
    stats_in[['PretRNAPct']] <- with(stats_in, 100 * (PretRNA / (PretRNA + csRiP)))
    stats_cs[['PhosEfficiency']] <- 1 / (stats_cs[['PretRNAPct']] / stats_in[['PretRNAPct']])
}

if (!is.null(mirna)) {
    in_mirna <- mcols(tss_in)[['name']][overlapsAny(tss_in, mirna)]
    stats_cs[['miRNA']] <- colSums(quant_inTSS_csRNA_m[quant_inTSS_csRNA[[1]] %in% in_mirna, ])[stats_cs[['Sample']]]
    stats_in[['miRNA']] <- colSums(quant_inTSS_input_m[quant_inTSS_input[[1]] %in% in_mirna, ])[stats_in[['Sample']]]
    stats_cs[['miRNADepletion']] <- 1 / ((stats_cs[['miRNA']] / stats_cs[['NuclearReads']]) / (stats_in[['miRNA']] / stats_in[['NuclearReads']]))
}

stats_cs[['Status']] <- ifelse(stats_cs[['csFRiP']] > snakemake@config[["qc"]][["min_cs_frip"]] & stats_cs[['PctNuclear']] > snakemake@config[["qc"]][["min_pct_nuclear"]], 'Ok', 'FAIL')

readr::write_tsv(stats_cs, OUT_CS)
readr::write_tsv(stats_in, OUT_IN)

cat('Printing summary stats:\n')
print(cbind(stats_cs[, c('Sample', 'Status', 'PctNuclear', 'csFRiP')], csFRiP_in = stats_in[['csFRiP']]))
