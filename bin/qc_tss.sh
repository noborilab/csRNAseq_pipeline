#!/usr/bin/env bash
#
# Author: Benjamin Jean-Marie Tremblay (benjamin.tremblay@tsl.ac.uk)
# Date created: 10 March 2025
# Date modified: 14 May 2025
#

# From quant files: RiPs
# From quant files + BED files: tRNAs, miRNAs
# Fromt stats files: total reads, nuclear/pt/mt reads, pos1 A freq

set -e

WD=`pwd`

PART1_DIR="."
OUT_DIR="."
TRNAS=
MIRNAS=
CSFRIP=0.7
PCTNUC=90

help() {
    echo "csRNA-seq post-processing part 2: quality control checks"
    echo "Benjamin Jean-Marie Tremblay (2025)"
    echo
    echo "Requirements: R (+rtracklayer, +readr)"
    echo
    echo "Options:"
    echo "-d   Path to dir containing outputs from part 1 (default=$PART1_DIR)"
    echo "-o   Output dir for final QC stats (default=$OUT_FILE)"
    echo "-c   csRNA-seq FRiP cutoff (default=$CSFRIP)"
    echo "-n   Percent nuclear reads cutoff (default=$PCTNUC)"
    echo "-t   Path to BED file containing locations of tRNAs (optional)"
    echo "-m   Path to BED file containing locations of miRNAs (optional)"
    echo "-h   Show this message"
    echo
    echo "If using containers, set the following variable in your env:"
    echo '$RSCRIPT'
    echo '(e.g. export RSCRIPT="singularity exec /path/to/containers/R.img")'
    exit 1
}

while getopts "i:d:o:c:n:t:m:h" opt; do
    case $opt in
        d) PART1_DIR=$OPTARG;;
        o) OUT_FILE=$OPTARG;;
        c) CSFRIP=$OPTARG;;
        n) PCTNUC=$OPTARG;;
        t) TRNAS=$OPTARG;;
        m) MIRNAS=$OPTARG;;
        h) help;;
        ?) help;;
    esac
done

STATS_CS=$PART1_DIR/stats_cs.txt
STATS_IN=$PART1_DIR/stats_in.txt
TSS_CS=$PART1_DIR/all_cs.tss_merged.bed
TSS_IN=$PART1_DIR/all_in.tss_merged.bed
QUANT_csTSS_csRNA=$PART1_DIR/all_cs.tss_merged_quant_cs.txt
QUANT_inTSS_csRNA=$PART1_DIR/all_cs.tss_merged_quant_in.txt
QUANT_csTSS_input=$PART1_DIR/all_in.tss_merged_quant_cs.txt
QUANT_inTSS_input=$PART1_DIR/all_in.tss_merged_quant_in.txt

echo "Aggregating QC stats ..."

$RSCRIPT Rscript -e "
library(rtracklayer)

tss_cs <- import(trimws(\"$TSS_CS\"))
tss_in <- import(trimws(\"$TSS_IN\"))

stats_cs <- read.table(trimws(\"$STATS_CS\"), header=TRUE, stringsAsFactors=FALSE)
stats_in <- read.table(trimws(\"$STATS_IN\"), header=TRUE, stringsAsFactors=FALSE)

quant_csTSS_csRNA <- readr::read_tsv(trimws(\"$QUANT_csTSS_csRNA\"))
quant_csTSS_input <- readr::read_tsv(trimws(\"$QUANT_inTSS_csRNA\"))
quant_inTSS_csRNA <- readr::read_tsv(trimws(\"$QUANT_csTSS_input\"))
quant_inTSS_input <- readr::read_tsv(trimws(\"$QUANT_inTSS_input\"))

f_mirna <- trimws(\"$MIRNAS\")
f_trna <- trimws(\"$TRNAS\")

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

stats_cs[['Status']] <- ifelse(stats_cs[['csFRiP']] > $CSFRIP & stats_cs[['PctNuclear']] > $PCTNUC, 'Ok', 'FAIL')

readr::write_tsv(stats_cs, \"$OUT_DIR/qc_cs.txt\")
readr::write_tsv(stats_in, \"$OUT_DIR/qc_in.txt\")

cat('Printing summary stats:\n')
print(cbind(stats_cs[, c('Sample', 'Status', 'PctNuclear', 'csFRiP')], csFRiP_in = stats_in[['csFRiP']]))
"

echo
echo "All done."

