suppressMessages(suppressPackageStartupMessages(library(rtracklayer)))

tss <- import(snakemake@input[["bed"]])

if (snakemake@params[["mask_srnas"]]) {
    srna <- import(snakemake@input[["srna"]])
    srna <- srna[!overlapsAny(srna, tss, ignore.strand = FALSE)]
    srna_p <- srna[as.character(strand(srna)) == '+']
    srna_n <- srna[as.character(strand(srna)) == '-']
} else {
    srna <- NULL
}

cs <- read.table(snakemake@params[["cs"]], header = FALSE)
cs <- structure(cs$V2, names = as.character(cs$V1))

nf <- as.data.frame(suppressMessages(readr::read_tsv(snakemake@input[["normfactors"]])))
rownames(nf) <- nf$sample

samplename <- basename(snakemake@output[["bw_p"]])
samplename <- gsub(".rpm.pos.bw$", "", samplename)

bg_p <- suppressWarnings(import(snakemake@input[["bg_p"]]))
bg_n <- suppressWarnings(import(snakemake@input[["bg_n"]]))
bg_p <- bg_p[as.character(seqnames(bg_p)) %in% as.character(1:5), ]
bg_n <- bg_n[as.character(seqnames(bg_n)) %in% as.character(1:5), ]
mcols(bg_p)[['score']] <- mcols(bg_p)[['score']] * nf[samplename, 'mult.per.million']
mcols(bg_n)[['score']] <- -abs(mcols(bg_n)[['score']]) * nf[samplename, 'mult.per.million']
seqlevels(bg_p) <- names(cs)
seqlevels(bg_n) <- names(cs)
seqlengths(bg_p) <- cs
seqlengths(bg_n) <- cs
if (!is.null(srna)) {
    bg_p <- bg_p[!overlapsAny(bg_p, srna_p)]
    bg_n <- bg_n[!overlapsAny(bg_n, srna_n)]
}
export.bw(bg_p, snakemake@output[["bw_p"]])
export.bw(bg_n, snakemake@output[["bw_n"]])
