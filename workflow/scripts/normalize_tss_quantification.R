suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(edgeR))))
suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

tss <- import(snakemake@input[["bed"]])

quant <- as.data.frame(suppressWarnings(suppressMessages(readr::read_tsv(trimws(snakemake@input[["quant"]])))))
tag_cols <- grep(" Tag Count", colnames(quant))
quant_m <- as.matrix(quant[, tag_cols])
colnames(quant_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_m)))
rownames(quant_m) <- quant[[1]]
readr::write_tsv(as.data.frame(cbind(TSS = rownames(quant_m), quant_m)),
    snakemake@output[["raw"]])

y <- DGEList(quant_m, group = gsub('_csrna[0-9]+$', '', colnames(quant_m)))
y <- calcNormFactors(y, method = 'TMMwsp')

nf <- y[['samples']]
nf[['sample']] <- rownames(nf)
nf[['final.factors']] <- nf[['norm.factors']] * nf[['lib.size']]
nf[['mult.per.million']] <- (nf[['final.factors']] / 1000000)^-1
readr::write_tsv(nf, snakemake@output[["normfactors"]])

y_cpm <- cpm(y)
readr::write_tsv(as.data.frame(cbind(TSS = rownames(y_cpm), y_cpm)),
    snakemake@output[["norm"]])
