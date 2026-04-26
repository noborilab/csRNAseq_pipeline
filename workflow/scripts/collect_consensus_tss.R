suppressMessages(suppressPackageStartupMessages(library(rtracklayer)))

input <- data.frame(row.names = NULL,
  id = NA,
  sample_name = strsplit(snakemake@params[["sample_names"]], " ", fixed = TRUE)[[1]],
  sample_type = strsplit(snakemake@params[["sample_types"]], " ", fixed = TRUE)[[1]],
  replicate = strsplit(snakemake@params[["replicates"]], " ", fixed = TRUE)[[1]]
)
input$id <- paste0(input$sample_name, "_", input$sample_type, input$replicate)
input <- input[input$sample_type == "csrna", ]

nreps <- snakemake@config[["filtering"]][["tss_min_reps"]]
if (nreps > 1) {
  samples <- unique(input$sample_name)
} else {
  samples <- input$id
}
tss_all <- vector('list', length(samples))
sampleTSSs <- function(s, fn, n) {
    if (length(fn) == 1) {
        tss <- sort(reduce(sort(import(fn))))
    } else {
        tss <- lapply(fn, import)
        for (i in seq_along(tss)) {
            mcols(tss[[i]])[['overlaps']] <- 0
            for (j in seq_along(tss)[-i]) {
                mcols(tss[[i]])[['overlaps']] <- mcols(tss[[i]])[['overlaps']] + as.integer(overlapsAny(tss[[i]], tss[[j]], ignore.strand=FALSE))
            }
        }
        for (i in seq_along(tss)) {
            tss[[i]] <- tss[[i]][mcols(tss[[i]])[['overlaps']] >= (n - 1)]
        }
        tss <- do.call(c, tss)
        tss <- sort(reduce(sort(tss)))
    }
    cat('    ', s, ':\t', length(tss), ' TSSs\n', sep = '')
    tss
}
for (i in seq_along(samples)) {
    if (nreps > 1) {
      input_i <- input[input$sample_name %in% samples[i], ]
    } else {
      input_i <- input[input$id %in% samples[i], ]
    }
    samples_i <- paste0(snakemake@params[["bed_dir"]], '/', input_i$id, '.tss.bed')
    tss_all[[i]] <- sampleTSSs(samples[i], samples_i, nreps)
}
cs <- read.table(snakemake@params[["chrom_sizes"]], header = FALSE)
valid_chroms <- as.character(cs$V1)

tss <- do.call(c, tss_all)
tss <- sort(reduce(sort(tss)))
tss <- tss[as.character(seqnames(tss)) %in% valid_chroms, ]
seqlevels(tss, pruning.mode = "coarse") <- intersect(valid_chroms, seqlevels(tss))
tss <- sort(tss)
tss[width(tss) < 150] <- resize(tss[width(tss) < 150], 150, 'center')
tssOvs <- findOverlaps(tss)
tssOvs <- tssOvs[queryHits(tssOvs) != subjectHits(tssOvs)]
if (length(tssOvs)) {
  tssOvs <- tssOvs[seq(1, length(tssOvs), by = 2)]
  if (length(tssOvs)) {
    cat('Adjusting overlapping TSSs ...\n')
    fix_ov_tss <- function(x) {
      wOv <- end(x[1]) - start(x[2])
      mvR <- wOv %/% 2
      mvL <- -(mvR + wOv %% 2)
      x[1] <- shift(x[1], mvL)
      x[2] <- shift(x[2], mvR + 1)
      x
    }
    for (i in seq_len(length(tssOvs))) {
      tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])] <-
        fix_ov_tss(tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])])
    }
  }
}
mirnas <- snakemake@params[["mirnas"]]
if (nzchar(mirnas)) {
    mi <- import(mirnas)
    n_before <- length(tss)
    tss <- subsetByOverlaps(tss, mi, ignore.strand = TRUE, invert = TRUE)
    cat('Removed ', n_before - length(tss), ' miRNA-overlapping TSSs\n', sep = '')
}
trnas <- snakemake@params[["trnas"]]
if (nzchar(trnas)) {
    tr <- import(trnas)
    n_before <- length(tss)
    tss <- subsetByOverlaps(tss, tr, ignore.strand = TRUE, invert = TRUE)
    cat('Removed ', n_before - length(tss), ' pre-tRNA-overlapping TSSs\n', sep = '')
}
mcols(tss)[['name']] <- paste0('TSS_', 1:length(tss));
cat('Final TSS count: ', length(tss), '\n', sep = '')
export.bed(tss, snakemake@output[[1]])
