suppressMessages(suppressPackageStartupMessages(library(rtracklayer)))

# Use pre-cleaned IDs passed from Snakefile directly to avoid ID reconstruction mismatch
all_ids    <- strsplit(snakemake@params[["ids"]],          " ", fixed = TRUE)[[1]]
all_names  <- strsplit(snakemake@params[["sample_names"]], " ", fixed = TRUE)[[1]]
all_types  <- strsplit(snakemake@params[["sample_types"]], " ", fixed = TRUE)[[1]]
all_reps   <- strsplit(snakemake@params[["replicates"]],   " ", fixed = TRUE)[[1]]

input <- data.frame(row.names = NULL,
  id = all_ids, sample_name = all_names,
  sample_type = all_types, replicate = all_reps
)
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

fix_ov_tss <- function(x) {
  wOv <- end(x[1]) - start(x[2])
  mvR <- wOv %/% 2
  mvL <- -(mvR + wOv %% 2)
  x[1] <- shift(x[1], mvL)
  x[2] <- shift(x[2], mvR + 1)
  x
}
max_iter <- 100
iter <- 0
while (iter < max_iter) {
  tssOvs <- findOverlaps(tss, ignore.strand = FALSE)
  tssOvs <- tssOvs[queryHits(tssOvs) != subjectHits(tssOvs)]
  if (!length(tssOvs)) break
  tssOvs <- tssOvs[seq(1, length(tssOvs), by = 2)]
  cat('Adjusting ', length(tssOvs), ' overlapping TSS pair(s) (iteration ', iter + 1, ') ...\n', sep = '')
  for (i in seq_len(length(tssOvs))) {
    tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])] <-
      fix_ov_tss(tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])])
  }
  iter <- iter + 1
}
if (iter == max_iter) {
  stop("ERROR: TSS overlap resolution did not converge after ", max_iter, " iterations")
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
mcols(tss)[['name']] <- paste0('TSS_', seq_along(tss))
cat('Final TSS count: ', length(tss), '\n', sep = '')
export.bed(tss, snakemake@output[[1]])
