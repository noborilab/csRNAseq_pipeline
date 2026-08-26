suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(rtracklayer))))

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

# The consensus set is a union, so a library that failed QC does not merely have poor
# numbers of its own: the spurious clusters it calls enter the shared set and everything
# downstream quantifies them. Report that either way, and drop those libraries when
# filtering.exclude_failed_from_consensus is on (off by default, which is the behaviour the
# pipeline has always had).
qc <- read.delim(trimws(snakemake@input[["qc_cs"]]), stringsAsFactors = FALSE)
failed <- qc$Sample[!is.na(qc$Status) & qc$Status == "FAIL"]
failed <- intersect(failed, input$id)
exclude_failed <- isTRUE(snakemake@params[["exclude_failed"]])
if (length(failed)) {
  if (exclude_failed) {
    cat("Excluding ", length(failed), " csRNA librar", if (length(failed) == 1) "y" else "ies",
        " that failed QC from the consensus: ", paste(failed, collapse = ", "), "\n", sep = "")
    input <- input[!input$id %in% failed, ]
    if (!nrow(input)) {
      stop("every csRNA library failed QC, so the consensus set would be empty. ",
           "Loosen qc/min_cs_frip and qc/min_pct_nuclear, or set ",
           "filtering.exclude_failed_from_consensus to false.")
    }
  } else {
    cat("WARNING: ", length(failed), " csRNA librar", if (length(failed) == 1) "y" else "ies",
        " failed QC and still contribute to the consensus union: ",
        paste(failed, collapse = ", "),
        "\n         Their spurious clusters enter the shared TSS set. Set ",
        "filtering.exclude_failed_from_consensus to true to drop them.\n", sep = "")
  }
}

# From params, not snakemake@config: the rule passes min_reps as a param already, and the
# config copy is the one Snakemake cannot see changing.
nreps <- as.integer(snakemake@params[["min_reps"]])
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

# Widening to 150 bp, and the shifting below, can push a cluster past either end of its
# chromosome. Slide such clusters back inside instead of trimming them, so every cluster
# keeps its width; a chromosome shorter than the cluster itself is left alone and reported.
chrom_len <- setNames(as.numeric(cs$V2), as.character(cs$V1))
nudge_into_chrom <- function(x) {
  len <- chrom_len[as.character(seqnames(x))]
  off <- ifelse(start(x) < 1, 1 - start(x),
         ifelse(!is.na(len) & end(x) > len, len - end(x), 0))
  off[is.na(off)] <- 0
  ok <- is.na(len) | width(x) <= len
  if (any(off != 0 & ok)) {
    cat('Sliding ', sum(off != 0 & ok), ' cluster(s) back inside the chromosome bounds\n', sep = '')
    x[off != 0 & ok] <- shift(x[off != 0 & ok], off[off != 0 & ok])
  }
  if (any(!ok)) {
    cat('WARNING: ', sum(!ok), ' cluster(s) are wider than their chromosome and left as they are\n', sep = '')
  }
  x
}
tss <- nudge_into_chrom(tss)

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
  # Every overlap is reported twice, as (i,j) and (j,i). Keep one of each pair by index
  # order: the previous "every other hit" form relied on the two members of a pair landing
  # adjacent in the Hits object, which holds for equal-width sorted ranges but is not
  # something findOverlaps promises.
  tssOvs <- tssOvs[queryHits(tssOvs) < subjectHits(tssOvs)]
  if (!length(tssOvs)) break
  cat('Adjusting ', length(tssOvs), ' overlapping TSS pair(s) (iteration ', iter + 1, ') ...\n', sep = '')
  for (i in seq_len(length(tssOvs))) {
    tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])] <-
      fix_ov_tss(tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])])
  }
  tss <- nudge_into_chrom(tss)
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
