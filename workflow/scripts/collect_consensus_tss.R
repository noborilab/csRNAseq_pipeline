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
    message("Excluding ", length(failed), " csRNA librar", if (length(failed) == 1) "y" else "ies",
            " that failed QC from the consensus: ", paste(failed, collapse = ", "))
    input <- input[!input$id %in% failed, ]
    if (!nrow(input)) {
      stop("every csRNA library failed QC, so the consensus set would be empty. ",
           "Loosen qc/min_cs_frip and qc/min_pct_nuclear, or set ",
           "filtering.exclude_failed_from_consensus to false.")
    }
  } else {
    message("WARNING: ", length(failed), " csRNA librar", if (length(failed) == 1) "y" else "ies",
            " failed QC and still contribute to the consensus union: ",
            paste(failed, collapse = ", "),
            "\n         Their calls enter the shared TSS set; a QC flag is not proof that every call is false. Set ",
            "filtering.exclude_failed_from_consensus to true to drop them.")
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
    if (length(fn) < n) {
        message("Skipping ", s, ": ", length(fn), " surviving libraries, but ", n, " required")
        return(GRanges())
    }
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
if (!length(tss)) stop("no TSSs meet tss_min_reps after QC exclusion and replicate filtering")
tss <- sort(reduce(sort(tss)))
tss <- tss[as.character(seqnames(tss)) %in% valid_chroms, ]
seqlevels(tss, pruning.mode = "coarse") <- intersect(valid_chroms, seqlevels(tss))
tss <- sort(tss)
tss[width(tss) < 150] <- resize(tss[width(tss) < 150], 150, 'center')

# Place every chromosome/strand jointly. Pairwise shifts can create another overlap,
# and moving a right-hand interval back inside its chromosome can undo the shift.
# Transform starts s_i to z_i = s_i - sum(widths before i). Nonoverlap is then
# exactly the constraint z_1 <= ... <= z_n. Isotonic regression finds the closest
# placement in squared displacement, and clamping its fit enforces both chromosome
# bounds. Rounding a monotone fit preserves nonoverlap and all original widths.
resolve_tss_overlaps <- function(x, chrom_len) {
  if (!length(x)) return(x)
  if (any(!as.character(strand(x)) %in% c("+", "-"))) {
    stop("TSS overlap resolution requires explicit + or - strands")
  }
  groups <- base::split(seq_len(length(x)),
                        paste(as.character(seqnames(x)), as.character(strand(x)), sep = "\t"))
  moved <- 0L
  for (indices in groups) {
    indices <- indices[order(start(x)[indices], end(x)[indices])]
    chromosome <- as.character(seqnames(x)[indices[1]])
    direction <- as.character(strand(x)[indices[1]])
    len <- unname(chrom_len[chromosome])
    if (length(len) != 1L || !is.finite(len) || len < 1 || len != floor(len)) {
      stop("Missing or invalid chromosome length for ", chromosome)
    }
    widths <- as.numeric(width(x)[indices])
    total_width <- sum(widths)
    if (any(widths < 1) || total_width > len) {
      stop("Cannot place ", length(indices), " TSS clusters on ", chromosome,
           " (", direction, "): combined width ", total_width,
           " bp exceeds chromosome length ", len,
           " bp, or a cluster has zero width. Nonoverlapping placement cannot preserve cluster widths.")
    }
    offsets <- c(0, head(cumsum(widths), -1L))
    desired <- as.numeric(start(x)[indices]) - offsets
    fitted <- if (length(indices) == 1L) desired else stats::isoreg(desired)$yf
    fitted <- pmax(1, pmin(len - total_width + 1, floor(fitted + 0.5)))
    new_starts <- fitted + offsets
    displacement <- new_starts - start(x)[indices]
    moved <- moved + sum(displacement != 0)
    x[indices] <- shift(x[indices], displacement)
  }
  cat('Placed TSS clusters without overlaps; shifted ', moved, ' cluster(s)\n', sep = '')
  sort(x)
}
chrom_len <- setNames(as.numeric(cs$V2), as.character(cs$V1))
tss <- resolve_tss_overlaps(tss, chrom_len)

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
