# Library complexity, from the tag directory rather than the bedGraphs.
#
# Complexity here is the number of distinct molecules a library holds, which
# without UMIs has to be approximated by distinct fragments. The unit is the
# (position, strand, length) triple HOMER already collapses each tag directory
# into, so both ends of a read define a fragment. That matters: measured on this
# panel, 5' position alone treats between a half and three quarters of the
# distinct fragments as duplicates, because an in-TSS position carries a median of
# about 3.7 distinct lengths. Two reads sharing both ends are far more likely to
# be copies of one molecule, which is the definition Picard applies to paired-end
# data.
#
# Restricted to the final TSS set, and that restriction is not cosmetic. Measured
# genome-wide, the two libraries in the panel that failed QC outright ranked first
# and second on complexity, because a library made of scattered background has more
# distinct positions than a good one. Restricted to called TSSs they fall to the
# middle and bottom. Organelle chromosomes are dropped for the same reason csFRiP
# drops them.
#
# Three numbers come out:
#
#   TSSDistinctFrag  distinct fragments observed inside the TSS set. Depends on how
#                    deeply the library was sequenced, so it is the denominator for
#                    the other two rather than a metric to compare across libraries.
#   TSSChao1         Chao1 extrapolation of the fragments the library could yield at
#                    infinite depth, from the singleton and doubleton counts. Needs
#                    no subsampling, which is why it is here: the shallowest library
#                    on the panel had 2.5M in-TSS reads, so a fixed-depth subsample
#                    would either discard 90% of a good library's data or return NA
#                    for the libraries most worth judging. Among libraries that pass
#                    the other gates it tracks a depth-matched count at rho = 0.88.
#   TSSSaturation    Good's coverage, 1 - singletons/reads: the estimated chance that
#                    the next read comes from a fragment already seen. This is the
#                    one to read when deciding whether more sequencing would pay.
#                    On the RNA input series it ran 0.936 at 50 ug and 0.994 at 1 ug.
#
# Neither TSSChao1 nor TSSSaturation should be gated on alone. Chao1 still inflates
# for background-heavy libraries, and Saturation runs at rho = 0.77 against csFRiP
# among passing libraries, so it partly restates purity. Read them beside csFRiP.

suppressWarnings(suppressMessages(suppressPackageStartupMessages({
    library(data.table)
    library(rtracklayer)
    library(GenomicRanges)
})))

sample_id <- snakemake@wildcards[["sample_cs"]]
tagdir    <- trimws(snakemake@input[["tagdir"]])
organelles <- trimws(unlist(strsplit(trimws(snakemake@params[["organelle_chroms"]]), "[[:space:]]+")))
organelles <- organelles[nzchar(organelles)]

out <- data.frame(Sample = sample_id, TSSDistinctFrag = NA_integer_,
                  TSSChao1 = NA_integer_, TSSSaturation = NA_real_)

tags <- list.files(tagdir, pattern = "\\.tags\\.tsv$", full.names = TRUE)
tags <- tags[!sub("\\.tags\\.tsv$", "", basename(tags)) %in% organelles]

if (length(tags)) {
    # HOMER writes one row per (position, strand, length) with a count, and a
    # leading empty field. Columns 2-6 are chr, position, strand (0 = +, 1 = -),
    # count and length.
    x <- rbindlist(lapply(tags, function(f)
        fread(f, header = FALSE, select = 2:6,
              col.names = c("chr", "pos", "strand", "count", "len"),
              colClasses = list(character = 2, integer = c(3, 4),
                                numeric = 5, integer = 6))))

    tss <- import(trimws(snakemake@input[["tss"]]))
    if (length(tss) && nrow(x)) {
        gr <- GRanges(x[["chr"]], IRanges(x[["pos"]], width = 1L),
                      strand = ifelse(x[["strand"]] == 0L, "+", "-"))
        # ignore.strand because the TSS set is stranded and a fragment counts as
        # being at a TSS on either strand, matching how csRiP is computed.
        keep <- overlapsAny(gr, tss, ignore.strand = TRUE)
        v <- as.integer(round(x[["count"]][keep]))
        v <- v[v > 0]

        if (length(v)) {
            N  <- sum(as.numeric(v))
            F  <- length(v)
            f1 <- sum(v == 1L)
            f2 <- sum(v == 2L)
            # Chao1, with the f2 = 0 correction. Both branches reduce to F when the
            # library has no singletons, which is what an exhausted library looks
            # like and is the right answer there.
            chao <- F + if (f2 > 0) f1^2 / (2 * f2) else f1 * (f1 - 1) / 2
            out[["TSSDistinctFrag"]] <- F
            out[["TSSChao1"]]        <- round(chao)
            out[["TSSSaturation"]]   <- round(1 - f1 / N, 4)
        }
    }
}

readr::write_tsv(out, snakemake@output[[1]])
