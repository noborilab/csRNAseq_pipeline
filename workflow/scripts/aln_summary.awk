# Pass-through SAM counter: prints every line unchanged and, at the end, writes a summary
# of what the alignment filter is about to discard.
#
#   ... | awk -f aln_summary.awk -v out=<file> -v mapq=30 -v excl=2308 -v incl=0 \
#             -v minlen=15 -v maxmm=999 -v nm_tag=NM | samtools view -b ...
#
# build_aln_cmd inserts this between the aligner and samtools view, because that is the
# only place the unfiltered alignment exists: it is a pipe, and the flagstat the pipeline
# keeps afterwards runs on the filtered BAM, so it reports 100% mapped for every library
# and carries no information. Without this, the reads lost between TrimmedReads and
# TotalReads (a median 74.9% per library on one Arabidopsis panel) are unaccounted for,
# and alignment failure cannot be told apart from multi-mapper loss, which have opposite
# implications.
#
# An inline filter rather than `tee >(...)` on purpose: process substitution races with
# shell exit and can truncate the summary, whereas a filter in the pipe cannot.
#
# Bit tests use integer arithmetic rather than and(), which is a gawk extension that mawk
# and BSD awk do not have.
#
# Variables (-v):
#   out          where to write the summary (required)
#   mapq         filtering/alignment_mapq            (samtools view -q)
#   excl         filtering/exclude_flags             (samtools view -F)
#   incl         filtering/include_flags             (samtools view -f)
#   minlen       filtering/min_alignment_length      (compared against the CIGAR
#                reference length, which is what samtools calls rlen)
#   maxmm        filtering/max_mismatch
#   nm_tag       NM for most aligners, nM for STAR
#   star_log     optional STAR Log.final.out. STAR keeps unmapped reads out of
#                Aligned.out.bam entirely, so the stream cannot see them and STAR's own
#                report is the only place its read accounting exists.
#   unmapped_fq  optional gzipped FASTQ to sample unmapped reads into
#   unmapped_n   how many unmapped reads to keep in that sample

function bits_any(flag, mask,   b) {
    # true when flag and mask share a set bit, i.e. samtools -F would drop the record
    for (b = 1; b <= mask; b *= 2)
        if (int(mask / b) % 2 == 1 && int(flag / b) % 2 == 1) return 1
    return 0
}

function bits_all(flag, mask,   b) {
    # true when every bit set in mask is set in flag, i.e. samtools -f would keep it
    for (b = 1; b <= mask; b *= 2)
        if (int(mask / b) % 2 == 1 && int(flag / b) % 2 == 0) return 0
    return 1
}

function ref_len(cigar,   i, c, num, total) {
    # length of the alignment on the reference: the CIGAR ops that consume reference
    total = 0
    num = ""
    for (i = 1; i <= length(cigar); i++) {
        c = substr(cigar, i, 1)
        if (c ~ /[0-9]/) {
            num = num c
            continue
        }
        if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") total += num + 0
        num = ""
    }
    return total
}

function mismatches(   i, prefix) {
    # -1 when the tag is absent, which samtools treats as failing the expression
    prefix = nm_tag ":i:"
    for (i = 12; i <= NF; i++)
        if (substr($i, 1, length(prefix)) == prefix)
            return substr($i, length(prefix) + 1) + 0
    return -1
}

function last_field(s,   n, part) {
    # STAR's log is "  Number of input reads |\t12345"
    n = split(s, part, "\t")
    return part[n] + 0
}

BEGIN {
    if (out == "") {
        print "aln_summary.awk: -v out=<path> is required" > "/dev/stderr"
        bad = 1
        exit 2
    }
    mapq += 0; excl += 0; incl += 0; minlen += 0; maxmm += 0; unmapped_n += 0
    if (unmapped_n > 0 && unmapped_fq != "") fq_cmd = "gzip -c > \"" unmapped_fq "\""
    else unmapped_n = 0
}

/^@/ { print; next }

{
    print
    records++
    flag = $2 + 0
    if (int(flag / 4) % 2) {
        unmapped++
        if (sampled < unmapped_n) {
            print "@" $1 "\n" $10 "\n+\n" $11 | fq_cmd
            sampled++
        }
        next
    }
    if (int(flag / 256) % 2)  { secondary++;     next }
    if (int(flag / 2048) % 2) { supplementary++; next }

    # One primary record per read from here on, so these counts are reads, not alignments,
    # and they partition: PrimaryMapped = BelowMapq + FailedOtherFilters + Passed.
    primary++
    q = $5 + 0
    mapq_hist[q]++
    if (q < mapq)                          { below_mapq++; next }
    if (excl > 0 && bits_any(flag, excl))  { filtered++;   next }
    if (incl > 0 && !bits_all(flag, incl)) { filtered++;   next }
    if (minlen > 0 && ref_len($6) < minlen) { filtered++;  next }
    nm = mismatches()
    if (nm < 0 || nm > maxmm)              { filtered++;   next }
    passed++
}

END {
    if (bad) exit 2
    if (star_log != "") {
        while ((getline line < star_log) > 0)
            if (line ~ /Number of input reads/) reads_seen = last_field(line)
        close(star_log)
        # Everything STAR read but did not put in the BAM: unmapped, too short, or mapped
        # to more loci than outFilterMultimapNmax allows.
        if (reads_seen > 0) unmapped = reads_seen - primary
    } else {
        reads_seen = unmapped + primary
    }
    printf "ReadsSeen\t%.0f\n",          reads_seen    > out
    printf "Records\t%.0f\n",            records       > out
    printf "Unmapped\t%.0f\n",           unmapped      > out
    printf "Secondary\t%.0f\n",          secondary     > out
    printf "Supplementary\t%.0f\n",      supplementary > out
    printf "PrimaryMapped\t%.0f\n",      primary       > out
    printf "BelowMapq\t%.0f\n",          below_mapq    > out
    printf "FailedOtherFilters\t%.0f\n", filtered      > out
    printf "Passed\t%.0f\n",             passed        > out
    printf "MinMapq\t%s\n",              mapq          > out
    printf "ExcludeFlags\t%s\n",         excl          > out
    printf "IncludeFlags\t%s\n",         incl          > out
    printf "MinAlignmentLength\t%s\n",   minlen        > out
    printf "MaxMismatch\t%s\n",          maxmm         > out
    printf "UnmappedSampled\t%.0f\n",    sampled       > out
    for (q = 0; q <= 255; q++)
        if (q in mapq_hist) printf "MAPQ\t%d\t%.0f\n", q, mapq_hist[q] > out
    close(out)
    if (unmapped_n > 0) close(fq_cmd)
}
