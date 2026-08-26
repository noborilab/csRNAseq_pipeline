# Rewrite HOMER's not-computed stats lines to "na".
#
#   awk -f fix_stats_placeholders.awk <sample>.stats.txt
#
# findcsRNATSS writes these two lines whether or not the input they need was supplied, so
# a placeholder reads as a measurement rather than as "not computed":
#
#   Fraction of stable transcript TSS clusters: 0.00%    with no -rna, which is every run
#                                                        since 0.9.0 dropped that arm
#   Fraction Promoter-Distal TSS clusters: 100.00%       with no annotation, since every
#                                                        cluster then counts as distal
#
# Each is only a placeholder when HOMER had nothing to compute it from, which it reports
# as an empty true-positive set on the relevant side: the first "Total TP (tss) regions"
# line is the input (annotation) side, the second the RNA-seq side. Annotation can come
# from -gtf or from a HOMER genome carrying its own, so the test is HOMER's own report
# rather than which flags the pipeline passed. Both the empty set and the placeholder
# value have to hold before a line is rewritten, so a genuine measurement that happens to
# land on 100.00% or 0.00% is left alone.
#
# The file is held in memory because the TP counts appear below the lines being rewritten.
# It runs to about fifty lines.

{ line[NR] = $0 }

/^[ \t]*Total TP \(tss\) regions:/ {
    n_tp++
    split($0, field, ":")
    tp[n_tp] = field[2] + 0
}

END {
    for (i = 1; i <= NR; i++) {
        l = line[i]
        if (n_tp >= 2 && tp[2] == 0 &&
            l ~ /^Fraction of stable transcript TSS clusters: 0\.00%$/) {
            l = "Fraction of stable transcript TSS clusters: na"
        } else if (n_tp >= 1 && tp[1] == 0 &&
                   l ~ /^Fraction Promoter-Distal TSS clusters: 100\.00%$/) {
            l = "Fraction Promoter-Distal TSS clusters: na"
        }
        print l
    }
}
