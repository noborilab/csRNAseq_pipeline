#!/usr/bin/env bash
# Sort a tab-delimited table on stable keys, keeping its single header line first.
#
#   bash sort_table.sh <file> <sort key arg>...
#
# annotatePeaks.pl and findcsRNATSS.pl both emit their rows in whatever order their
# threads finish in, so two runs on the same data produced tables holding identical
# numbers in a different order. Nothing here ever read those tables positionally, so no
# result depended on it, but it did leave a run unreproducible byte for byte, and it left
# cor() adding the same pairs up in a different order each time, which moved the reported
# Pearson correlation in the last bit of a double. Sorting the output here rather than
# asking for a flag upstream keeps HOMER itself untouched, which is the same choice
# fix_stats_placeholders.awk makes.
#
# The separator is fixed to a tab because these tables carry values with spaces in them
# (annotation names, gene descriptions), and sort's default whitespace splitting would
# not line up with the columns the caller is naming. LC_ALL=C because the collation of a
# string key otherwise follows the machine's locale, which would make the output
# reproducible on one host but not between two.
set -euo pipefail
f="$1"
shift
head -n 1 "$f"
tail -n +2 "$f" | LC_ALL=C sort -t$'\t' "$@"
