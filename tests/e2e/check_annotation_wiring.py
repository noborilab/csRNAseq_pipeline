#!/usr/bin/env python3
"""
Check that the optional GTF is actually reaching HOMER and changing what it does.

    python tests/e2e/check_annotation_wiring.py <e2e output dir>

-gtf is easy to wire up in a way that looks fine and does nothing: the pipeline would still
run, the outputs would still appear, and the only sign it never arrived would be a
placeholder line buried in a per-sample log. So assert on HOMER's own report rather than on
the option being present in the config.
"""

import glob
import os
import sys

def fail(msg):
    print(f"  FAIL: {msg}", file=sys.stderr)
    sys.exit(1)

def ok(msg):
    print(f"  ok  {msg}")

def main(outdir):
    stats = sorted(glob.glob(os.path.join(outdir, "tss", "*_csrna*.stats.txt")))
    if not stats:
        fail(f"no csRNA *.stats.txt under {outdir}/tss")
    text = open(stats[0]).read()
    logs = sorted(glob.glob(os.path.join(outdir, "qc", "*_csrna*.find_tss.log")))
    log = open(logs[0]).read() if logs else ""

    # Do not assert on "Skipping TSS assignment": that comes from the annotatePeaks calls
    # inside findcsRNATSS, which want a full HOMER genome directory for their own
    # TSS-distance annotation, and it appears whether or not -gtf was given. The signal that
    # -gtf actually took effect is that HOMER could build true and false positive sets and
    # therefore chose the enrichment threshold from the data instead of falling back to
    # -defaultLog2Fold.
    if "Custom annotation GTF file:" not in log:
        fail("HOMER never logged reading a custom annotation GTF, so -gtf did not arrive")
    ok("HOMER read the custom annotation GTF")

    tp = [int(l.split(":")[1]) for l in text.splitlines()
          if l.strip().startswith("Total TP (tss) regions:")]
    fp = [int(l.split(":")[1]) for l in text.splitlines()
          if l.strip().startswith("Total FP (exon) regions:")]
    if not tp or max(tp) == 0:
        fail(f"HOMER found no annotated-TSS true positives ({tp}), so the annotation-based "
             "threshold could not be computed and it fell back to -defaultLog2Fold")
    if not fp or max(fp) == 0:
        fail(f"HOMER found no exonic false positives ({fp}); the threshold calculation "
             "needs both sets")
    ok(f"HOMER built its TP/FP sets from the annotation (TP {max(tp)}, FP {max(fp)})")

    if "Maximum CDF difference: -10000000000" in text.split("vs. RNA-seq")[0]:
        fail("input-side Maximum CDF difference is the no-op sentinel, so the "
             "annotation-based threshold selection did not run")
    ok("HOMER chose the input enrichment threshold from the data, not the fallback")

    if "cmd = " not in text:
        fail("no cmd line in the stats file")
    cmd = next(l for l in text.splitlines() if l.startswith("cmd = "))
    if "-gtf" not in cmd:
        fail(f"-gtf missing from the HOMER command line: {cmd}")
    ok("HOMER was invoked with -gtf")

    print()
    print("All annotation-wiring assertions passed.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1])
