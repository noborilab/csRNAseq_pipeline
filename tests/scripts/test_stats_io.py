"""Exercise the stats writer and the production R stats-reading expressions."""
import csv
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class StatsIOTests(unittest.TestCase):
    def test_writer_missing_and_zero_frequency(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tagdir = root / "sample"
            tagdir.mkdir()
            (tagdir / "tagInfo.txt").write_text("genome=test\t1\t30\nchr1\t1\t30\n")
            (root / "sample.trimming.txt").write_text("Processed 40 reads\nOutput 30 reads\n")
            (root / "sample.aln.raw.txt").write_text(
                "PrimaryMapped\t30\nBelowMapq\t0\nFailedOtherFilters\t0\n")
            for strand, count in [("pos", 10), ("neg", 20)]:
                with gzip.open(root / f"sample.raw.{strand}.bedGraph.gz", "wt") as handle:
                    handle.write(f"chr1\t0\t1\t{count}\n")
            output = root / "stats.txt"
            env = dict(os.environ, OUT_PREFIX=str(root), BG_PREFIX=str(root),
                       LOG_PREFIX=str(root), OUT=str(output), ORGANELLE_CHROMS="chrC chrM")
            for content, expected in [("1\t0.25\n", "NA"), ("0\t\n", "NA"),
                                      ("0\t\t0.5\n", "NA"),
                                      ("0\t0\n", "0"), ("0\t0.25\n", "0.25")]:
                with self.subTest(content=content):
                    (tagdir / "tagFreq.txt").write_text(content)
                    subprocess.run(["bash", str(ROOT / "workflow/scripts/gather_stats.sh"),
                                    "sample"], env=env, check=True)
                    with output.open() as handle:
                        row, = csv.DictReader(handle, delimiter="\t")
                    self.assertEqual(row["Freq1A"], expected)
                    self.assertEqual(row["PosReads"], "10")
                    self.assertEqual(row["NegReads"], "20")
                    self.assertTrue(all(value not in (None, "") for value in row.values()))

    def test_both_readers_preserve_columns_and_missing_values(self):
        # Evaluate the actual reader assignments, without loading unrelated genomic
        # inputs. Test each stats input independently, including legacy empty fields.
        subprocess.run(["Rscript", "--vanilla", "-e", r'''
            setClass("StatsIOTest", slots = c(input = "list"))
            path <- tempfile(fileext = ".tsv")
            writeLines(c("Sample\tOrganelleReads\tFreq1A\tPosReads\tNegReads",
                         "s1\t0\t\t10\t20", "s2\t0\tNA\t11\t21",
                         "s3\t0\tna\t12\t22", "s4\t0\t0\t13\t23",
                         "s5\t0\t0.25\t14\t24"), path)
            for (script in c("qc_initial_tss.R", "qc_final_tss.R")) {
                env <- new.env()
                env$STATS_CS <- env$STATS_IN <- path
                env$snakemake <- new("StatsIOTest", input = list(stats_cs=path, stats_in=path))
                exprs <- parse(file.path("workflow/scripts", script))
                readers <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) && is.symbol(x[[2]]) &&
                    as.character(x[[2]]) %in% c("stats_cs", "stats_in"), as.list(exprs))
                stopifnot(length(readers) == 2L)
                for (expr in readers) eval(expr, env)
                for (tab in list(env$stats_cs, env$stats_in)) {
                    stopifnot(nrow(tab) == 5L, ncol(tab) == 5L,
                              identical(tab$Sample, paste0("s", 1:5)),
                              all(is.na(tab$Freq1A[1:3])),
                              identical(tab$Freq1A[4:5], c(0, 0.25)),
                              identical(tab$PosReads, 10:14),
                              identical(tab$NegReads, 20:24))
                }
            }
            unlink(path)
        '''], cwd=ROOT, check=True)


if __name__ == "__main__":
    unittest.main()
