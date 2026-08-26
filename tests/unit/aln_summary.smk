"""
Unit test for workflow/scripts/aln_summary.awk

The awk sits in the alignment pipe, so two things have to hold: the stream it passes on
must be byte-identical to the stream it received, and its counts must partition the reads
correctly. The SAM fixture holds one record for each branch (unmapped, secondary,
supplementary, below MAPQ, too short, too many mismatches, no NM tag at all), which the
synthetic e2e libraries cannot exercise.

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/aln_summary.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_aln_summary")


rule test_aln_summary:
    input:
        sam=FIXTURES + "/aln_summary.sam",
    output:
        passthrough=TEST_OUTDIR + "/passthrough.sam",
        summary=TEST_OUTDIR + "/aln.raw.txt",
        unmapped=TEST_OUTDIR + "/unmapped.sample.fastq.gz",
    params:
        awk=workflow.basedir + "/../../workflow/scripts/aln_summary.awk",
    shell:
        # Gates chosen so every branch fires: maxmm=2 fails one record, minlen=15 fails
        # another, mapq=30 fails two, and unmapped_n=2 keeps two of the three unmapped
        # reads so the cap is exercised rather than just the write.
        """
            awk -f {params.awk} -v out={output.summary} \
                -v mapq=30 -v excl=2308 -v incl=0 -v minlen=15 -v maxmm=2 -v nm_tag=NM \
                -v unmapped_fq={output.unmapped} -v unmapped_n=2 \
                {input.sam} > {output.passthrough}
        """
