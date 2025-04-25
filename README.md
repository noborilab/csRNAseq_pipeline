# csRNA-seq pipeline

Author: Benjamin Jean-Marie Tremblay

## Requirements

* [homer](http://homer.ucsd.edu/homer/index.html) (with the tair10 genome; tested with v5.1)
* [samtools](https://www.htslib.org) (tested with v1.20)
* [bfqutils](https://github.com/noborilab/bfqutils) (tested with v1.0)
* [bwa](https://github.com/lh3/bwa) (tested with v0.7.19)
* [bedtools](https://github.com/arq5x/bedtools2) (tested with v2.31.1)
* [R](https://cran.r-project.org) (with packages: rtracklayer, readr, edgeR; tested with v4.4.1)

## The pipeline

The scripts included in this pipeline allow for going from raw reads to a final set of merged TSSs, alongside quantification data and normalized bigWigs. All scripts always check whether an output exists already and will only proceed unless it is absent, or if the appropriate option is set to override this behaviour. This makes it easy to add samples to an existing project without having to process different batches in separate directories and needing to manually integrate them every time.

### Part 1: `process_reads.sh`

This first script processes the raw reads and outputs HOMER tag directories, which are used by all other downstream scripts. Single-end reads are trimmed with bfqtrimse, and paired-end reads are merged using bfqmerge. Trimmed reads are aligned to a genome using bwa-aln and filtered, finally allowing for the creation of tag directories.

### Part 2: `identify_tss.sh`

This second script goes through all samples and identifies TSSs as well as creating raw bedGraph files. It also creates a few additional outputs used by the next step in the pipeline, quality control. These include sRNA peaks from the input samples, and quantification data from a combined set of all identified TSSs.

### Part 3: `qc_tss.sh`

Using the outputs from the last part of the pipeline, this script calculates key quality control statistics including the percentage of nuclear reads and the fraction of reads in peaks. At the end it will output a final sample quality indicator (OK or FAIL).

### Part 4: `finalize_tss.sh`

After performing quality control, it may be that some samples will need to be discarded before assembling a final set of TSSs and performing quantification. This script will take selected high quality samples and generate a final set of TSSs and quantify them. It will then use this quantification data to create normalized bigWig files.

