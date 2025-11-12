## Running Snakemake

Recommended `--cores`: 6

## Required config entries

### `sample_table`

This must be a TSV consisting of at least three columns, titled: `sample_name`, `sample_type`, and `read_r1`. Additional optional columns include: `replicate` and `read_r2`. By default the pipeline assumes a file called "samples.tsv" is present in the current working directory if this parameter is not configured.

- `sample_name`: This must be the name of the sample, not including any replicate number or library type information. This is so the pipeline can match corresponding csRNA-seq and input libraries based on their matching `sample_name` and `replicate` columns. During read quantification, this information is also used for normalization; all samples with a matching `sample_name` are considered to be of the same group.
- `replicate`: Should be an integer representing the replicate number. Make sure the corresponding csRNA-seq and input libraries have identical values. If this column is not present, all samples are assumed to be present without replicates (and consequently, all sample names outside of their `sample_type` combinations must be unique).
- `sample_type`: Can be one of two values representing the type of library: `csrna` and `input`.
- `read_r1`: Paths to read 1 FASTQ files. If a sample has multiple files, separate them using commas.
- `read_r2`: Paths to read 2 FASTQ files if paired-end. Multiple files can be separated using commas. This column is optional. If a sample list contains a mix of single-end and paired-end data, then this column can simply be left blank for the single-end samples (but be sure to keep the tab character). An example is shown below:

```
sample_name	sample_type	read_r1	read_r2
sample1	csrna	sample1_1.csrna.r1.fq.gz,sample1_2.csrna.r1.fq.gz	
sample1	input	sample1.input.r1.fq.gz	sample1.input.r2.fq.gz
```

### `chrom_sizes`

A regular `chrom.sizes` file for the organism in question. If one is not readily available, a fasta index can instead be used, which can be created using the `samtools index` command (using the genome fasta as input). A convenient way to obtain the genome fasta is to simply use the one provided by HOMER.

### `program / genome_index`

The expected genome index path by whichever alignment tool is being used. If not present, but a genome fasta has been provided (`genome_fasta`) then the pipeline will attempt to create the index (which will be stored in the `intermediate_dir`). If you wish for an index to be created for you, please do not use the same file path for both `genome_index` and `genome_fasta` (even if they can be the same, e.g. in the case of the bwa index) since the pipeline will then assume the index is already present. The consequence of needing to use different paths means the genome fasta will be duplicated unnecessarily for a bwa index, but if this is a concern then simply create the index manually beforehand (the `genome_fasta` is otherwise not needed).

### `program / homer / genome`

Some steps in the pipeline require the use of programs included in HOMER. To use HOMER, it must be configured to have available the appropriate genome. The name of the genome as it is called by HOMER must be provided in the config file.


