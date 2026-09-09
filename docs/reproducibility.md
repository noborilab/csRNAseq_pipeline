# Reproducing an analysis

## Software and reference snapshots

The workflow uses the tracked `workflow/envs/full_env.yaml` for `--use-conda`, resolved
relative to `workflow/Snakefile`. `minimal_env.yaml` is for manually activated bwa-only
runs. Snakemake 9.13.3 automatic conda deployment requires the `conda` executable;
with a micromamba-only installation, activate or `micromamba run` the prepared environment
and omit `--use-conda`. Core analysis package versions are pinned in both specs. Optional tools and
transitive dependencies are not completely locked; a successful solve is not evidence
that every platform has been tested. HOMER requires manual installation on macOS, and
`bfqutils` and `bam2td` require manual installation on all platforms.

After testing an environment, archive its exact conda package build list:

```bash
micromamba list -n csRNAseq --explicit > environment.explicit.txt
# Recreate on the same OS/architecture:
micromamba create -n csRNAseq_replay -f environment.explicit.txt
```

Also preserve manually installed tools or their source revisions and build instructions.
A conda export does not include them. Keep any extra R packages installed outside conda,
and confirm the replay with the synthetic suite. For automated per-rule environment
locking, Snakemake supports platform-specific explicit pin files; see its
[deployment documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

`tsl/singularity.def` now uses a public AlmaLinux base, fixed revisions for bfqutils,
bam2td and bwa, a bedtools release, a specified HOMER/genome version, and a fixed
Bioconductor release. It copies the local workflow instead of cloning a moving branch.
Build from a clean source snapshot at the repository root so that local config/data
files are not bundled inadvertently. This is a bwa-focused Linux recipe, not an image
containing all optional aligners.

OS repositories, the installer and some R/transitive packages can still change. The
recipe is **not a complete immutable build lock**. Archive the tested SIF image and its
SHA256 alongside the analysis rather than rebuilding it and assuming equivalence:

```bash
apptainer build csrnaseq.sif tsl/singularity.def
sha256sum csrnaseq.sif > csrnaseq.sif.sha256
```

The revised recipe has not been built on the macOS development host; validate it on
Linux before distribution. Exact version pins do not replace that check.

## What `run_info.txt` records

The manifest depends on the final count/QC outputs and both bigWig strands. It includes:

- Resolved configuration, pipeline version, commit and tracked-change status.
- SHA256 fingerprints of the sample table, provided FASTA, chromosome sizes, GTF and
  QC BED annotations, and all discovered aligner index components.
- For an installed HOMER genome name, its installation `config.txt` and files under the
  configured genome directory. Failure to resolve an alias is recorded explicitly.
- Workflow source and executable fingerprints, exact installed Python package versions,
  conda package build identifiers, and R package versions/session information.

Hashes identify content; they do not archive that content. Keep the referenced files,
FASTQs, sample table and software snapshot with the manifest. FASTQs are not hashed by
this rule; retain their sequencing-provider checksums separately. Reference hashing can
add I/O for large genomes and indexes.

Use immutable reference paths. Replacing an existing FASTA or annotation in place is
not a supported way to switch genome versions: the manifest can describe current files
without proving that every existing intermediate was regenerated. For a changed
reference, use a new index path and a fresh output directory. Existing index sentinels
are now alignment inputs, but the workflow does not validate every index component's
completeness or infer which FASTA built an externally supplied index.

## Normalization and bigWig scale

The count matrix uses edgeR **TMMwsp**, the singleton-pairing variant of TMM. Library
sizes are recomputed from the retained TSS count matrix. For sample `s`, every raw
bedGraph value is multiplied by:

```
1e6 / (TMMwsp_factor[s] * sum_of_counts_in_retained_TSSs[s])
```

The same factor scales both strands; negative-strand values are made negative for
display. The `.rpm.*.bw` filenames are retained for compatibility. Their values are
relative to the TMMwsp-adjusted retained-TSS library size, not total raw, mapped or
nuclear read depth. Changing catalogue membership can therefore rescale unchanged
loci. Optional masking removes input-called regions that do not overlap final TSSs;
it does not refit the multiplier after masking.

These tracks support relative comparisons under normalization assumptions. They do not
establish absolute global changes in transcription or cap capture. For those questions,
use suitable external controls or a validated stable reference and account for changes
in RNA composition. See the [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).
