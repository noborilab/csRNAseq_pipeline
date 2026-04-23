Consensus TSS coordinates
--------------------------

BED6 file of the final consensus transcription start sites retained after
requiring detection in at least ``tss_min_reps`` replicates per condition.

Post-processing applied:

- TSSs narrower than 150 bp are symmetrically padded to 150 bp
- Any remaining overlapping TSSs are split by shifting them apart equally
- Organelle chromosomes are excluded (controlled by the ``chrom_sizes`` file)

The ``name`` column contains stable identifiers of the form ``TSS_N``.
