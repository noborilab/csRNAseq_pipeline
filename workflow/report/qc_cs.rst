csRNA library QC summary
-------------------------

Per-sample quality metrics for csRNA (capped) libraries. Samples with
``Status = FAIL`` did not meet one or both of the configured thresholds
(``min_cs_frip`` and ``min_pct_nuclear``).

Key columns:

- **csFRiP** — fraction of nuclear reads falling within csRNA-called TSSs; should be > 0.9 for a good library
- **PctNuclear** — percentage of reads originating from nuclear chromosomes
- **csEnrichment** — ratio of csRNA FRiP to input FRiP; measures how selectively the capped fraction is enriched at TSSs
- **sDepletion** — inverse ratio of small-RNA signal; measures depletion of uncapped reads
- **PretRNAPct** / **PhosEfficiency** — pre-tRNA contamination and 5′-phosphate removal efficiency (if ``trnas`` BED is configured)
- **miRNADepletion** — miRNA depletion ratio as an independent phosphorylation efficiency check (if ``mirnas`` BED is configured)
