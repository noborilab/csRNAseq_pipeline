Input library QC summary
-------------------------

Per-sample quality metrics for input (total RNA) libraries. These metrics
characterise the background signal used by HOMER for TSS calling.

Key columns: **csFRiP** (fraction of input reads in csRNA TSSs),
**sFRiP** (fraction in input-library TSSs), **PctNuclear**.
A low ``csFRiP`` relative to the matched csRNA library confirms that
the csRNA enrichment step worked correctly.
