TMM normalisation factors
--------------------------

edgeR ``TMMwsp`` normalisation factors and derived per-million multipliers
for each csRNA sample.

- **norm.factors** — edgeR TMM factor
- **lib.size** — total library size (raw tag count)
- **final.factors** — ``norm.factors × lib.size`` (the effective library size)
- **mult.per.million** — ``(final.factors / 1 000 000)⁻¹``; multiply raw bedGraph
  scores by this value to obtain RPM-normalised bigWig tracks
