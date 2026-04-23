Normalised TSS counts (CPM)
----------------------------

TMM-normalised counts per million (CPM) for each consensus TSS across all
csRNA samples, computed using edgeR ``TMMwsp`` normalisation.

Rows are consensus TSSs (``TSS_N`` identifiers matching ``tss.final.bed``).
Columns are csRNA sample IDs.

These values are suitable for downstream differential expression analysis
or visualisation. For differential testing, use the raw counts in
``tss.final.raw.txt`` together with the edgeR normalisation factors in
``norm_factors.txt``.
