import pandas as pd
import numpy as np

counts = pd.read_csv(
    "results/expression/gse147986_count_matrix.tsv",
    sep="\t",
    index_col=0
)

library_sizes = counts.sum(axis=0)

cpm = counts.div(
    library_sizes,
    axis=1
) * 1_000_000

log_cpm = np.log2(cpm + 1)

cpm.to_csv(
    "results/normalised_expression/gse147986_cpm.tsv",
    sep="\t"
)

log_cpm.to_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t"
)

print("Counts:", counts.shape)
print("CPM:", cpm.shape)
print("logCPM:", log_cpm.shape)
