import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

expr = pd.read_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t",
    index_col=0
)

groups = pd.read_csv(
    "results/normalised_expression/mnda_groups.tsv",
    sep="\t"
)

with open("metadata/gene_sets/pyhin_interferon_target_genes.txt") as f:
    genes = [line.strip() for line in f if line.strip()]

genes = [g for g in genes if g in expr.index]

samples = groups.sort_values("group")["sample_id"].tolist()

mat = expr.loc[genes, samples]

mat_z = mat.sub(mat.mean(axis=1), axis=0)
mat_z = mat_z.div(mat.std(axis=1), axis=0)

out_dir = Path("figures/interferon_heatmap")
out_dir.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 7))

im = ax.imshow(
    mat_z,
    aspect="auto",
    interpolation="nearest"
)

ax.set_yticks(np.arange(len(genes)))
ax.set_yticklabels(genes)

ax.set_xticks([])
ax.set_xlabel("DLBCL samples ordered by MNDA expression group")
ax.set_title("PYHIN and interferon-associated genes in MNDA-high and MNDA-low DLBCL")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Row z-score of log2(CPM + 1)")

fig.tight_layout()

fig.savefig(
    out_dir / "mnda_high_low_pyhin_interferon_heatmap.png",
    dpi=300
)

plt.close(fig)

print("Finished.")
