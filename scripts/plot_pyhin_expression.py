import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

input_file = "results/normalised_expression/pyhin_log_cpm_by_sample.tsv"

output_dir = Path("figures/pyhin")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

genes = ["IFI16", "AIM2", "MNDA", "PYHIN1"]
groups = ["ABC", "GCB", "Type_3"]

rng = np.random.default_rng(42)

for gene in genes:

    data = [
        df.loc[df["coo"] == group, gene].dropna().values
        for group in groups
    ]

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.boxplot(
        data,
        tick_labels=groups,
        showfliers=False
    )

    for i, values in enumerate(data, start=1):

        jitter = rng.normal(
            loc=0,
            scale=0.05,
            size=len(values)
        )

        x = np.full(len(values), i) + jitter

        ax.scatter(
            x,
            values,
            alpha=0.6,
            s=18
        )

    ax.set_ylabel("log2(CPM + 1)")
    ax.set_xlabel("Cell-of-origin subgroup")
    ax.set_title(f"{gene} expression in DLBCL")

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{gene}_expression_by_coo.png",
        dpi=300
    )

    plt.close(fig)

print("Finished.")
