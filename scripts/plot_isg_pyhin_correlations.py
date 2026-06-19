import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path

scores = pd.read_csv(
    "results/normalised_expression/isg_scores.tsv",
    sep="\t"
)

pyhin = pd.read_csv(
    "results/normalised_expression/pyhin_log_cpm_by_sample.tsv",
    sep="\t"
).rename(columns={"sample": "sample_id"})

df = scores.merge(
    pyhin[["sample_id", "IFI16", "AIM2", "MNDA", "PYHIN1"]],
    on="sample_id",
    how="left"
)

out_dir = Path("figures/interferon_scores")
out_dir.mkdir(parents=True, exist_ok=True)

genes = ["IFI16", "AIM2", "MNDA", "PYHIN1"]
groups = ["ABC", "GCB", "Type_3"]

for gene in genes:
    rho, p = spearmanr(df["isg_score"], df[gene], nan_policy="omit")

    fig, ax = plt.subplots(figsize=(5, 4))

    for group in groups:
        subset = df[df["coo"] == group]
        ax.scatter(
            subset[gene],
            subset["isg_score"],
            label=group,
            alpha=0.7,
            s=28
        )

    ax.set_xlabel(f"{gene} expression, log2(CPM + 1)")
    ax.set_ylabel("ISG programme score")
    ax.set_title(f"{gene} and ISG programme activity")
    ax.text(
        0.05,
        0.95,
        f"Spearman rho = {rho:.2f}\np = {p:.2e}",
        transform=ax.transAxes,
        verticalalignment="top"
    )
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(
        out_dir / f"{gene}_isg_score_correlation.png",
        dpi=300
    )
    plt.close(fig)

print("Finished.")
