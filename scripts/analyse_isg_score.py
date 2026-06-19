import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from pathlib import Path

scores = pd.read_csv(
    "results/normalised_expression/isg_scores.tsv",
    sep="\t"
)

pyhin = pd.read_csv(
    "results/normalised_expression/pyhin_log_cpm_by_sample.tsv",
    sep="\t"
)

pyhin = pyhin.rename(columns={"sample": "sample_id"})

df = scores.merge(
    pyhin[["sample_id", "IFI16", "AIM2", "MNDA", "PYHIN1"]],
    on="sample_id",
    how="left"
)

abc = df.loc[df["coo"] == "ABC", "isg_score"]
gcb = df.loc[df["coo"] == "GCB", "isg_score"]

stat, p = mannwhitneyu(
    abc,
    gcb,
    alternative="two-sided"
)

group_stats = (
    df.groupby("coo")["isg_score"]
    .agg(["count", "median", "mean", "std"])
    .reset_index()
)

comparison = pd.DataFrame([{
    "comparison": "ABC_vs_GCB",
    "ABC_median": abc.median(),
    "GCB_median": gcb.median(),
    "difference_ABC_minus_GCB": abc.median() - gcb.median(),
    "p_value": p
}])

correlations = []

for gene in ["IFI16", "AIM2", "MNDA", "PYHIN1"]:
    rho, p_corr = spearmanr(
        df["isg_score"],
        df[gene],
        nan_policy="omit"
    )

    correlations.append({
        "gene": gene,
        "spearman_rho_with_isg_score": rho,
        "p_value": p_corr
    })

correlations = pd.DataFrame(correlations)

out_dir = Path("results/normalised_expression")
out_dir.mkdir(parents=True, exist_ok=True)

group_stats.to_csv(
    out_dir / "isg_score_group_summary.tsv",
    sep="\t",
    index=False
)

comparison.to_csv(
    out_dir / "isg_score_abc_vs_gcb_mannwhitney.tsv",
    sep="\t",
    index=False
)

correlations.to_csv(
    out_dir / "pyhin_isg_score_correlations.tsv",
    sep="\t",
    index=False
)

print("Group summary:")
print(group_stats)
print()
print("ABC vs GCB:")
print(comparison)
print()
print("PYHIN correlations:")
print(correlations)
