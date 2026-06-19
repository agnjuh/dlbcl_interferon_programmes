import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

expr = pd.read_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t",
    index_col=0
)

groups = pd.read_csv(
    "results/normalised_expression/mnda_groups.tsv",
    sep="\t"
)

high = groups.loc[
    groups["group"] == "MNDA_high",
    "sample_id"
].tolist()

low = groups.loc[
    groups["group"] == "MNDA_low",
    "sample_id"
].tolist()

results = []

for gene in expr.index:

    high_values = expr.loc[gene, high]
    low_values = expr.loc[gene, low]

    stat, p = mannwhitneyu(
        high_values,
        low_values,
        alternative="two-sided"
    )

    results.append({
        "gene": gene,
        "high_median": high_values.median(),
        "low_median": low_values.median(),
        "logCPM_difference":
            high_values.median() - low_values.median(),
        "p_value": p
    })

res = pd.DataFrame(results)

res["fdr"] = multipletests(
    res["p_value"],
    method="fdr_bh"
)[1]

res = res.sort_values("fdr")

res.to_csv(
    "results/normalised_expression/mnda_high_vs_low_de.tsv",
    sep="\t",
    index=False
)

print(res.head(25))
