import pandas as pd
from scipy.stats import mannwhitneyu

expr = pd.read_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t",
    index_col=0
)

meta = pd.read_csv(
    "metadata/gse147986_coo_metadata.tsv",
    sep="\t"
)

coo_map = dict(
    zip(meta["gsm"], meta["coo"])
)

rows = []

for sample in expr.columns:

    gsm = sample.split("_")[0]

    rows.append({
        "sample": sample,
        "gsm": gsm,
        "coo": coo_map[gsm],
        "IFI16": expr.loc["IFI16", sample],
        "AIM2": expr.loc["AIM2", sample],
        "MNDA": expr.loc["MNDA", sample],
        "PYHIN1": expr.loc["PYHIN1", sample]
    })

df = pd.DataFrame(rows)

df.to_csv(
    "results/normalised_expression/pyhin_log_cpm_by_sample.tsv",
    sep="\t",
    index=False
)

genes = [
    "IFI16",
    "AIM2",
    "MNDA",
    "PYHIN1"
]

stats_rows = []

for gene in genes:

    abc = df.loc[
        df["coo"] == "ABC",
        gene
    ]

    gcb = df.loc[
        df["coo"] == "GCB",
        gene
    ]

    stat, p = mannwhitneyu(
        abc,
        gcb,
        alternative="two-sided"
    )

    stats_rows.append({
        "gene": gene,
        "ABC_median": abc.median(),
        "GCB_median": gcb.median(),
        "difference_ABC_minus_GCB":
            abc.median() - gcb.median(),
        "p_value": p
    })

stats = pd.DataFrame(stats_rows)

stats.to_csv(
    "results/normalised_expression/pyhin_abc_vs_gcb_mannwhitney.tsv",
    sep="\t",
    index=False
)

print(stats)
