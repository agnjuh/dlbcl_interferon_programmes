import pandas as pd

expr = pd.read_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t",
    index_col=0
)

mnda = expr.loc["MNDA"]

q25 = mnda.quantile(0.25)
q75 = mnda.quantile(0.75)

groups = []

for sample, value in mnda.items():

    if value <= q25:
        groups.append({
            "sample_id": sample,
            "mnda_expression": value,
            "group": "MNDA_low"
        })

    elif value >= q75:
        groups.append({
            "sample_id": sample,
            "mnda_expression": value,
            "group": "MNDA_high"
        })

groups = pd.DataFrame(groups)

groups.to_csv(
    "results/normalised_expression/mnda_groups.tsv",
    sep="\t",
    index=False
)

print(groups["group"].value_counts())
print()
print(groups.head())
