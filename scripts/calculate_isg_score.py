import pandas as pd
from pathlib import Path

expr = pd.read_csv(
    "results/normalised_expression/gse147986_log_cpm.tsv",
    sep="\t",
    index_col=0
)

meta = pd.read_csv(
    "metadata/gse147986_coo_metadata.tsv",
    sep="\t"
)

with open("metadata/gene_sets/interferon_stimulated_genes.txt") as f:
    genes = [line.strip() for line in f if line.strip()]

available_genes = [g for g in genes if g in expr.index]
missing_genes = sorted(set(genes) - set(available_genes))

print(f"Using {len(available_genes)} ISG genes")

if missing_genes:
    print(f"Missing genes: {', '.join(missing_genes)}")

scores = expr.loc[available_genes].mean(axis=0)

result = pd.DataFrame({
    "sample_id": scores.index,
    "gsm": [s.split("_")[0] for s in scores.index],
    "isg_score": scores.values
})

result = result.merge(
    meta[["gsm", "coo"]],
    on="gsm",
    how="left"
)

Path("results/normalised_expression").mkdir(
    parents=True,
    exist_ok=True
)

result.to_csv(
    "results/normalised_expression/isg_scores.tsv",
    sep="\t",
    index=False
)

print()
print(result.groupby("coo")["isg_score"].median().sort_values(ascending=False))
