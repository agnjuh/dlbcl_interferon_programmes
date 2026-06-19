import pandas as pd
from pathlib import Path

de_file = "results/normalised_expression/mnda_high_vs_low_de.tsv"
gene_set_file = "metadata/gene_sets/pyhin_interferon_target_genes.txt"
out_file = "results/normalised_expression/mnda_high_vs_low_pyhin_interferon_hits.tsv"

de = pd.read_csv(de_file, sep="\t")

with open(gene_set_file) as f:
    target_genes = [line.strip() for line in f if line.strip()]

hits = de[de["gene"].isin(target_genes)].copy()
hits = hits.sort_values("fdr")

Path("results/normalised_expression").mkdir(parents=True, exist_ok=True)

hits.to_csv(out_file, sep="\t", index=False)

print(hits)
