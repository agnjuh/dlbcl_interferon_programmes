import pandas as pd
from scipy.stats import fisher_exact
from pathlib import Path

de = pd.read_csv(
    "results/normalised_expression/mnda_high_vs_low_de.tsv",
    sep="\t"
)

with open("metadata/gene_sets/pyhin_interferon_target_genes.txt") as f:
    interferon_genes = {line.strip() for line in f if line.strip()}

background_genes = set(de["gene"])
interferon_genes = interferon_genes & background_genes

significant = set(
    de.loc[
        (de["fdr"] < 0.05) &
        (de["logCPM_difference"] > 0),
        "gene"
    ]
)

a = len(significant & interferon_genes)
b = len(significant - interferon_genes)
c = len(interferon_genes - significant)
d = len(background_genes - significant - interferon_genes)

table = [[a, b], [c, d]]

odds_ratio, p_value = fisher_exact(
    table,
    alternative="greater"
)

result = pd.DataFrame([{
    "gene_set": "PYHIN_interferon_target_genes",
    "background_genes": len(background_genes),
    "gene_set_size": len(interferon_genes),
    "significant_upregulated_genes": len(significant),
    "overlap": a,
    "odds_ratio": odds_ratio,
    "p_value": p_value
}])

out_dir = Path("results/enrichment")
out_dir.mkdir(parents=True, exist_ok=True)

result.to_csv(
    out_dir / "mnda_high_interferon_gene_enrichment.tsv",
    sep="\t",
    index=False
)

overlap_genes = sorted(significant & interferon_genes)

pd.DataFrame({
    "gene": overlap_genes
}).to_csv(
    out_dir / "mnda_high_interferon_gene_overlap.tsv",
    sep="\t",
    index=False
)

print(result)
print()
print("Overlap genes:")
print(", ".join(overlap_genes))
