import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

de = pd.read_csv(
    "results/normalised_expression/mnda_high_vs_low_de.tsv",
    sep="\t"
)

de["neglog10_fdr"] = -np.log10(
    de["fdr"].replace(0, 1e-300)
)

highlight = {
    "MNDA",
    "IFI16",
    "AIM2",
    "PYHIN1",
    "STAT1",
    "STAT2",
    "IRF7",
    "IRF9",
    "OAS1",
    "OAS2",
    "OAS3",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "MX1",
    "MX2",
    "ISG15",
    "RSAD2",
    "DDX58",
    "CXCL11"
}

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(
    de["logCPM_difference"],
    de["neglog10_fdr"],
    alpha=0.5,
    s=10
)

hits = de[de["gene"].isin(highlight)]

ax.scatter(
    hits["logCPM_difference"],
    hits["neglog10_fdr"],
    s=30
)

for _, row in hits.iterrows():
    ax.text(
        row["logCPM_difference"],
        row["neglog10_fdr"],
        row["gene"],
        fontsize=8
    )

ax.axhline(
    -np.log10(0.05),
    linestyle="--"
)

ax.set_xlabel(
    "Median log2(CPM + 1) difference\n(MNDA_high − MNDA_low)"
)

ax.set_ylabel(
    "-log10(FDR)"
)

ax.set_title(
    "MNDA-high versus MNDA-low DLBCL"
)

fig.tight_layout()

Path(
    "figures/differential_expression"
).mkdir(
    parents=True,
    exist_ok=True
)

fig.savefig(
    "figures/differential_expression/mnda_high_vs_low_volcano.png",
    dpi=300
)

plt.close()

print("Finished.")
