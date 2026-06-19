import gzip
import pandas as pd

def clean(x):
    return x.replace('"', '').strip()

with gzip.open(
    "data/metadata/GSE147986_series_matrix.txt.gz",
    "rt"
) as f:
    lines = f.readlines()

titles = None
gsms = None
coo = None

for line in lines:

    if line.startswith("!Sample_title"):
        titles = [clean(x) for x in line.strip().split("\t")[1:]]

    elif line.startswith("!Sample_geo_accession"):
        gsms = [clean(x) for x in line.strip().split("\t")[1:]]

    elif (
        line.startswith("!Sample_characteristics_ch1")
        and "coo subgroups:" in line
    ):
        coo = [
            clean(x).replace("coo subgroups: ", "")
            for x in line.strip().split("\t")[1:]
        ]

meta = pd.DataFrame({
    "gsm": gsms,
    "sample": titles,
    "coo": coo
})

meta.to_csv(
    "metadata/gse147986_coo_metadata.tsv",
    sep="\t",
    index=False
)

print(meta["coo"].value_counts())
