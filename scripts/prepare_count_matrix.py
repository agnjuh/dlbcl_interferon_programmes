from pathlib import Path
import pandas as pd

raw_dir = Path("data/raw")
out_dir = Path("results/expression")
out_dir.mkdir(parents=True, exist_ok=True)

count_files = sorted(raw_dir.glob("GSM*.count.txt.gz"))

if not count_files:
    raise FileNotFoundError("No count files found in data/raw")

tables = []

for file in count_files:
    sample_id = file.name.replace(".count.txt.gz", "")
    df = pd.read_csv(file, sep="\t", header=None, names=["gene", sample_id])
    df = df.set_index("gene")
    tables.append(df)

count_matrix = pd.concat(tables, axis=1)
count_matrix = count_matrix.fillna(0).astype(int)

count_matrix.to_csv(out_dir / "gse147986_count_matrix.tsv", sep="\t")

sample_table = pd.DataFrame({
    "sample_id": count_matrix.columns,
    "gsm": [x.split("_")[0] for x in count_matrix.columns],
    "source_sample": ["_".join(x.split("_")[1:]) for x in count_matrix.columns],
})

sample_table.to_csv("metadata/gse147986_samples.tsv", sep="\t", index=False)

print(f"Count matrix written: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")
print("Sample table written: metadata/gse147986_samples.tsv")
