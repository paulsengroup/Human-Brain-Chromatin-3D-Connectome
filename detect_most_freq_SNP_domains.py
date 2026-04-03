import bioframe as bf
import pandas as pd
import argparse
from scipy.spatial.distance import pdist, squareform
import os



parser = argparse.ArgumentParser(description="Process SNP")
parser.add_argument("--model", required=True, help="Path to SNP table")
args = parser.parse_args()
sample_name = os.path.basename(args.model).replace(".cmm.tsv", "")

model = bf.read_table(args.model)

snp = pd.read_csv("snp_coordinates.csv", delimiter=",")
snp["end"] = snp["pos"] + 1

snp = snp.iloc[:, [1, 2, 7, 4, 3]]
snp.columns = snp.columns.astype(str)

snp.columns.values[0] = "chrom"
snp.columns.values[1] = "start"
snp.columns.values[2] = "end"
snp.columns.values[3] = "snp"

snp['chrom'] = snp['chrom'].str[3:]

model = model.iloc[ : , :6]
model.columns = model.columns.astype(str)
model.columns.values[0] = "chrom"
model.columns.values[1] = "start"
model.columns.values[2] = "end"
model.columns.values[3] = "x"
model.columns.values[4] = "y"
model.columns.values[5] = "z"

overlap = bf.overlap(model, snp, suffixes=('_1','_2'))

overlap = (
    overlap
    .groupby(['chrom_1', 'start_1', 'end_1', 'x_1', 'y_1', 'z_1'], dropna=False)
    .agg({'snp_2': lambda x: ','.join(x.dropna().astype(str)), 
         'trait_2': lambda x: ','.join(x.dropna().astype(str))})
    .reset_index()
)


overlap['region'] = overlap['chrom_1'].astype(str) + ':' + overlap['start_1'].astype(str) + '-' + overlap['end_1'].astype(str)
overlap = overlap[['trait_2', 'region', 'x_1', 'y_1', 'z_1']]


overlap["count_immune"] = (
    overlap["trait_2"].str.count(r"\bMS\b") +
    overlap["trait_2"].str.count(r"\bIBD\b")
)

overlap["count_neurodeg"] = (
    overlap["trait_2"].str.count(r"\bAD\b") +
    overlap["trait_2"].str.count(r"\bPD\b")
)

top_number = 30

overlap = overlap.sort_values(by="count_immune", ascending=False)
threshold_immune = overlap["count_immune"].iloc[top_number]

overlap = overlap.sort_values(by="count_neurodeg", ascending=False)
threshold_neurodeg = overlap["count_neurodeg"].iloc[top_number]

filtered = overlap[
    (overlap["count_immune"] >= threshold_immune) |
    (overlap["count_neurodeg"] >= threshold_neurodeg)
]

all_distances = squareform(
    pdist(filtered[["x_1", "y_1", "z_1"]], metric="euclidean")
)

all_distances = pd.DataFrame(
    all_distances,
    index=filtered["region"],
    columns=filtered["region"]
)
all_distances.to_csv(f"{sample_name}_top30.tsv", sep="\t", index=False)
filtered.to_csv(f"{sample_name}_top30_beads.tsv", sep="\t", index=False)
