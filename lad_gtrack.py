import bioframe as bf
import numpy as np
import sys
import pandas as pd

gtrack = bf.read_table(sys.argv[1])
bed = bf.read_table(sys.argv[2])
bed = bed.loc[ : , :2]
bed.columns = ["chrom", "start", "end"]

gtrack.columns = ["chrom", "start", "end", "id", "radius", "periphery", "edges"]
org_chrom_names = gtrack["chrom"].to_numpy()


overlap = bf.overlap(gtrack, bed, suffixes=('_1','_2'), how='left', return_overlap = True)
overlap['overlap_length'] = overlap['overlap_end'] - overlap['overlap_start']
result = (
    overlap.groupby("id_1", as_index=False)
    .agg({
        "chrom_1": "first",
        "start_1": "first",
        "end_1": "first",
        "radius_1": "first",
        "periphery_1": "first",
        "edges_1": "first",
        "overlap_length": "sum",
    })
)


result['overlap_percentage'] = (result['overlap_length'] / (result['end_1'] - result['start_1']))
result["periphery_1"] = np.where(result['overlap_percentage'].values >= 0.5, 1, ".")
result = result.iloc[:, [1, 2, 3, 0, 4, 5, 6]]

result.columns.values[0:9] = ["###seqid", "start", "end", "id", "radius", "periphery", "edges"]

chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
result["###seqid"] = pd.Categorical(result["###seqid"], categories=chrom_order, ordered=True)
result = result.sort_values(by=["###seqid", "start", "end"])

result.to_csv(sys.argv[3], index=False, sep="\t")
