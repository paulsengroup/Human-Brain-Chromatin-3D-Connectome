import bioframe as bf
import pandas as pd
import sys
import os

bed_file = sys.argv[1]
base = os.path.basename(bed_file).replace("_subcom.bed", "")

model_path = f"sim_all/{base}_1.cmm.tsv"
output_path = f"sim_all/{base}_chrom3d_subcompartments.tsv"

sub_compartments = bf.read_table(bed_file)
model = bf.read_table(model_path)


sub_compartments = sub_compartments.iloc[:, [0, 1, 2, 3]]
sub_compartments.columns = sub_compartments.columns.astype(str)
sub_compartments.columns.values[0] = "chrom"
sub_compartments.columns.values[1] = "start"
sub_compartments.columns.values[2] = "end"
sub_compartments.columns.values[3] = "sub_com"
sub_compartments['sub_com'] = sub_compartments['sub_com'].str[:5]

sub_compartments.replace("A.1.1","A3", inplace=True)  
sub_compartments.replace("A.1.2","A2", inplace=True)
sub_compartments.replace("A.2.1","A1", inplace=True)
sub_compartments.replace("A.2.2","A0", inplace=True)
sub_compartments.replace("B.1.1","B0", inplace=True)
sub_compartments.replace("B.1.2","B1", inplace=True)
sub_compartments.replace("B.2.1","B2", inplace=True)
sub_compartments.replace("B.2.2","B3", inplace=True)
sub_compartments['chrom'] = sub_compartments['chrom'].str[3:]

sub_compartments['length'] = sub_compartments['end'] - sub_compartments['start']

model = model.iloc[ : , :3]
model.columns = model.columns.astype(str)
model.columns.values[0] = "chrom"
model.columns.values[1] = "start"
model.columns.values[2] = "end"

df_subset = bf.overlap(model, sub_compartments, suffixes=('_1','_2'))

df_res = df_subset.groupby(['chrom_1', 'start_1', 'end_1', "sub_com_2"], as_index=False, dropna=False).agg({'length_2':'sum'})
df_res = df_res.loc[df_res.groupby(['chrom_1', 'start_1', 'end_1'], dropna=False)['length_2'].idxmax()]

df_res['ID'] = df_res['chrom_1'] + ":" + df_res['start_1'].apply(str) + "-" + df_res['end_1'].apply(str)

df_res = df_res[["ID", "sub_com_2"]]
df_res.to_csv(output_path, index=False, sep="\t", header=False)
