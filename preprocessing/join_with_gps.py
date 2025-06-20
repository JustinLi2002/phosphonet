from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

GPS = Path("/home/FCAM/juli/HRP/data/Final/GPS/GPS5_2020-02-26_all_matrix.csv")

"""
To get the kinase-specific scores for every potential phosphorylation site in the proteome
"""
df = pd.read_csv(GPS)
df

"""
Different residue types  have different kinase families. Only model Y-sites 
"""
df_Y = df.loc[df["site"].str.startswith("Y")] # Filters rows where the site column starts with “Y”.
df_Y

"""
Loops through every kinase column, checks if the column is entirely “-”, collects those names, then drops them. We are doing this because columns with zero information hurt both runtime and statistics. They are useless for analysis
"""
kinases = list(df.columns)[5:] # all kinase columns
droped = []
for kinase in kinases:
    if np.all(df_Y[kinase] == "-"): # column is **all** dashes
        droped.append(kinase)
df_Y = df_Y.drop(droped, axis=1)
df_Y

"""
load the raw dataset
"""

phos_df = pd.read_csv(
    "/home/FCAM/juli/HRP/Phosphorylation_site_dataset", 
    sep="\t",
)
phos_df

"""
Keep human tyrosine sites only
"""

phos_df = phos_df.loc[phos_df["MOD_RSD"].str.startswith("Y")]
phos_df = phos_df.loc[phos_df["ORGANISM"] == "human"]
phos_df

"""
Converts Y106-p to Y106. Makes sites comparable with GPS notation.
"""
phos_df["MOD_RSD"] = phos_df["MOD_RSD"].str.removesuffix("-p")
phos_df

# Renames columns 

dst = df_Y.rename({"substrate_acc":"ACC_ID", "site": "MOD_RSD"}, axis=1)
dst = dst.set_index(["ACC_ID", "MOD_RSD"])
dst

"""
Adds the 65 kinase-score columns + 4 GPS metadata cols onto the PSP scaffold. Produces a table for analysis
"""
merged = src.join(dst)
merged

merged.to_csv("phossite_gps.csv")

def get_embedding(phospho_group_id, *, df = merged):
    """ return embedding as a numpy array, based on phospho group id
    """

    cols = ['ABL1', 'ABL2', 'ALK', 'AXL', 'BAZ1B', 'BLK', 'BMX', 'BTK', 'CSF1R',
       'CSK', 'DDR2', 'EGFR', 'EPHA2', 'EPHA3', 'EPHA4', 'EPHA8', 'EPHB1',
       'EPHB2', 'ERBB2', 'ERBB4', 'FER', 'FES', 'FGFR1', 'FGFR2', 'FGFR3',
       'FGFR4', 'FGR', 'FLT1', 'FLT3', 'FLT4', 'FRK', 'FYN', 'HCK', 'IGF1R',
       'INSR', 'ITK', 'JAK1', 'JAK2', 'JAK3', 'KDR', 'KIT', 'LCK', 'LYN',
       'MATK', 'MERTK', 'MET', 'MST1R', 'MUSK', 'NTRK1', 'NTRK2', 'PDGFRA',
       'PDGFRB', 'PTK2', 'PTK2B', 'PTK6', 'RET', 'SRC', 'SYK', 'TEC', 'TEK',
       'TNK2', 'TXK', 'TYK2', 'YES1', 'ZAP70']

    rows = df.loc[df["SITE_GRP_ID"] == phospho_group_id]
    
    if len(rows) == 0:
        raise ValueError("Invalid group id")

    features = rows[cols].to_numpy().astype(float)

    features = features.mean(axis = 0)

    return features
