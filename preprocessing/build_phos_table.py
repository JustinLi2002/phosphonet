import numpy as np
import pandas as pd
from pathlib import Path


## Read Node2Vec output into dataframe
"""
Read the node2vec embedding text file
First line = "18357 128"  → network size + embedding dim
Remaining lines: 9606.ENSP... 128 floats
"""
lut = {}                                               # lookup-table: EnsemblID → 128-dim vector

with open("/home/FCAM/juli/HRP/results/ppi_node2vec.emb.txt") as f:
    _ = f.readline()                                   # skip the count header
    for line in f:                                     # iterate over every embedding row
        items = line.split(" ")                        # whitespace-split 130 items (ID + 128 floats + '\n')
        pid    = items[0]                              # e.g. "9606.ENSP00000297591"
        values = np.array(items[1:], dtype=float)      # convert 128 strings → float array
        lut[pid] = values                              # store in dict

print(list(lut.keys())[:10])
print(len(lut))

## convert lut into a dataframe
"""
Build a DataFrame whose index is the bare ENSP ID
"""

keys = [k[5:] for k in lut.keys()]  # strip the "9606." prefix → "ENSP00000297591"

nodevec = pd.DataFrame(
    {
        "Ensembl_ID": keys,
        "Embedding": lut.values(),
    }
)

nodevec


## Read the Phosphsite table

# The table is pre-joined with GPS to have Kinase scores

phos_df = pd.read_csv("phossite_gps.csv")
phos_df

## Join the PhosSite table with the node2vec table
"""
The ID used are not the same. We will rely on mapping table to match IDs in these

PhosSite table has UniProt ID
nodevec table has Emsembl_ID

- Add an Emsembl_ID column in PhosSite Table (join with mapping table)
- Add embedding values into PhosSite Table (join with nodevec) rely on emsembl_ID
"""

src = phos_df.set_index("ACC_ID")
src

dst = mapping.rename({"Entry":"ACC_ID"}, axis=1)
dst = dst.set_index("ACC_ID")
dst

np.array(mapping.loc[mapping["Entry"] == "Q69YN4", "From"]).astype(str)[0]

## Create a LUT dict. Key: PhosGrpID, Value: Nodevec embeddings

from tqdm import tqdm

nodevec_dict = {}
for grp_id in tqdm(np.unique(phos_df["SITE_GRP_ID"])):
    embeddings = []
    for uniprot_id in phos_df.loc[phos_df["SITE_GRP_ID"] == grp_id, "ACC_ID"]:
        ensp_ids = mapping.loc[mapping["Entry"] == uniprot_id, "From"]

        for ensp_id in np.array(ensp_ids).astype(str):
            try:
                embeddings.append(lut["9606." + ensp_id])
            except:
                pass
    
    if len(embeddings) > 0:
        nodevec_dict[grp_id] = np.mean(embeddings, axis=0)
      
print(f"length: {len(nodevec_dict)}")
print(nodevec_dict[3426383])

## Create new dictionary mapping PhosGrpID to full embedding
"""
The full embedding is two-part embedding:

- nodevec id based on protein id
- Kinase scores based on Kinpred 
"""

def get_embedding(phospho_group_id, *, df = phos_df):
    """ return embedding values as a numpy array, based on phospho group id
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

grp_id_lut = {}

for grp_id, nodevec_embedding in tqdm(nodevec_dict.items()):
    kinpred_scores = get_embedding(grp_id)

    if not np.any(np.isnan(kinpred_scores)):
        grp_id_lut[str(grp_id)] = np.r_[kinpred_scores, nodevec_embedding].tolist()

print(len(grp_id_lut))
print(grp_id_lut["3426383"])


# save as a json file

import json

with open("grp_id_to_embedding.json", "w") as f:
    json.dump(grp_id_lut, f)
