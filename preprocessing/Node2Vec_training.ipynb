import os, gzip, pandas as pd, networkx as nx, numpy as np
from tqdm import tqdm
from node2vec import Node2Vec


DATA_FILE   = "/home/FCAM/juli/HRP/data/9606.protein.physical.links.full.v12.0.txt.gz"
OUT_DIR     = "/home/FCAM/juli/HRP/results"
CUTOFF      = 200

DIMS        = 128
WALK_LENGTH = 80
NUM_WALKS   = 20
SEED        = 42

GRID_PQ     = [(1, 0.5), (1, 1), (4, 1), (4, 4)]

os.makedirs(OUT_DIR,           exist_ok=True)
TMP_DIR = os.path.join(OUT_DIR, "tmp")
os.makedirs(TMP_DIR,           exist_ok=True)


chunks = []
with gzip.open(DATA_FILE, "rt") as fh:
    reader = pd.read_csv(
        fh,
        sep=r"\s+",
        comment="#",                     
        header=0,                       
        usecols=["protein1", "protein2", "combined_score"],
        dtype={
            "protein1":        "category",
            "protein2":        "category",
            "combined_score":  "int32",
        },
        chunksize=2_000_000,
        low_memory=False,
    )
    for chunk in tqdm(reader, desc="Loading"):
        chunk = chunk[chunk.combined_score >= CUTOFF]
        chunks.append(chunk)

df = pd.concat(chunks, ignore_index=True)
print(f"✔ kept {len(df):,} edges ≥{CUTOFF}")


# Check graph size

# Note: the df edge list shows each edge twice because the graph is un-ddierctional

n_nodes = len(np.unique(df["protein1"]))

print(f"Number of edges: {len(df)}")
print(f"Number of Nodes: {n_nodes}")

G = nx.from_pandas_edgelist(
    df, "protein1", "protein2", edge_attr="combined_score", create_using=nx.Graph
)
print(f"✔ Graph: {G.number_of_nodes():,} nodes • {G.number_of_edges():,} edges")

best_model, best_score, best_pq = None, -np.inf, None
edges_list = list(G.edges())
sample_pos = edges_list[:10_000] if len(edges_list) >= 10_000 else edges_list

def mean_similarity(model, pairs):
    """ 
    Inputs:
        model: Node2Vec class represent the learnin model
        pairs: [(u1, v1), (u2, v2), ... (un, vn)]
    """
    sims = [model.wv.similarity(u, v) for u, v in pairs]
    return float(np.mean(sims))

for p, q in GRID_PQ:
    print(f"→ training node2vec with p={p}, q={q}")
    n2v = Node2Vec(
        G,
        dimensions=DIMS,
        walk_length=WALK_LENGTH,
        num_walks=NUM_WALKS,
        p=p, q=q,
        weight_key="combined_score",
        workers=os.cpu_count(),
        temp_folder=TMP_DIR,
        seed=SEED,
        quiet=True,
    )
    model = n2v.fit(window=10, min_count=1, batch_words=2048, epochs=1)

    score = mean_similarity(model, sample_pos)
    print(f"   mean edge-sim {score:.4f}")

    if score > best_score:
        best_model, best_score, best_pq = model, score, (p, q)

print(f"✔ best (p, q) = {best_pq}  |  mean edge-sim = {best_score:.4f}")


txt_path = os.path.join(OUT_DIR, "ppi_node2vec.emb.txt")
bin_path = os.path.join(OUT_DIR, "ppi_node2vec.kv")

best_model.wv.save_word2vec_format(txt_path)   
best_model.wv.save(bin_path)                  

print("✔ saved text   ➜", txt_path)
print("✔ saved binary ➜", bin_path)

node0 = best_model.wv.index_to_key[0]
print("vector(", node0, ")[:5] =", best_model.wv[node0][:5])
