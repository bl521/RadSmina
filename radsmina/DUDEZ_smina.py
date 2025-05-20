#!/usr/bin/env python
# coding: utf-8
"""Dock Goldilocks SMILES against ROCK1 with HNSW + smina."""

# ── standard library ────────────────────────────────────────────────────────
import json
import os
import pickle
import random
import sys
import time
from inspect import getsource  # if you actually need it

# ── third-party ─────────────────────────────────────────────────────────────
import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm
from smina.dock import dock_with_smina
from smina.utils import (
    build_library_index,
    get_ligand_block,
    retrieve_mol2_file,
)

# ── local / first-party ─────────────────────────────────────────────────────
from rad.construction import getGraphs
from rad.traversal import traverseHNSW
from utils.paths import data_path, radsmina_path

# ---------------------------------------------------------------------------

RDLogger.DisableLog("rdApp.*")  # keep rdkit quiet



# ### Load the SMILES dataset

# In[2]:


input_pickle = "goldilocks_smiles.pkl"
with open(input_pickle, 'rb') as f:
    dudez_data = pickle.load(f)


# In[3]:


print(f"Loaded {len(dudez_data)} molecules from {input_pickle}.")


# ### Set parameters for fingerprints and generate them

# In[4]:


FP_LENGTH = 1024
FP_RADIUS = 2


# In[5]:


failed_smiles = []  # List to store failed SMILES
successful_count = 0  # Counter for successful molecules


# In[6]:


dudez_fps = [] ## To store molecular fingerprints in a compact binary form (packed bits).
node_id = 0 ## A unique identifier for each molecule (node).

# Store a mapping from node_id -> zID
id_to_zid = []

for zid in tqdm(dudez_data, total=len(dudez_data), desc="Generating Fingeprints"):
    smi = dudez_data[zid] ## SMILES string for the molecule.

    # Some smiles will fail molecule generation. We just skip them
    mol = Chem.MolFromSmiles(smi) ## The RDKit library converts the SMILES string into a molecular object
    if mol is None:
        ## Log failed SMILES
        failed_smiles.append(smi)
        continue
    
    # If successful, process as usual
    successful_count += 1
    # Convert rdkit bit vect fingerprint to numpy array
    arr = np.zeros((1,), dtype=np.uint8)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=FP_RADIUS, nBits=FP_LENGTH) ## Generate Molecular Fingerprints
    DataStructs.ConvertToNumpyArray(fp, arr) ## Convert Morgan Fingerprints to NumPy Array

    # IMPORTANT: Make sure to pack bit fingerprints - it vastly speeds up HNSW construction
    dudez_fps.append(np.packbits(arr)) 
    ## Groups every 8 bits into a single byte.
    ## Reduces the length of the array by a factor of 8 (e.g., from 1024 bits to 128 bytes).
    
    id_to_zid.append(zid) # Record that node_id corresponds to this zID
    
    node_id += 1 ## Every molecule gets a unique identifier

dudez_fps = np.array(dudez_fps) ## Convert Fingerprints List to NumPy Array
# Output results
print(f"Total SMILES processed: {len(dudez_data)}")
print(f"Successful molecule generation: {successful_count}")
print(f"Failed SMILES count: {len(failed_smiles)}")
print("Failed SMILES examples:", failed_smiles[:10])  # Print a few failed SMILES


# ### Set parameters for HNSW and construct it

# In[7]:


EF_CONSTRUCTION = 400 ## graph_quality
M = 16 ## max_neighbours_node


# In[8]:


hnsw_layer_graphs, hnsw_index = getGraphs(dudez_fps, ef_construction=EF_CONSTRUCTION, M=M)
element_levels = hnsw_index["element_levels"]
max_level = hnsw_index["max_level"]

for lvl in range(max_level + 1):
    num_nodes_in_layer = np.sum(element_levels >= lvl)
    print(f"Number of nodes in layer {lvl}: {num_nodes_in_layer}")


# ### Build file index for quicker traversal


# In[10]:


folder_dir = data_path("super_goldilocks")
big_index = build_library_index(folder_dir)
print(len(big_index))


# ### Define the receptor and receptor ligand

# In[21]:


RECEPTOR = "ROCK1"
NUM_TO_TRAVERSE = 10000 # Maximum number of molecules to score


# In[33]:


receptor_pdb = os.path.join(data_path("receptor_files"),f"{RECEPTOR}.pdb")
reclig_pdb = os.path.join(data_path("reclig_files"),f"{RECEPTOR}-lig.pdb")


# In[34]:


scores_by_node = {}  # node_id -> docking score
def score_fn(node_id):

    # 1) bounds check + get ligand ID
    if node_id < 0 or node_id >= len(id_to_zid):
        return np.inf
    ligand_id = id_to_zid[node_id]
#     print(f"[DEBUG] node_id={node_id} => ligand_id={ligand_id}")
    
#     # Check if in big_index
    if ligand_id not in big_index:
#         print(f"[DEBUG] {ligand_id} not in big_index.")
        return np.inf
    smi = dudez_data[ligand_id]
    mol2_path, offset = big_index[ligand_id]
#     print(f"[DEBUG] Found offset={offset} in {mol2_path} for {ligand_id}")

    # 2) retrieve temp .mol2
    mol2_temp = retrieve_mol2_file(ligand_id, big_index)
    
    if not mol2_temp:
        # either not found or offset read failed
        return np.inf

    # 3) dock
    best_score,mol2_text = dock_with_smina(mol2_temp, receptor_pdb, reclig_pdb)
    best_pose_output = os.path.join(radsmina_path("best_pose_output"), f"{ligand_id}_bestpose.mol2")
    with open(best_pose_output, "w") as f:
        f.write(mol2_text)
    # remove the temp file
    os.remove(mol2_temp)

    if best_score is None:
        print(f"[DEBUG] dock_with_smina_mol2 returned None for {ligand_id}")
        return np.inf
    
    scores_by_node[ligand_id] = [node_id, smi, best_score]

    
    return best_score


# In[ ]:


traversed_nodes = traverseHNSW(hnsw_layer_graphs, score_fn, NUM_TO_TRAVERSE)


# ### Store the Smina docking socres

# In[ ]:


with open("m16ef400_scores_rock1.json", "w") as f:
    json.dump(scores_by_node, f, indent=2)

