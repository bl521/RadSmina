import pickle
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import time
from tqdm import tqdm
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

from rad.construction import getGraphs
from rad.traversal import traverseHNSW

from smina.dock import dock_with_smina
from smina.utils import build_offset_index_for_mol2, get_ligand_block, build_goldilocks_index, retrieve_mol2_file

# Disable rdkit error logging to keep output clean
RDLogger.DisableLog('rdApp.*')  

### Load the DUDEZ DOCK dataset
input_pickle = "goldilocks_smiles.pkl"
with open(input_pickle, 'rb') as f:
    dudez_data = pickle.load(f)
print(f"Loaded {len(dudez_data)} molecules from {input_pickle}.")


### Set parameters for fingerprints and generate them
FP_LENGTH = 1024
FP_RADIUS = 2
failed_smiles = []  # List to store failed SMILES
successful_count = 0  # Counter for successful molecules
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


### Set parameters for HNSW and construct it
EF_CONSTRUCTION = 400 ## graph_quality
M = 16 ## max_neighbours_node
hnsw_layer_graphs = getGraphs(dudez_fps, ef_construction=EF_CONSTRUCTION, M=M)
for layer_idx, layer in enumerate(hnsw_layer_graphs):
    print(f"Layer {layer_idx} has {len(layer)} nodes.")

dudez_dir = "/rds/general/user/bl521/home/smina_docking/dudez"
big_index = build_goldilocks_index(dudez_dir)
print(len(big_index))

RECEPTOR = "ROCK1"
NUM_TO_TRAVERSE = 10000 # Maximum number of molecules to score
# Load the JSON that maps receptor names -> rec.crg.pdb paths
json_file = "/rds/general/user/bl521/home/rad/examples/pdb_receptors.json"
with open(json_file, "r") as f:
    receptor_map = json.load(f)
# Load the JSON that maps target names -> xtal-lig.pdb paths
json_file = "/rds/general/user/bl521/home/rad/examples/pdb_tar_lig.json"
with open(json_file, "r") as f:
    target_map = json.load(f)
if RECEPTOR not in receptor_map:
    raise ValueError(f"[ERROR] Receptor '{RECEPTOR}' not found in {json_file}")
receptor_pdb = receptor_map[RECEPTOR]
target_pdb = target_map[RECEPTOR]

scores_by_node = {}  # node_id -> docking score
def score_fn(node_id):
    start_all = time.perf_counter()

    # 1) bounds check + get ligand ID
    stage1_start = time.perf_counter()
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
    stage1_end = time.perf_counter()
#     print(f"[DEBUG] Found offset={offset} in {mol2_path} for {ligand_id}")

    # 2) retrieve temp .mol2
    stage2_start = time.perf_counter()
    mol2_temp = retrieve_mol2_file(ligand_id, big_index)
    stage2_end = time.perf_counter()
    
#     block = get_ligand_block(mol2_path, offset)
    # Write and debug-print the block
#     mol2_temp1 = write_block_to_temp_mol2(block, ligand_id)
    if not mol2_temp:
        # either not found or offset read failed
        return np.inf

    # 3) dock
    stage3_start = time.perf_counter()
    best_score, _ = dock_with_smina(mol2_temp, receptor_pdb, target_pdb)
    stage3_end = time.perf_counter()
    # remove the temp file
    stage4_start = time.perf_counter()
    os.remove(mol2_temp)
    stage4_end = time.perf_counter()

    if best_score is None:
        print(f"[DEBUG] dock_with_smina returned None for {ligand_id}")
        return np.inf

    scores_by_node[ligand_id] = [node_id, smi, best_score]
    
    return best_score

traversed_nodes = traverseHNSW(hnsw_layer_graphs, score_fn, NUM_TO_TRAVERSE)

with open("zid_scores_rock1.json", "w") as f:
    json.dump(scores_by_node, f, indent=2)