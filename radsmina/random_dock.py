#!/usr/bin/env python3

import os
import json
import random
import subprocess
import tempfile
from tqdm import tqdm  # <-- import tqdm
from smina.dock import dock_with_smina

def get_ligand_blocks(mol2_path):
    """
    Yields blocks (list of lines) for each ligand in mol2_path.
    A block starts at lines containing '########## Name:' 
    and continues until the next such line or EOF.
    """
    with open(mol2_path, 'r') as f:
        current_block = []
        for line in f:
            if (line.startswith("##########")
                and "Name:" in line
                and "Long Name:" not in line):
                if current_block:
                    yield current_block
                current_block = [line]
            else:
                current_block.append(line)
        if current_block:
            yield current_block

def extract_id_and_smiles(block):
    """
    From a block of lines, extract:
      - ligand_id from 'Name:  ZINC0000..'
      - SMILES from 'SMILES:  C1=CC..'
    Return (ligand_id, smiles).
    """
    ligand_id = None
    smiles = None
    for line in block:
        if "Name:" in line and "Long Name:" not in line:
            parts = line.split("Name:")
            if len(parts) > 1:
                ligand_id = parts[1].strip()
        if "SMILES:" in line:
            parts = line.split("SMILES:")
            if len(parts) > 1:
                smiles = parts[1].strip()
    return (ligand_id, smiles)

def build_ligand_dict(folder):
    """
    Parses all 'super_*.mol2' files in 'folder', 
    storing { ligand_id: (smiles, block_of_lines) }.
    """
    big_dict = {}
    for fname in os.listdir(folder):
        if not fname.endswith(".mol2"):
            continue
        mol2_path = os.path.join(folder, fname)
        if not os.path.isfile(mol2_path):
            continue
        
        for block in get_ligand_blocks(mol2_path):
            lid, smi = extract_id_and_smiles(block)
            if lid and smi:
                # store the entire block so we can dock
                big_dict[lid] = (smi, block)
    
    return big_dict

def write_block_to_temp_mol2(block):
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".mol2", delete=False)
    tmp_path = tmp.name
    tmp.close()
    with open(tmp_path, "w") as f:
        f.writelines(block)
    return tmp_path

def main():
    # 1) Build big dictionary from 'super_goldilocks/super_*.mol2'
    super_folder = "/rds/general/user/bl521/home/smina_docking/super_goldilocks"
    big_dict = build_ligand_dict(super_folder)
    print(f"[INFO] Found {len(big_dict)} ligands across all super_*.mol2 files.")
    
    # 2) Randomly sample 10,000 
    sample_size = 10000
    all_ids = list(big_dict.keys())
    chosen_ids = random.sample(all_ids, min(sample_size, len(all_ids)))
    
    # 3) Choose the receptor and ligand
    receptor_pdb = "/rds/general/user/bl521/home/smina_docking/dudez/ROCK1/rec.crg.pdb"
    target_pdb   = "/rds/general/user/bl521/home/smina_docking/dudez/ROCK1/xtal-lig.pdb"
    
    # 4) Dock them, store results in dictionary
    results = {}  # { ligand_id: [SMILES, best_score] }
    for lid in tqdm(chosen_ids, desc="Docking ligands"):
        smi, block = big_dict[lid]
        tmp_mol2 = write_block_to_temp_mol2(block)
        best_score = dock_with_smina(tmp_mol2, receptor_pdb, target_pdb)
        os.remove(tmp_mol2)
        if best_score is not None:
            results[lid] = [smi, best_score]
        else:
            results[lid] = [smi, None] 
    
    # 5) Save results to JSON
    out_file = "random_10k_dock_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"[INFO] Docked {len(chosen_ids)} molecules, wrote {out_file}")

if __name__ == "__main__":
    main()
