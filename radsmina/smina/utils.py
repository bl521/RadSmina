#!/usr/bin/env python3

import os
import tempfile

def build_offset_index_for_mol2(mol2_path):
    """
    Return { ligand_id: offset_in_file } for each block in mol2_path.
    """
    offset_index = {}
    with open(mol2_path, 'rb') as f:
        while True:
            current_offset = f.tell()
            line_bytes = f.readline()
            if not line_bytes:
                break
            line_str = line_bytes.decode('utf-8', errors='replace')
            if line_str.startswith("##########") and "Name:" in line_str and "Long Name:" not in line_str:
                parts = line_str.split("Name:")
                if len(parts) > 1:
                    ligand_id = parts[1].strip()
                    offset_index[ligand_id] = current_offset
    return offset_index

def get_ligand_block(mol2_path, offset):
    """
    Reads lines from 'mol2_path' starting at 'offset' until next "########## Name:" or EOF.
    Returns the list of lines.
    """
    block_lines = []
    with open(mol2_path, 'rb') as f:
        f.seek(offset)
        while True:
            line_bytes = f.readline()
            if not line_bytes:
                break  # EOF
            line_str = line_bytes.decode('utf-8', errors='replace')

            # If we see the next "########## Name:" and we already read at least 1 line,
            # that means a new block is starting -> stop
            if (line_str.startswith("##########") and "Name:" in line_str 
                and "Long Name:" not in line_str
                and len(block_lines) > 0):
                break

            block_lines.append(line_str)

    return block_lines

def build_goldilocks_index(folder):
    """
    Look through all the .mol2 file in the folder.
    Create offset index dictionary for each molecule in the mol2 file.
    """
    big_index = {}
    for fname in os.listdir(folder):
        if not fname.endswith(".mol2"):
            continue
        mol2_path = os.path.join(folder, fname)
        if not os.path.isfile(mol2_path):
            continue

        offsets = build_offset_index_for_mol2(mol2_path)
        for lid, off in offsets.items():
            big_index[lid] = (mol2_path, off)

    return big_index

def retrieve_mol2_file(ligand_id, big_index):
    """
    1) Look up (mol2_path, offset) from 'big_index[ligand_id]'.
    2) Extract lines from .mol2 file at that offset (using get_ligand_block).
    3) Write them to a temporary .mol2 file.
    4) Return path to that temp file or None if fails.
    """
    if ligand_id not in big_index:
        return None
    mol2_path, offset = big_index[ligand_id]

    block = get_ligand_block(mol2_path, offset)  # from your existing function
    if not block:
        return None

    # Write to a temp file
    with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as tmp:
        tmp_mol2 = tmp.name
    with open(tmp_mol2, "w") as out:
        out.writelines(block)
    return tmp_mol2