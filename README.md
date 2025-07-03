# RAD-SMINA 📈🧬  
*Retrieval-Augmented Docking with SMINA on the DUDE-Z “Goldilocks” library*

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)
[![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 1 · Project motivation
The original **RAD** workflow couples an HNSW graph with DOCK 3.7 to avoid brute-force docking of giga-scale libraries.  
This repository replaces DOCK 3.7 with **SMINA**, re-applying **RAD** on a 132 k-molecule subset of the DUDE-Z “Goldilocks” set.

## 2 · Repository layout
```bash
RadSmina/                     # project root ── an installable Python package
│
├─ plots/                     # all figure-generation code & results
│   ├─ enrichment_results/       # examples of the enrichment plots generated
│   ├─ performance_results/      # examples of the performance graphs generated
│   ├─ correlation_plots/        # examples of the correlation graphs generated
│   ├─ enrichment_plot.py       # → draws early-recall curves
│   ├─ performance.py           # → box-plots of docking score distributions
│   └─ correlation.py           # → pairwise Tanimoto similarities vs. score difference plots
│
├─ rad/                       # upstream RAD (Hall & Keiser) for reference
│
├─ radsmina/                  # **this project’s actual library code**
│   │
│   ├─ data/                  # example data to reproduce the results
│   │   ├─ receptor_files/       # PDBQT receptors used by SMINA
│   │   ├─ reclig_files/         # reference ligands (define autobox)
│   │   ├─ super_goldilocks/     # trimmed 3-D mol2s (132 k)
│   │   │                        # - *Download required* see relevant scripts for instructions
│   │   └─ goldilocks_smiles.pkl # text SMILES used for HNSW construction
│   │
│   ├─ smina/                 # **thin Python wrapper around SMINA**
│   │   ├─ temp_output/          # docking poses & log files (auto-cleaned)
│   │   ├─ __init__.py
│   │   ├─ dock.py              # run dock_with_smins(), parse scores
│   │   └─ utils.py             # helper functions for smina docking
│   │
│   ├─ scores/                # example docked socres ready for use
│   │
│   ├─ best_pose_output/      # optional: where `dock.py --keep_poses` saves .mol2
│   │
│   ├─ dockingjob.sh          # hpc example script (128 cores, 28 h wall)
│   ├─ DUDEZ_smina_with_scores.ipynb*  # **FAST replay** – reuse pre-computed SMINA scores → traverse graph → emit JSONs for plotting
│   ├─ DUDEZ_smina.ipynb*                # **FULL pipeline** – build HNSW → RAD traversal → run SMINA docking on-the-fly → write score JSONs
│   └─ DUDEZ_smina.py*                # same as DUDEZ_smina.ipynb, refactored as a Python script for HPC batch jobs
│
├─ utils/                     # project-agnostic helpers
│   └─ paths.py               # centralises folder & file paths (edit here once)
│
├─ .gitignore
├─ environment.yml
├─ README.md
└─ setup.py
```
### **Main RAD+SMINA pipeline**
| Path                                 | Purpose                                                                                                                      |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------- |
| **DUDEZ\_smina.ipynb**               | **End-to-end demo**. Builds the HNSW, performs RAD traversal, invokes SMINA live, and saves `<setting>_scores_<target>.json`.|
| **DUDEZ\_smina\_with\_scores.ipynb** | **Fast replay**. Loads a pre-computed SMINA score json, reruns the traversal without docking, and writes the same JSON output for the plotting scripts.    |
| **DUDEZ\_smina.py**                  | Command-line version of `DUDEZ_smina.ipynb`.  Designed for submission to HPC via `dockingjob.sh`—no Jupyter kernel required. |
| **dockingjob.sh**                    | Example HX1 PBS script: requests 128 cores, 28 h wall-time, and executes `python3 DUDEZ_smina.py`.                            |

> **Tip:**
> 1. Use `DUDEZ_smina_with_scores.ipynb` if you only want to regenerate the plots without waiting for docking; run `DUDEZ_smina.ipynb` for the full RAD-SMINA workflow.
> 2. Use `DUDEZ_smina.py` + `dockingjob.sh` for large production runs on HX1; use the notebooks for interactive experimentation and debugging.

## 3 · Quick start
```bash
# 1  Clone
git clone https://github.com/bl521/RadSmina.git
cd RadSmina

# 2  Create conda env
conda env create -f environment.yml
conda activate radsmina

# 3 Install the code in editable mode
pip install -e .
```

## 4 · What Gets Saved After a Traversal?
The `traverseHNSW()` helper walks the graph, calls the `score_fn`, and returns a dictionary `scores_by_node`.\
At the end of each run we persist that dictionary as a single JSON file:
```python
traversed_nodes = traverseHNSW(hnsw_layer_graphs,
                               score_fn,
                               NUM_TO_TRAVERSE)

with open("zid_scores_rock1.json", "w") as f:
    json.dump(scores_by_node, f, indent=2)
```
### File location & naming
Default path → `RadSmina/radsmina/` (where your jupyter notebook is located).

Recommended pattern → `<setting>_scores_<target>_.json`, e.g.`m16ef400_scores_rock1.json`.

### JSON schema
| Key (str)   | Value (list)                       | Example                                           |
| ----------- | ---------------------------------- | ------------------------------------------------- |
| `ligand_id` | `[node_id, SMILES, docking_score]` | `"ZINC00012345": [8261, "CCN(C)C(=O)...", -11.2]` |

`node_id` – integer index of the ligand in the HNSW graph

`SMILES` – canonical SMILES used for fingerprinting

`docking_score` – best score returned by SMINA (`np.inf` if docking failed)

### Re-using the files
The plotting utilities (`enrichment_plot.py`, `performance.py`,
`correlation.py`) assume exactly this layout; no edits are required as
long as the filenames and schema are unchanged.

If you wish to keep multiple traversals for the same target, simply give
each JSON a unique filename (e.g. `m32ef400_scores_rock1_repeat2.json`);
the scripts accept a list of paths.


## 5 · Reproducing the paper plots
You can regenerate every figure that appears in the report from the command line — no notebook editing required.
All plot drivers live under `plots/` and point to the JSON score files produced by radsmina (or the bundled examples in `radsmina/scores/`).\
If you want to redirect the scripts to your own results, change all the json files in the scripts and make them to your own result paths.
| Script (`plots/…`)   | Figure(s) it produces             | What it shows                                                                                                      |        
| -------------------- | --------------------------------- | ------------------------------------------------------------------------------------------------------------------ | 
| `enrichment_plot.py` | `coverage_<setting>.png`          | Early-recall/enrichment curves (RAD vs. random traversal)                                                          |        |                                |
| `performance.py`     | `performance_<setting>_plot.png` | Box-plots of the full docking-score distributions, highlighting the 10 lowest scores with a frequency colour-scale  |        |                                |
| `correlation.py`     | `correlationplot<num>/hexbinplot<num>.png`          | Pair-wise Tanimoto similarity vs. Score scatter & Pearson-r annotation                           |

## 6 · Known issues / TODO
| ID     | Item                                                                                                  | Current status                                                                                                                                   | Planned action                                                                                                                                                                                                                                                                                                                                               |
| ------ | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **#1** | **`super_goldilocks/` omitted** – the 132 k trimmed `.mol2` files are *not* in the repository (≈700 MB) |SOLVED| N/A |


