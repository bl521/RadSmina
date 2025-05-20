# RAD-SMINA 📈🧬  
*Retrieval-Augmented Docking with SMINA on the DUDE-Z “Goldilocks” library*

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
│   └─ correlation.py           # → pocket / ligand property vs recall
│
├─ rad/                       # upstream RAD (Hall & Keiser) for reference
│
├─ radsmina/                  # **this project’s actual library code**
│   │
│   ├─ data/                  # example data to reproduce the results
│   │   ├─ receptor_files/       # PDBQT receptors used by SMINA
│   │   ├─ reclig_files/         # reference ligands (define autobox)
│   │   ├─ super_goldilocks/     # trimmed 3-D mol2s (132 k)
│   │   └─ goldilocks_smiles.pkl # text SMILES used for HNSW construction
│   │
│   ├─ smina/                 # **thin Python wrapper around SMINA**
│   │   ├─ temp_output/          # docking poses & log files (auto-cleaned)
│   │   ├─ __init__.py
│   │   ├─ dock.py              # run_smina(), parse scores
│   │   └─ utils.py             # helper functions for smina docking
│   │
│   ├─ scores/                # example docked socres ready for use
│   │
│   ├─ best_pose_output/      # optional: where `dock.py --keep_poses` saves .mol2
│   │
│   ├─ dockingjob.sh          # hpc example script (128 cores, 28 h wall)
│   ├─ DUDEZ_smina_with_scores.ipynb*  # **FAST replay** – reuse pre-computed SMINA scores → traverse graph → emit JSONs for plotting
│   └─ DUDEZ_smina.ipynb*                # **FULL pipeline** – build HNSW → RAD traversal → run SMINA docking on-the-fly → write score JSONs
│
├─ rds/                       # lightweight fork of Hall & Keiser’s RAD utilities
│
├─ utils/                     # project-agnostic helpers
│   ├─ paths.py               # centralises folder & file paths (edit here once)
│   └─ __pycache__/
│
├─ .gitignore
├─ environment.yml
├─ README.md
└─ setup.py
```
### **Main RAD+SMINA pipeline**
| Path                                 | Purpose                                                                                                                                                  |
| ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **DUDEZ\_smina.ipynb**               | **End-to-end demo**. Builds the HNSW, performs RAD traversal, calls SMINA in real time, and saves the resulting `<setting>_scores_<target>.json` files.                 |
| **DUDEZ\_smina\_with\_scores.ipynb** | **Fast replay**. Loads a pre-computed SMINA score table, reruns the traversal without docking, and writes the same JSON output for the plotting scripts. |


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

## 4 · Configuration
| Parameter         | Purpose                             | Typical value |
| ----------------- | ----------------------------------- | ------------- |
| `M`               | max # neighbours per HNSW node      | **16**        |
| `ef_construction` | breadth when building HNSW          | **400**       |
| `FP_LENGTH`       | *ECFP bit-vector length*            | **1024**      |
| `FP_RADIUS`       | *ECFP bond radius*                  | **2**         |

## 5 · What Gets Saved After a Traversal?
The traverseHNSW() helper walks the graph, calls the score_fn, and returns a dictionary scores_by_node.
At the end of each run we persist that dictionary as a single JSON file:

## 6 · Reproducing the paper plots
You can regenerate every figure that appears in the report from the command line — no notebook editing required.
All plot drivers live under plots/ and point to the JSON score files produced by radsmina (or the bundled examples in scores/).
| Script (`plots/…`)   | Figure(s) it produces             | What it shows                                                                                                      |        
| -------------------- | --------------------------------- | ------------------------------------------------------------------------------------------------------------------ | 
| `enrichment_plot.py` | `coverage_<setting>.png`          | Early-recall/enrichment curves (RAD vs. random traversal)                                                          |        |                                |
| `performance.py`     | `performance_<setting>_plot.png` | Box-plots of the full docking-score distributions, highlighting the 10 lowest scores with a frequency colour-scale  |        |                                |
| `correlation.py`     | `correlationplot<num>/hexbinplot<num>.png`          | Pair-wise Tanimoto similarity vs. Score scatter & Pearson-r annotation                           |



