# RAD + SMINA Pipeline README
A quick-start guide for Retrieval-Augmented Docking with SMINA

## 1  Overview
Retrieval-Augmented Docking (RAD) couples an HNSW similarity graph with an external docking engine so that only a small, high-priority subset of a large library is scored.
This repository adapts the original RAD workflow (which used DOCK 3.7) to SMINA, demonstrating that comparable enrichment can be achieved while docking far fewer molecules.
```text
          SMILES  ──▶  fingerprints
                            │
                            ▼
          ┌───────────────────────────────┐
          │   HNSW build   (M, ef_c)      │   ❶ run once
          └───────────────────────────────┘
                            ▲
                            │ greedy traversal ❷
 node_id ─▶ ligand_id ─▶ tmp .mol2 ─▶ SMINA ─▶ best score   ❸ score_fn

