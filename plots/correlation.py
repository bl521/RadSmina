import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from itertools import combinations
from scipy.stats import pearsonr
from tqdm import tqdm
from utils.paths import scores_path

def load_json_file(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)

    entries = []
    for mol_id, arr in data.items():
        # array = [node_id, SMILES, docking_score]
        smiles = arr[1]
        score = arr[2]
        entries.append((smiles, score))

    return entries

def pairwise_similarity_vs_score_diff(entries, radius=2, nBits=1024):

    # Convert SMILES -> RDKit Molecule -> fingerprint
    fps = []
    scores = []
    print("Generating fingerprints...")
    for (smi, scr) in tqdm(entries, desc="Molecules"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            # If SMILES fails to parse, skip or handle differently
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        fps.append(fp)
        scores.append(scr)

    # Create arrays of (similarity, score_diff) across all pairs
    similarities = []
    score_diffs  = []

    n = len(fps)
    pair_count = n*(n-1)//2
    print(f"Kept {n} valid molecules after parsing SMILES.")

    # itertools.combinations => pairwise combos
    for i, j in tqdm(combinations(range(n), 2), total=pair_count, desc="Pairwise combos"):
        sim_ij = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        diff_ij = abs(scores[i] - scores[j])

        similarities.append(sim_ij)
        score_diffs.append(diff_ij)

    similarities = np.array(similarities)
    score_diffs  = np.array(score_diffs)

    # Compute Pearson correlation (r, p-value)
    corr, pval = pearsonr(similarities, score_diffs)

    plt.figure(figsize=(7,5))
    # plt.hexbin(
    #     similarities,
    #     score_diffs,
    #     gridsize=50,   # number of hex bins along each dimension
    #     cmap="viridis"   # color map
    # )
    sns.regplot(x=similarities, y=score_diffs, scatter_kws={'alpha':0.3}, line_kws={'color':'red'})
    # plt.colorbar(label="Count of Molecule Pairs")  # create color legend
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Absolute Score Difference (kcal/mol)")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Absolute Score Difference (kcal/mol)")
    plt.title(f"Pairwise Tanimoto vs. Score Diff\nr={corr:.3f}, p={pval:.2e}")
    plt.grid(alpha=0.3, linestyle="--")
    plt.tight_layout()
    plt.show()
    # plt.title("Hexbin Plot of Similarity vs. Score Difference")
    plt.savefig("correlationplot100.png", dpi=300)
    # plt.show()

    return corr, pval

def main():
    json_path = scores_path("m16ef400_scores_rock1.json") 
    entries = load_json_file(json_path)

    entries.sort(key=lambda x: x[1])  # x[1] is the score
    lowest_entries = entries[:100]

    # Calculate correlation
    corr, pval = pairwise_similarity_vs_score_diff(lowest_entries)
    print(f"Pairwise Tanimoto vs. |score_i - score_j| correlation = {corr:.3f}, p-value = {pval:.3g}")

if __name__ == "__main__":
    main()

