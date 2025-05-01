import json
import numpy as np
import matplotlib.pyplot as plt

# Global cutoff for "virtual active" definition
VIRTUAL_ACTIVE_SCORE_CUTOFF = -10.0

def collect_smina_virtual_actives(smina_json_paths, score_cutoff=-10.0):
    virtual_actives = set()

    for path in smina_json_paths:
        with open(path, 'r') as f:
            data = json.load(f)  # normal dict with insertion order in Py3.7+
        # each value is like [node_id, SMILES, docking_score]
        for mol_id, arr in data.items():
            score = arr[-1]  # last element is the docking score
            if score < score_cutoff:
                virtual_actives.add(mol_id)

    return virtual_actives

def load_traversal_order(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    traversal_list = list(data.keys())
    return traversal_list

def main():
    # 1) Gather the "Smina" virtual actives from multiple JSONs
    smina_jsons = [
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores0.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores1.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores2.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores3.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores4.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores5.json",
        "/rds/general/user/bl521/home/rad/examples/goldilocks_scores6.json"
    ]
    smina_virtual_actives = collect_smina_virtual_actives(smina_jsons, score_cutoff=-10.0)
    print(f"Number of Smina-based virtual actives: {len(smina_virtual_actives)}")

    # 2) For each method, load the json files
    rad_json = "/rds/general/user/bl521/home/rad/examples/m16ef400_scores_rock1.json"
    random_json = "/rds/general/user/bl521/home/rad/examples/random_10k_dock_results.json"

    rad_traversal = load_traversal_order(rad_json)
    rad_ids = set(rad_traversal)
    common_rad = smina_virtual_actives.intersection(rad_ids)
    print("Overlap with RAD traversal:", len(common_rad))

    random_traversal = load_traversal_order(random_json)
    random_ids = set(random_traversal)
    common_random = smina_virtual_actives.intersection(random_ids)
    print("Overlap with random10k traversal:", len(common_random))

    plot_mols_traversed = np.linspace(0, len(rad_ids), 50, dtype=int)
    
    rad_virtual_actives_recovered = []
    for n in plot_mols_traversed:
        # partial set: the first n molecules
        rad_subset = set(rad_traversal[:n])
        rad_recovered_count = len(rad_subset.intersection(smina_virtual_actives))
        rad_fraction = (rad_recovered_count / len(smina_virtual_actives))*100
        rad_virtual_actives_recovered.append(rad_fraction)
    print(f"Percentage of virtual hits retrieved by RAD: {rad_fraction:.2f}%")

    rnd_virtual_actives_recovered = []
    for n in plot_mols_traversed:
        # partial set: the first n molecules
        rnd_subset = set(random_traversal[:n])
        rnd_recovered_count = len(rnd_subset.intersection(smina_virtual_actives))
        rnd_fraction = (rnd_recovered_count / len(smina_virtual_actives))*100
        rnd_virtual_actives_recovered.append(rnd_fraction)
    print(f"Percentage of virtual hits retrieved by Random10k: {rnd_fraction:.2f}%")

    # 4) Plot
    plt.figure(figsize=(7,5))
    plt.plot(plot_mols_traversed, rad_virtual_actives_recovered, label="Smina+RAD10k", color='blue')
    plt.plot(plot_mols_traversed, rnd_virtual_actives_recovered, label="Random10k", color='green')

    # Add more lines for additional methods
    plt.xlabel("Number of Molecules Traversed")
    plt.ylabel("Percentage of Smina-based Virtual Actives Recovered")
    plt.title("Enrichment Plot for RAD of ROCK1 (score < -10.0 in Smina)")
    plt.ylim(0, 100)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("coverage_m16ef00.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
