import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.lines import Line2D

def load_rad_scores(rad_json_path):
    with open(rad_json_path, "r") as f:
        rad_data = json.load(f)
    return [entry[-1] for entry in rad_data.values()]

def load_smina_scores(json_paths):
    smina_scores = []
    for path in json_paths:
        with open(path, "r") as f:
            data = json.load(f)
        smina_scores.extend(entry[-1] for entry in data.values())
    return smina_scores

def load_random_scores(random_json_path):
    with open(random_json_path, "r") as f:
        random_data = json.load(f)
    return [entry[-1] for entry in random_data.values()]

def compare_boxplots(rad_json_path, smina_json_paths, random_json_path, lowest_n=10):
    """Boxplot for RAD vs Smina docking scores, highlighting the lowest few scores 
       with a color scale representing how many molecules share each score.
    """
    # 1) Load the scores
    rad_scores = load_rad_scores(rad_json_path)
    smina_scores = load_smina_scores(smina_json_paths)
    random_scores = load_random_scores(random_json_path)

    # 2) Find the lowest N
    rad_lowest = sorted(rad_scores)[:lowest_n]
    smina_lowest = sorted(smina_scores)[:lowest_n]
    random_lowest = sorted(random_scores)[:lowest_n]

    # 3) Count unique scores
    rad_unique_scores, rad_counts = np.unique(rad_lowest, return_counts=True)
    smina_unique_scores, smina_counts = np.unique(smina_lowest, return_counts=True)
    random_unique_scores, random_counts = np.unique(random_lowest, return_counts = True)

    # 4) Create a global colormap scale based on min/max count
    all_counts = np.concatenate([rad_counts, smina_counts, random_counts])
    min_count, max_count = all_counts.min(), all_counts.max()

    cmap = cm.Purples
    norm = mcolors.Normalize(vmin=min_count, vmax=max_count)

    # --- Create a figure and axes explicitly ---
    fig, ax = plt.subplots(figsize=(7, 5))

    # 5) Draw the boxplot 
    bp = ax.boxplot(
        [rad_scores, smina_scores, random_scores],
        medianprops=dict(color="black"),
        patch_artist=True,
        tick_labels=["RAD+Smina(10k)", "Smina(132k)", "Random(10k)"]
    )

    # Color the boxes differently
    box_colors = ["lightgreen", "lightcoral", "lightblue"]
    for box, color in zip(bp['boxes'], box_colors):
        box.set_facecolor(color)

    # 6) Scatter the lowest scores, shading by frequency
    #    - rad at x=0.9, smina at x=2.1
    for score, count in zip(rad_unique_scores, rad_counts):
        color = cmap(norm(count))
        ax.scatter(1.0, score, color=color, edgecolors='k', s=80, zorder=3)
    for score, count in zip(smina_unique_scores, smina_counts):
        color = cmap(norm(count))
        ax.scatter(2.0, score, color=color, edgecolors='k', s=80, zorder=3)
    for score, count in zip(random_unique_scores, random_counts):
        color = cmap(norm(count))
        ax.scatter(3.0, score, color=color, edgecolors='k', s=80, zorder=3)

    # 7) Add a colorbar, specifying which figure and axes to use
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Just a dummy array to let colorbar know the scale
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Number of molecules at this score")

    # 8) Legend:
    lowest_legend = Line2D([0],[0], marker='o', color='none',
                           markeredgecolor='k', label='Lowest scored molecules',
                           markersize=10, markerfacecolor='none')

    ax.legend(handles=[lowest_legend],
              loc='upper left',
              bbox_to_anchor=(0.01, 0.99))

    # Basic styling
    ax.set_title("Comparison of Docking Score Distributions")
    ax.set_ylabel("Docking Score (kcal/mol)")
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)

    fig.tight_layout()
    fig.savefig("performance_comparison_plot.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    rad_json = scores("m16ef400_scores_rock1.json")
    smina_jsons = [
        scores("goldilocks_scores0.json"),
        scores("goldilocks_scores1.json"),
        scores("goldilocks_scores2.json"),
        scores("goldilocks_scores3.json"),
        scores("goldilocks_scores4.json"),
        scores("goldilocks_scores5.json"),
        scores("goldilocks_scores6.json")
    ]
    random_json = scores("random_10k_scores.json")
    compare_boxplots(rad_json, smina_jsons, random_json, lowest_n=10)


