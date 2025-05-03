from pathlib import Path

# project_root is the folder that contains this script (adjust .parent level if needed)
project_root = Path(__file__).resolve().parent          # …/RadSmina/rad
# or: project_root = Path(__file__).resolve().parents[1]   if this file is nested deeper

# now define resources relative to that root
json_path   = project_root / "examples" / "zid_scores_rock1.json"
model_path  = project_root / "models"   / "checkpoint.pt"
output_dir  = project_root / "results"  / "run_01"

# use the paths normally
with json_path.open() as f:
    data = f.read()
