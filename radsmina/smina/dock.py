import subprocess

def dock_with_smina(ligand, receptor, rec_lig):

    docked_output = ligand.replace(".mol2", "_docked.mol2")

    cmd_dock = [
        "smina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--out", docked_output,
        "--autobox_ligand", rec_lig,
        "--autobox_add", "4",
        "--exhaustiveness", "8",
        "--num_modes", "1"
    ]
    docking_process = subprocess.run(cmd_dock, capture_output=True, text=True)
    if docking_process.returncode != 0:
        raise RuntimeError(f"Smina docking failed:\n{docking_process.stderr}")

    output_text = docking_process.stdout
    # print("Full docking output:\n", docking_process.stdout)

    best_score = None
    for line in output_text.splitlines():
        if line.strip().startswith("1 "):
            parts = line.split()
            try:
                best_score = float(parts[1])
                break
            except:
                pass
    # Read back the docked .mol2 file that smina generated
    with open(docked_output, 'r') as f:
        mol2_text = f.read()

    return best_score, mol2_text


