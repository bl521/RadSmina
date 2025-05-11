from pathlib import Path

#  /rds/general/user/bl521/home/RadSmina/
project_root = Path(__file__).resolve().parent.parent

def radsmina_path(*parts):
    return project_root / "radsmina" / Path(*parts)

def scores_path(*parts):
    return project_root / "radsmina" / "scores"/ Path(*parts)
