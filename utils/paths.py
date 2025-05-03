from pathlib import Path

#  /rds/general/user/bl521/home/RadSmina_old/
project_root = Path(__file__).resolve().parent.parent

def examples(*parts):
    return project_root / "examples" / Path(*parts)

def scores(*parts):
    return project_root / "examples" / "scores"/ Path(*parts)
