import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
zip_path = Path(r"C:\Users\harik\Downloads\pubchem_invitrodb_v4_3_19AUG2026.zip")
out = ROOT / "Data" / "Raw" / "aeid_files"
out.mkdir(parents=True, exist_ok=True)

aeids = [786, 787, 789, 790, 2053, 2054, 2119, 2120]

with zipfile.ZipFile(zip_path) as z:
    for name in z.namelist():
        base = name.split("/")[-1]
        for a in aeids:
            if base.startswith(f"invitrodb_v4.3_aeid{a}_for_pubchem"):
                z.extract(name, out)
                print("extracted", base)