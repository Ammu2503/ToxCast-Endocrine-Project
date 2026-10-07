import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
f = pd.read_csv(ROOT / "results" / "smiles_failed.csv")
print("failed:", len(f))
blank = f["SMILES"].isna() | f["SMILES"].astype(str).str.strip().eq("")
print("no SMILES at all:", int(blank.sum()))
print("has SMILES but RDKit failed:", int((~blank).sum()))
print("NTP Mix entries:", int(f["PREFERRED_NAME"].str.contains("NTP Mix", na=False).sum()))
print(f.loc[~blank, ["DTXSID", "PREFERRED_NAME", "SMILES"]].to_string())