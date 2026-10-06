import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
files = sorted((ROOT / "Data" / "Raw" / "aeid_files").rglob("*.xlsx"))

for f in files:
    d = pd.read_excel(f, sheet_name="Results")
    d = d.iloc[3:].copy()   # first 3 rows only describe the columns
    d["HITC"] = pd.to_numeric(d["HITC"], errors="coerce")
    print("=====", f.name)
    print("rows:", len(d), "| unique IDs:", d["PUBCHEM_EXT_DATASOURCE_REGID"].nunique())
    print(d["PUBCHEM_ACTIVITY_OUTCOME"].value_counts(dropna=False).to_string())
    print("HITC missing:", d["HITC"].isna().sum(), "| HITC >= 0.9:", (d["HITC"] >= 0.9).sum())
    if "aeid786_" in f.name:
        print(d.head(5).to_string())