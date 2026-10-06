import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
files = sorted((ROOT / "Data" / "Raw" / "aeid_files").rglob("*.xlsx"))

frames = []
for f in files:
    aeid = int(f.name.split("aeid")[1].split("_")[0])
    d = pd.read_excel(f, sheet_name="Results")
    d = d[d["PUBCHEM_ACTIVITY_OUTCOME"].notna()].copy()
    d["aeid"] = aeid
    d["sample_id"] = d["PUBCHEM_EXT_DATASOURCE_REGID"]
    d["dtxsid"] = d["PUBCHEM_ACTIVITY_URL"].astype(str).str.extract(r"(DTXSID\d+)")
    d["outcome"] = d["PUBCHEM_ACTIVITY_OUTCOME"].astype(int)
    d["hitc"] = pd.to_numeric(d["HITC"], errors="coerce")
    d["ac50"] = pd.to_numeric(d["AC50"], errors="coerce")
    frames.append(d[["aeid", "sample_id", "dtxsid", "outcome", "hitc", "ac50"]])

long = pd.concat(frames, ignore_index=True)
out = ROOT / "Data" / "Processed"
out.mkdir(parents=True, exist_ok=True)
long.to_csv(out / "er_antagonist_long.csv", index=False)

print("total rows:", len(long))
print("rows with no DTXSID:", long["dtxsid"].isna().sum())
for a, g in long.groupby("aeid"):
    print(a, "rows", len(g), "| unique DTXSID", g["dtxsid"].nunique(),
          "| DTXSID with more than 1 row", int((g["dtxsid"].value_counts() > 1).sum()))
print(long.groupby("outcome")["hitc"].describe().to_string())
print("unique DTXSID across all endpoints:", long["dtxsid"].nunique())
