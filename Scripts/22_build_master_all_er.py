import re
import zipfile
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ZIP = Path(r"C:\Users\harik\Downloads\pubchem_invitrodb_v4_3_19AUG2026.zip")
OUT_DIR = ROOT / "Data" / "Raw" / "aeid_files_all_er"
name_col = "assay_component_endpoint_name"

# 1. Endpoint list from script 21 (30 ER endpoints, Rev-ErbA and ch1/ch2 removed)
av = pd.read_csv(ROOT / "results" / "er_assay_availability.csv")
av.loc[av[name_col].str.startswith("ACEA_ER"), "mechanism_guess"] = "cell growth (guess)"
aeids = sorted(av["aeid"].tolist())

# 2. Extract only those files from the zip
OUT_DIR.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(ZIP) as z:
    for info in z.infolist():
        m = re.search(r"aeid(\d+)_", info.filename)
        if m and int(m.group(1)) in aeids and not info.is_dir():
            target = OUT_DIR / Path(info.filename).name
            if not target.exists():
                target.write_bytes(z.read(info))

# 3. Long table (same logic as scripts 08 and 09)
frames = []
for f in sorted(OUT_DIR.glob("*.xlsx")):
    aeid = int(re.search(r"aeid(\d+)_", f.name).group(1))
    d = pd.read_excel(f, sheet_name="Results")
    d = d[d["PUBCHEM_ACTIVITY_OUTCOME"].notna()].copy()
    d["aeid"] = aeid
    d["dtxsid"] = d["PUBCHEM_ACTIVITY_URL"].astype(str).str.extract(r"(DTXSID\d+)")
    d["outcome"] = d["PUBCHEM_ACTIVITY_OUTCOME"].astype(int)
    frames.append(d[["aeid", "dtxsid", "outcome"]])
    print("read", f.name, len(d), "rows")
long = pd.concat(frames, ignore_index=True)
long = long.dropna(subset=["dtxsid"])
long["active"] = (long["outcome"] == 2).astype(int)  # 1 and 3 = not active; any active sample = active

master = long.pivot_table(index="dtxsid", columns="aeid", values="active", aggfunc="max")
master.columns = [f"aeid_{c}" for c in master.columns]
master.to_csv(ROOT / "Data" / "Processed" / "er_master_all.csv")

# 4. Summary per endpoint, with category and cell line
tested = master.notna().sum()
active = master.sum()
info = av.set_index("aeid")
summary = pd.DataFrame({
    "endpoint": [info.loc[int(c.split("_")[1]), name_col] for c in master.columns],
    "category": [info.loc[int(c.split("_")[1]), "mechanism_guess"] for c in master.columns],
    "cell_line": [info.loc[int(c.split("_")[1]), "cell_short_name"] for c in master.columns],
    "signal_direction": [info.loc[int(c.split("_")[1]), "signal_direction"] for c in master.columns],
    "tested": tested.values,
    "missing_pct": (master.isna().mean() * 100).round(1).values,
    "active": active.astype(int).values,
    "active_pct_of_tested": (active / tested * 100).round(1).values,
}, index=master.columns)
summary.to_csv(ROOT / "results" / "er_master_all_missingness.csv")

print("\nchemicals:", len(master), "| endpoints:", master.shape[1])
print("overall missing %:", round(master.isna().values.mean() * 100, 1))
print(summary.to_string())

# 5. How many categories does each chemical have data for?
cat = {c: summary.loc[c, "category"] for c in master.columns}
per_cat = pd.DataFrame({k: master[[c for c in master.columns if cat[c] == k]].notna().any(axis=1)
                        for k in sorted(set(cat.values()))})
print("\nChemicals with data in each category:")
print(per_cat.sum().to_string())
print("\nNumber of categories tested per chemical:")
print(per_cat.sum(axis=1).value_counts().sort_index().to_string())