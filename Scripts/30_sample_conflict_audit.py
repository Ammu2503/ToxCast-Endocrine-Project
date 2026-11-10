"""How often do repeated samples of one chemical disagree inside one assay?
Reads the same files as script 22 (Data/Raw/aeid_files_all_er/*.xlsx, sheet Results).
Outcome codes: 2 = active, 1 = inactive, 3 = wrong direction (grouped with inactive).
Writes results/sample_conflict_audit.csv."""
import sys
from pathlib import Path
for _p in Path(__file__).resolve().parents:   # find er_common.py even if the file is nested
    for _c in (_p, _p / "Scripts"):
        if (_c / "er_common.py").exists():
            sys.path.insert(0, str(_c))
import re
import pandas as pd
from er_common import ROOT
IN_DIR = ROOT / "Data" / "Raw" / "aeid_files_all_er"
rows = []
for f in sorted(IN_DIR.glob("*.xlsx")):
    aeid = int(re.search(r"aeid(\d+)_", f.name).group(1))
    d = pd.read_excel(f, sheet_name="Results")
    d = d[d["PUBCHEM_ACTIVITY_OUTCOME"].notna()].copy()
    d["dtxsid"] = d["PUBCHEM_ACTIVITY_URL"].astype(str).str.extract(r"(DTXSID\d+)")[0]
    d["outcome"] = d["PUBCHEM_ACTIVITY_OUTCOME"].astype(int)
    d["active"] = (d["outcome"] == 2).astype(int)
    g = d.groupby("dtxsid")
    n = g["active"].size()
    k = g["active"].sum()
    inc = g["outcome"].apply(lambda s: (s == 3).any())
    multi = n > 1
    conflict = multi & (k > 0) & (k < n)
    rows.append({
        "assay": f"aeid_{aeid}", "chemicals": len(n), "rows": len(d),
        "multi_sample": int(multi.sum()),
        "pct_multi": round(100 * multi.mean(), 1),
        "conflict": int(conflict.sum()),
        "pct_conflict_of_multi": round(100 * conflict.sum() / max(multi.sum(), 1), 1),
        "any_wrong_direction": int(inc.sum()),
        "active_any": int((k > 0).sum()),
        "active_majority": int((k / n > 0.5).sum()),
        "active_all": int((k == n).sum()),
        "active_any_single_only": int(((k > 0) & ~multi).sum()),
        "flip_any_vs_majority": int(((k > 0) & ~(k / n > 0.5)).sum()),
    })
    print("done", aeid)

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "sample_conflict_audit.csv", index=False)
cols = ["assay", "chemicals", "pct_multi", "pct_conflict_of_multi",
        "any_wrong_direction", "active_any", "active_majority", "flip_any_vs_majority"]
print(res[cols].to_string(index=False))
tot = res[["chemicals", "conflict", "multi_sample"]].sum()
print("\nAll assays: chemical-assay pairs", int(tot["chemicals"]))
print("with several samples:", int(tot["multi_sample"]))
print("with conflicting samples:", int(tot["conflict"]))