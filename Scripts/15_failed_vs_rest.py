import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
m = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)
failed = pd.read_csv(ROOT / "results" / "smiles_failed.csv")["DTXSID"]
m["no_structure"] = m.index.isin(failed)
print(m["no_structure"].value_counts().rename({True: "no structure", False: "has structure"}))

cols = [c for c in m.columns if c != "no_structure"]
rows = []
for c in cols:
    for flag, g in m.groupby("no_structure"):
        tested = g[c].notna().sum()
        rows.append({"endpoint": c,
                     "group": "no_structure" if flag else "has_structure",
                     "tested_pct": round(tested / len(g) * 100, 1),
                     "active_pct_of_tested": round(g[c].sum() / tested * 100, 1) if tested else None})
res = pd.DataFrame(rows).pivot(index="endpoint", columns="group",
                               values=["tested_pct", "active_pct_of_tested"])
print(res.to_string())

prim = ["aeid_786", "aeid_789", "aeid_2053", "aeid_2119"]
for flag, g in m.groupby("no_structure"):
    tested = g[prim].notna().any(axis=1)
    anyact = (g[prim] == 1).any(axis=1)
    print("no structure" if flag else "has structure",
          "| active in at least 1 primary endpoint, % of chemicals tested in any:",
          round(anyact.sum() / tested.sum() * 100, 1))

res.to_csv(ROOT / "results" / "failed_structure_comparison.csv")