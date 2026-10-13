import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
m = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)
fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")

A = ["aeid_786", "aeid_787", "aeid_789", "aeid_790"]
B = ["aeid_2053", "aeid_2054"]
used_cols = A + B

u = m[used_cols].copy()
u["missing_groupA"] = u[A].isna().all(axis=1)
u["missing_groupB"] = u[B].isna().all(axis=1)
u["has_structure"] = u.index.isin(fp.index)

tested_any = u[used_cols].notna().any(axis=1)
keep = tested_any & u["has_structure"]

print("Starting chemicals:", len(m))
print("Dropped, no result in the 6 ERa endpoints:", int((~tested_any).sum()))
print("Dropped, no structure:", int((~u["has_structure"]).sum()))
print("Dropped, either reason:", int((~keep).sum()))
print("Kept chemicals:", int(keep.sum()), "| endpoints:", len(used_cols))
print("Percent missing in kept table:",
      round(u.loc[keep, used_cols].isna().mean().mean() * 100, 1))

u.loc[keep].to_csv(ROOT / "Data" / "Processed" / "er_antagonist_master_used.csv")
u.loc[~keep].to_csv(ROOT / "results" / "master_used_dropped.csv")