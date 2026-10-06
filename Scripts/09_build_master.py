import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
long = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_long.csv")

# Rule 1: active = outcome 2. Outcomes 1 and 3 count as not active.
long["active"] = (long["outcome"] == 2).astype(int)

o3 = long[long["outcome"] == 3]
print("outcome 3 rows:", len(o3), "| hitc below -0.5:", int((o3["hitc"] < -0.5).sum()))

# How often do repeat samples of one chemical disagree?
print("\nRepeat samples that disagree:")
for a, g in long.groupby("aeid"):
    s = g.groupby("dtxsid")["active"].agg(["min", "max", "count"])
    multi = s[s["count"] > 1]
    conflict = int((multi["min"] != multi["max"]).sum())
    print(a, "chemicals with repeats:", len(multi), "| disagree:", conflict)

# Rule 2: active if any sample is active. Blank means never tested.
master = long.pivot_table(index="dtxsid", columns="aeid", values="active", aggfunc="max")
master.columns = [f"aeid_{c}" for c in master.columns]
master.to_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv")

tested = master.notna().sum()
active = master.sum()
summary = pd.DataFrame({
    "tested": tested,
    "missing_pct": (master.isna().mean() * 100).round(1),
    "active": active.astype(int),
    "active_pct_of_tested": (active / tested * 100).round(1),
})
print("\nchemicals:", len(master), "| endpoints:", master.shape[1])
print(summary.to_string())
print("overall missing %:", round(master.isna().values.mean() * 100, 1))
print("chemicals tested in all 8 endpoints:", int(master.notna().all(axis=1).sum()))

(ROOT / "results").mkdir(exist_ok=True)
summary.to_csv(ROOT / "results" / "er_antagonist_missingness.csv")