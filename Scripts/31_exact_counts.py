"""Exact counts per test from the master table: tested, active, not active, blank."""
import sys
from pathlib import Path
for _p in Path(__file__).resolve().parents:   # find er_common.py even if the file is nested
    for _c in (_p, _p / "Scripts"):
        if (_c / "er_common.py").exists():
            sys.path.insert(0, str(_c))
import pandas as pd
from er_common import ROOT, load_master

m = load_master()
rows = []
for c in m.columns:
    n = int(m[c].notna().sum())
    a = int((m[c] == 1).sum())
    rows.append({"assay": c, "tested": n, "active": a, "not_active": n - a,
                 "blank": len(m) - n, "pct_blank": round(100 * (len(m) - n) / len(m), 1),
                 "pct_active": round(100 * a / n, 1) if n else None})
res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "exact_counts_per_test.csv", index=False)
print(res.to_string(index=False))
print("Chemicals:", len(m))