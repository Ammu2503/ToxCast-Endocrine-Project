import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
long = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_long.csv")

pairs = {"ERa_BLA": (786, 787), "ERa_LUC_0.5nM": (789, 790),
         "ERa_LUC_0.1nM": (2053, 2054), "ERb_BLA": (2119, 2120)}
FOLD = 3   # antagonist AC50 must be at least 3 times lower than viability AC50

rows = []
for name, (a, v) in pairs.items():
    ant = long[(long.aeid == a) & (long.outcome == 2)].groupby("dtxsid")["ac50"].min()
    via_ids = set(long[(long.aeid == v) & (long.outcome == 2)]["dtxsid"])
    via_ac50 = long[(long.aeid == v) & (long.outcome == 2)].groupby("dtxsid")["ac50"].min()
    d = ant.to_frame("ant_ac50")
    d["via_active"] = d.index.isin(via_ids)
    d["via_ac50"] = via_ac50.reindex(d.index)
    clean = ~d["via_active"]
    selective = d["via_active"] & (d["ant_ac50"] * FOLD <= d["via_ac50"])
    cyto = d["via_active"] & ~selective
    rows.append({"assay": name, "antagonist_active": len(d),
                 "no_viability_hit": int(clean.sum()),
                 "viability_hit_but_selective": int(selective.sum()),
                 "likely_cytotoxic": int(cyto.sum()),
                 "kept_pct": round((clean.sum() + selective.sum()) / len(d) * 100, 1)})
out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(ROOT / "results" / "er_antagonist_cytotox_rule.csv", index=False)