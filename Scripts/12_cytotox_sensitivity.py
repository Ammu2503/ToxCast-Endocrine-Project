import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
long = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_long.csv")

pairs = {"ERa_BLA": (786, 787), "ERa_LUC_0.5nM": (789, 790),
         "ERa_LUC_0.1nM": (2053, 2054), "ERb_BLA": (2119, 2120)}
rows = []
for fold in [1, 2, 3, 5, 10]:
    for name, (a, v) in pairs.items():
        ant = long[(long.aeid == a) & (long.outcome == 2)].groupby("dtxsid")["ac50"].min()
        via = long[(long.aeid == v) & (long.outcome == 2)].groupby("dtxsid")["ac50"].min()
        d = ant.to_frame("ant")
        d["via"] = via.reindex(d.index)
        has_via = d.index.isin(via.index)
        selective = has_via & (d["ant"] * fold <= d["via"])
        kept = (~has_via).sum() + selective.sum()
        rows.append({"assay": name, "fold": fold,
                     "kept_pct": round(kept / len(d) * 100, 1),
                     "missing_ac50": int(d["ant"].isna().sum())})
res = pd.DataFrame(rows).pivot(index="assay", columns="fold", values="kept_pct")
print(res.to_string())
print("antagonist hits with a missing AC50:",
      int(sum(r["missing_ac50"] for r in rows if r["fold"] == 1)))
res.to_csv(ROOT / "results" / "er_antagonist_cytotox_sensitivity.csv")