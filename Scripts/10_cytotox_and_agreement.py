import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
m = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)

A = ["aeid_786", "aeid_787", "aeid_789", "aeid_790"]
B = ["aeid_2053", "aeid_2054", "aeid_2119", "aeid_2120"]
for name, g in [("group A", A), ("group B", B)]:
    same = (m[g].notna().all(axis=1) | m[g].isna().all(axis=1)).all()
    print(name, "tested together as a block:", same)

pairs = {"ERa_BLA": (786, 787), "ERa_LUC_0.5nM": (789, 790),
         "ERa_LUC_0.1nM": (2053, 2054), "ERb_BLA": (2119, 2120)}
rows = []
for name, (a, v) in pairs.items():
    ant, via = m[f"aeid_{a}"], m[f"aeid_{v}"]
    ok = ant.notna() & via.notna()
    act = (ant == 1) & ok
    both = act & (via == 1)
    rows.append({"assay": name, "tested": int(ok.sum()),
                 "antagonist_active": int(act.sum()),
                 "also_viability_active": int(both.sum()),
                 "active_only_in_antagonist": int((act & (via == 0)).sum()),
                 "pct_possible_cytotox": round(both.sum() / max(act.sum(), 1) * 100, 1)})
cyto = pd.DataFrame(rows)
print(cyto.to_string(index=False))
cyto.to_csv(ROOT / "results" / "er_antagonist_cytotox_check.csv", index=False)

prim = m[["aeid_786", "aeid_789", "aeid_2053", "aeid_2119"]]
print("\nAgreement between primary endpoints (correlation):")
print(prim.corr().round(2).to_string())
print("Active in all 4 primary endpoints:", int((prim == 1).all(axis=1).sum()))
print("Active in at least 1 primary endpoint:", int((prim == 1).any(axis=1).sum()))

m.index.to_series().to_csv(ROOT / "Data" / "Processed" / "dtxsid_list.txt",
                           index=False, header=False)