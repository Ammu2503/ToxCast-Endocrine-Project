import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold

ROOT = Path(__file__).resolve().parent.parent
m = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)
A = ["aeid_786", "aeid_787", "aeid_789", "aeid_790"]
B = ["aeid_2053", "aeid_2054", "aeid_2119", "aeid_2120"]
m["missA"] = m[A].isna().all(axis=1)
m["missB"] = m[B].isna().all(axis=1)

print("missing group A only:", int((m.missA & ~m.missB).sum()))
print("missing group B only:", int((m.missB & ~m.missA).sum()))
print("missing both groups:", int((m.missA & m.missB).sum()))
print("tested in both:", int((~m.missA & ~m.missB).sum()))

# 1. Do chemicals missing one group look different in the other group?
both = m[~m.missA & ~m.missB]
for name, g, cols in [("missing A, compared on group B", m[m.missA], B),
                      ("missing B, compared on group A", m[m.missB], A)]:
    print("\n" + name)
    for c in cols:
        print(c, "active %  missing-group:", round(g[c].mean() * 100, 1),
              "| tested in both:", round(both[c].mean() * 100, 1))

# 2. Can chemical structure predict which group is missing?
fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")
d = m.join(fp, how="inner")
X = d[fp.columns]
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
rows = []
print("\nChemicals with a fingerprint:", len(d))
for label, y in [("missing group A", d["missA"]), ("missing group B", d["missB"])]:
    auc = cross_val_score(LogisticRegression(max_iter=2000), X, y.astype(int),
                          cv=cv, scoring="roc_auc")
    print(label, "AUC from fingerprints:", round(auc.mean(), 2), "+/-", round(auc.std(), 2))
    rows.append({"target": label, "auc_mean": round(auc.mean(), 3), "auc_std": round(auc.std(), 3)})
pd.DataFrame(rows).to_csv(ROOT / "results" / "missingness_check.csv", index=False)