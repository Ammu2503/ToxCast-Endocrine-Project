"""Does a chemical's structure predict whether an endpoint is blank for it?
AUC near 0.5 = structure cannot tell (fits random gaps). AUC well above 0.5 = the gaps
depend on what the chemical is. Replaces script 17 for the 30 endpoint table."""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from er_common import ROOT, load_master, load_summary, load_fingerprint_summary, SEED

MIN_PER_CLASS = 100          # skip endpoints with fewer tested or fewer blank chemicals

m = load_master()
summ = load_summary()
S = load_fingerprint_summary()
common = m.index.intersection(S.index)
print(f"Chemicals in master: {len(m)} | with a structure: {len(common)} | "
      f"without: {len(m) - len(common)}")
m, S = m.loc[common], S.loc[common]
X = (S - S.mean()) / S.std().replace(0, 1)

rows = []
for c in m.columns:
    y = m[c].isna().astype(int)
    n_blank, n_tested = int(y.sum()), int((1 - y).sum())
    if min(n_blank, n_tested) < MIN_PER_CLASS:
        rows.append({"endpoint": c, "name": summ.loc[c, "endpoint"], "n_blank": n_blank,
                     "n_tested": n_tested, "auc": np.nan, "note": "too few in one group"})
        continue
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    p = cross_val_predict(LogisticRegression(max_iter=2000), X, y, cv=cv,
                          method="predict_proba")[:, 1]
    rows.append({"endpoint": c, "name": summ.loc[c, "endpoint"], "n_blank": n_blank,
                 "n_tested": n_tested, "auc": round(roc_auc_score(y, p), 3), "note": ""})

res = pd.DataFrame(rows).sort_values("auc", ascending=False)
res.to_csv(ROOT / "results" / "missingness_vs_structure.csv", index=False)
print(res.to_string(index=False))
print("\nHow to read: 0.50 = no link to structure. Above about 0.65 = clear link. "
      "A clear link means NOT missing completely at random.")