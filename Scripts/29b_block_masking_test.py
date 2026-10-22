"""Block masking test. For each screening program, hide ALL of that program's results
for 20 percent of the chemicals it tested, let MICE guess them from the OTHER programs
plus structure, and compare with the truth. This copies how real gaps look (whole
programs missing), so it is a harder and fairer test than hiding random cells.
Writes results/block_masking_test.csv."""
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, recall_score, balanced_accuracy_score
from er_common import ROOT, load_master, load_fingerprint_summary, SEED

M = 3
HIDE_FRAC = 0.2
MIN_ACTIVE = 30
PROGRAMS = {
    "Tox21": [785, 786, 787, 788, 789, 790, 1852, 2053, 2054, 2115, 2116, 2119, 2120],
    "ATG": [75, 117],
    "ACEA": [2],
    "OT": [742, 743, 744, 745, 746, 747, 750, 751],
    "NVS binding": [708, 714, 725],
}

m = load_master()
S = load_fingerprint_summary()
allc = [f"aeid_{a}" for p in PROGRAMS.values() for a in p]
cols = [c for c in allc if c in m.columns]
d = m[cols].join(S, how="inner")
rng = np.random.RandomState(SEED)
rows = []

for prog, aeids in PROGRAMS.items():
    pc = [f"aeid_{a}" for a in aeids if f"aeid_{a}" in d.columns]
    tested = d.index[d[pc].notna().any(axis=1)]
    if len(pc) == 0 or len(tested) < 200:
        continue
    ids = rng.choice(tested, size=int(len(tested) * HIDE_FRAC), replace=False)
    masked = d.copy()
    masked.loc[ids, pc] = np.nan
    fills = []
    for s in range(M):
        imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=s)
        arr = imp.fit_transform(masked)
        out = pd.DataFrame(arr, index=masked.index, columns=masked.columns)
        fills.append(out)
    mice = sum(fills) / M
    print(prog, "done:", len(ids), "chemicals hidden")
    for c in pc:
        truth = d.loc[ids, c].dropna()
        if len(truth) < 30 or truth.sum() < MIN_ACTIVE // 3:
            continue
        y = truth.astype(int)
        tr = masked.index[masked[c].notna()]
        prev = masked.loc[tr, c].mean()
        clf = LogisticRegression(max_iter=2000).fit(S.loc[tr], masked.loc[tr, c].astype(int))
        struct = pd.Series(clf.predict_proba(S.loc[y.index])[:, 1], index=y.index)
        for name, sc in {"MICE": mice.loc[y.index, c].clip(0, 1),
                         "structure only": struct}.items():
            pred = (sc >= np.quantile(sc, 1 - prev)).astype(int)
            rows.append({"program": prog, "assay": c, "method": name, "n_hidden": len(y),
                         "active_rate": round(y.mean(), 3),
                         "avg_precision": round(average_precision_score(y, sc), 3),
                         "recall": round(recall_score(y, pred), 3),
                         "balanced_acc": round(balanced_accuracy_score(y, pred), 3)})
        rows.append({"program": prog, "assay": c, "method": "all inactive",
                     "n_hidden": len(y), "active_rate": round(y.mean(), 3),
                     "avg_precision": round(y.mean(), 3), "recall": 0.0,
                     "balanced_acc": 0.5})

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "block_masking_test.csv", index=False)
print(res.groupby(["program", "method"])[["avg_precision", "recall", "balanced_acc"]]
      .mean().round(3))
print("\nIf MICE no longer beats structure only here,")
print("the random-cell masking test was too easy.")