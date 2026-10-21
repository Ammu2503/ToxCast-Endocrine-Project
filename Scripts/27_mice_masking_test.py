"""Masking test. Hide some KNOWN results, let MICE guess them, compare with the truth.
Baselines: (a) call everything inactive, (b) structure only. Final assay list is read
from results/final_assays.txt (see er_common.py). Writes results/mice_masking_test.csv."""
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, recall_score, f1_score,
                             balanced_accuracy_score, mean_absolute_error)
from er_common import ROOT, load_master, load_fingerprint_summary, final_assays, SEED

M = 5                 # MICE runs with different random seeds
MASK_FRAC = 0.2       # share of each assay's known results that is hidden
MIN_ACTIVE = 30       # skip assays with fewer known active results than this

m = load_master()
assays = final_assays(m)
S = load_fingerprint_summary()
d = m[assays].join(S, how="inner")
rng = np.random.RandomState(SEED)

masked = d.copy()
hidden = {}
for c in assays:
    known = d.index[d[c].notna()]
    if d.loc[known, c].sum() < MIN_ACTIVE:
        continue
    ids = rng.choice(known, size=int(len(known) * MASK_FRAC), replace=False)
    hidden[c] = ids
    masked.loc[ids, c] = np.nan
print(f"Testing {len(hidden)} assays. Skipped (too few actives): "
      f"{[a for a in assays if a not in hidden]}")

fills = []
for s in range(M):
    imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=s)
    out = pd.DataFrame(imp.fit_transform(masked), index=masked.index, columns=masked.columns)
    fills.append(out)
    print("MICE run", s + 1, "of", M, "done")
mice = sum(fills) / M

rows = []
for c, ids in hidden.items():
    y = d.loc[ids, c].astype(int)
    train = masked.index[masked[c].notna()]
    prev = masked.loc[train, c].mean()
    clf = LogisticRegression(max_iter=2000).fit(S.loc[train], masked.loc[train, c].astype(int))
    struct = pd.Series(clf.predict_proba(S.loc[ids])[:, 1], index=ids)
    scores = {"MICE": mice.loc[ids, c].clip(0, 1), "structure only": struct}
    for name, sc in scores.items():
        cut = np.quantile(sc, 1 - prev)
        pred = (sc >= cut).astype(int)
        rows.append({"assay": c, "method": name, "n_hidden": len(ids),
                     "active_rate": round(y.mean(), 3),
                     "mae": round(mean_absolute_error(y, sc), 3),
                     "avg_precision": round(average_precision_score(y, sc), 3),
                     "recall": round(recall_score(y, pred), 3),
                     "f1": round(f1_score(y, pred, zero_division=0), 3),
                     "balanced_acc": round(balanced_accuracy_score(y, pred), 3)})
    rows.append({"assay": c, "method": "all inactive", "n_hidden": len(ids),
                 "active_rate": round(y.mean(), 3), "mae": round(y.mean(), 3),
                 "avg_precision": round(y.mean(), 3), "recall": 0.0, "f1": 0.0,
                 "balanced_acc": 0.5})

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "mice_masking_test.csv", index=False)
print(res.to_string(index=False))
print("\nMean by method (lower mae is better; higher for the others):")
print(res.groupby("method")[["mae", "avg_precision", "recall", "f1", "balanced_acc"]]
      .mean().round(3))
print("\nMICE earns its place only if it beats BOTH baselines.")