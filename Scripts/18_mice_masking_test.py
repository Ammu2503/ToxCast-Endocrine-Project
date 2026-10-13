import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, recall_score, balanced_accuracy_score

ROOT = Path(__file__).resolve().parent.parent
GROUPS = {
    "A": ["aeid_786", "aeid_787", "aeid_789", "aeid_790"],
    "B": ["aeid_2053", "aeid_2054", "aeid_2119", "aeid_2120"],
}
N_COMP = 20      # fingerprint summary columns given to MICE
M = 5            # number of MICE runs (different random seeds)
MASK_FRAC = 0.2  # share of fully tested chemicals whose group is hidden
SEED = 0

m = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)
fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")
d = m.join(fp, how="inner")
assays = GROUPS["A"] + GROUPS["B"]
fpcols = list(fp.columns)

# 20 summary columns from the 2048 fingerprint bits
svd = TruncatedSVD(n_components=N_COMP, random_state=SEED)
S = pd.DataFrame(svd.fit_transform(d[fpcols].values), index=d.index,
                 columns=[f"fp{i}" for i in range(N_COMP)])
full = d[assays].join(S)

rng = np.random.RandomState(SEED)
complete = d[assays].dropna().index        # chemicals tested in all 8
rows = []

for hide_name, hide_cols in GROUPS.items():
    test_ids = rng.choice(complete, size=int(len(complete) * MASK_FRAC), replace=False)
    masked = full.copy()
    masked.loc[test_ids, hide_cols] = np.nan   # hide the whole group, like real data
    train_ids = masked.index.difference(test_ids)

    # MICE: M runs with different seeds, average the filled-in values
    fills = []
    for s in range(M):
        imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=s)
        out = pd.DataFrame(imp.fit_transform(masked), index=masked.index, columns=masked.columns)
        fills.append(out.loc[test_ids, hide_cols])
    mice = sum(fills) / M

    # Baseline 1: always say inactive
    # Baseline 2: structure only (logistic regression on the 20 fingerprint columns)
    for c in hide_cols:
        y = full.loc[test_ids, c].astype(int)
        known = full.loc[train_ids, c].dropna().astype(int)
        prev = known.mean()
        clf = LogisticRegression(max_iter=2000).fit(S.loc[known.index], known)
        struct = pd.Series(clf.predict_proba(S.loc[test_ids])[:, 1], index=test_ids)

        scores = {"MICE": mice[c], "structure only": struct}
        for name, sc in scores.items():
            cut = np.quantile(sc, 1 - prev)    # call the top 'prevalence' share active
            pred = (sc >= cut).astype(int)
            rows.append({"hidden group": hide_name, "endpoint": c, "method": name,
                         "n_test": len(test_ids), "active_rate_test": round(y.mean(), 3),
                         "avg_precision": round(average_precision_score(y, sc), 3),
                         "recall": round(recall_score(y, pred), 3),
                         "balanced_acc": round(balanced_accuracy_score(y, pred), 3)})
        rows.append({"hidden group": hide_name, "endpoint": c, "method": "all inactive",
                     "n_test": len(test_ids), "active_rate_test": round(y.mean(), 3),
                     "avg_precision": round(y.mean(), 3), "recall": 0.0, "balanced_acc": 0.5})

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "mice_masking_test.csv", index=False)
print(res.to_string(index=False))
print("\nMean by method:")
print(res.groupby("method")[["avg_precision", "recall", "balanced_acc"]].mean().round(3))