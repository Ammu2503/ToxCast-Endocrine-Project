"""Does imputation start to help as more assays are added?
Runs deletion vs MICE vs native-blank model for growing assay sets, same target, same splits.
Writes results/assay_set_sweep.csv, assay_set_sweep_summary.csv and assay_set_sweep.png.
Change TARGET only after agreeing with Dr. Kattoum."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict
from sklearn.metrics import (average_precision_score, recall_score, f1_score,
                             balanced_accuracy_score)
from er_common import ROOT, load_master, load_fingerprint_summary, SEED

TARGET = "aeid_786"      # TO BE AGREED
N_SPLITS = 5             # repeated random 80/20 splits
N_IMPUTE = 2             # MICE completions per split
MIN_TRAIN_ROWS = 100

# Assay sets, each one larger than the one before. The target is removed automatically.
# The viability endpoint of the target assay (aeid_787) is left out on purpose, because it
# comes from the same plate as the target.
CORE = ["aeid_785", "aeid_786", "aeid_788", "aeid_789", "aeid_2053"]
ATG = ["aeid_75", "aeid_117"]
ACEA = ["aeid_2"]
OT = ["aeid_742", "aeid_743", "aeid_744", "aeid_745", "aeid_746", "aeid_747",
      "aeid_750", "aeid_751"]
BINDING = ["aeid_708", "aeid_714", "aeid_725"]
SETS = {
    "1 core": CORE,
    "2 +ATG": CORE + ATG,
    "3 +ACEA": CORE + ATG + ACEA,
    "4 +OT": CORE + ATG + ACEA + OT,
    "5 +binding": CORE + ATG + ACEA + OT + BINDING,
}

m = load_master()
S = load_fingerprint_summary()


def best_cut(y, p):
    cuts = np.unique(np.quantile(p, np.linspace(0.05, 0.95, 91)))
    scores = [balanced_accuracy_score(y, (p >= c).astype(int)) for c in cuts]
    return cuts[int(np.argmax(scores))]


def metrics(y, p, cut):
    pred = (p >= cut).astype(int)
    return {"avg_precision": average_precision_score(y, p),
            "recall": recall_score(y, pred, zero_division=0),
            "f1": f1_score(y, pred, zero_division=0),
            "balanced_acc": balanced_accuracy_score(y, pred)}


def cv_probs(model, X, y):
    cv = StratifiedKFold(3, shuffle=True, random_state=SEED)
    return cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]


def rf():
    return RandomForestClassifier(300, class_weight="balanced_subsample",
                                  min_samples_leaf=3, n_jobs=-1, random_state=SEED)


def hgb():
    return HistGradientBoostingClassifier(class_weight="balanced", random_state=SEED)


rows = []
for set_name, cols in SETS.items():
    features = [c for c in cols if c != TARGET and c in m.columns]
    d = m[[TARGET] + features].join(S, how="inner")
    d = d[d[TARGET].notna()]
    y_all = d[TARGET].astype(int)
    X_all = d[features + list(S.columns)]
    n_complete = int(X_all.notna().all(axis=1).sum())
    print(f"\n{set_name}: {len(features)} assays, {len(d)} chemicals, "
          f"{n_complete} with no blanks ({n_complete / len(d) * 100:.1f}%)")
    for split in range(N_SPLITS):
        tr, te = train_test_split(d.index, test_size=0.2, stratify=y_all, random_state=split)
        Xtr, Xte, ytr, yte = X_all.loc[tr], X_all.loc[te], y_all.loc[tr], y_all.loc[te]
        ctr = Xtr.notna().all(axis=1)
        cte = Xte.notna().all(axis=1)
        base = {"set": set_name, "n_assays": len(features), "split": split}

        # 1. Deletion
        if ctr.sum() >= MIN_TRAIN_ROWS and ytr[ctr].nunique() == 2 and yte[cte].nunique() == 2:
            cut = best_cut(ytr[ctr], cv_probs(rf(), Xtr[ctr], ytr[ctr]))
            p = rf().fit(Xtr[ctr], ytr[ctr]).predict_proba(Xte[cte])[:, 1]
            rows.append({**base, "method": "deletion", "test_set": "complete rows",
                         "n_train": int(ctr.sum()), "n_test": int(cte.sum()),
                         **metrics(yte[cte], p, cut)})
        else:
            rows.append({**base, "method": "deletion", "test_set": "not possible",
                         "n_train": int(ctr.sum()), "n_test": int(cte.sum())})

        # 2. MICE, imputer fitted on training chemicals only
        ps, pcv = [], []
        for k in range(N_IMPUTE):
            imp = IterativeImputer(max_iter=10, sample_posterior=True,
                                   random_state=100 * split + k)
            Xtr_f = pd.DataFrame(imp.fit_transform(Xtr), index=Xtr.index, columns=Xtr.columns)
            Xte_f = pd.DataFrame(imp.transform(Xte), index=Xte.index, columns=Xte.columns)
            ps.append(rf().fit(Xtr_f, ytr).predict_proba(Xte_f)[:, 1])
            pcv.append(cv_probs(rf(), Xtr_f, ytr))
        p_all = pd.Series(np.mean(ps, axis=0), index=Xte.index)
        cut = best_cut(ytr, np.mean(pcv, axis=0))
        rows.append({**base, "method": "MICE", "test_set": "all rows",
                     "n_train": len(Xtr), "n_test": len(Xte), **metrics(yte, p_all, cut)})
        if cte.sum() > 10 and yte[cte].nunique() == 2:
            rows.append({**base, "method": "MICE", "test_set": "complete rows",
                         "n_train": len(Xtr), "n_test": int(cte.sum()),
                         **metrics(yte[cte], p_all[cte], cut)})

        # 3. Native handling of blanks
        p_nat = pd.Series(hgb().fit(Xtr, ytr).predict_proba(Xte)[:, 1], index=Xte.index)
        cut = best_cut(ytr, cv_probs(hgb(), Xtr, ytr))
        rows.append({**base, "method": "native NaN", "test_set": "all rows",
                     "n_train": len(Xtr), "n_test": len(Xte), **metrics(yte, p_nat, cut)})
        print("  split", split + 1, "of", N_SPLITS, "done")

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "assay_set_sweep.csv", index=False)
cols = ["n_train", "n_test", "avg_precision", "recall", "f1", "balanced_acc"]
summary = res.groupby(["set", "method", "test_set"])[cols].agg(["mean", "std"]).round(3)
summary.to_csv(ROOT / "results" / "assay_set_sweep_summary.csv")
print("\n", summary.to_string())

# Figure: left = ranking quality, right = chemicals deletion can still train on
sets = list(SETS)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
lines = [("deletion", "complete rows"), ("MICE", "complete rows"),
         ("MICE", "all rows"), ("native NaN", "all rows")]
for meth, ts in lines:
    sub = res[(res.method == meth) & (res.test_set == ts)]
    g = sub.groupby("set")["avg_precision"].agg(["mean", "std"]).reindex(sets)
    ax1.errorbar(range(len(sets)), g["mean"], yerr=g["std"], marker="o", capsize=3,
                 label=f"{meth} ({ts})")
ax1.set_xticks(range(len(sets)))
ax1.set_xticklabels(sets, rotation=20)
ax1.set_ylabel("Average precision (higher is better)")
ax1.set_title("Prediction quality by assay set")
ax1.legend(fontsize=8)
dele = res[(res.method == "deletion")].groupby("set")["n_train"].mean().reindex(sets)
ax2.bar(range(len(sets)), dele.values)
ax2.set_xticks(range(len(sets)))
ax2.set_xticklabels(sets, rotation=20)
ax2.set_ylabel("Training chemicals kept by deletion")
ax2.set_title("What deletion keeps")
fig.tight_layout()
fig.savefig(ROOT / "results" / "assay_set_sweep.png", dpi=150)
plt.close(fig)
print("Saved assay_set_sweep.csv, assay_set_sweep_summary.csv, assay_set_sweep.png")