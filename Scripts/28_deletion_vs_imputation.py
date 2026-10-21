"""Main test: does filling blanks with MICE help a prediction model?
Same split, same model type, three ways of handling blanks:
  1. deletion   : train only on chemicals with no blank in the inputs
  2. MICE       : fill training blanks (imputer fitted on training only), 3 imputations
  3. native NaN : a model that accepts blanks directly (extra control)
Test set is untouched until the end. Change TARGET only after agreeing with Dr. Kattoum."""
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict
from sklearn.metrics import (average_precision_score, recall_score, f1_score,
                             balanced_accuracy_score)
from er_common import ROOT, load_master, load_fingerprint_summary, final_assays, SEED

TARGET = "aeid_786"     # TO BE AGREED: the one result the model predicts
N_SPLITS = 10           # repeated random 80/20 splits
N_IMPUTE = 3            # MICE completions per split
MIN_TRAIN_ROWS = 100    # deletion arm needs at least this many complete training rows

m = load_master()
assays = final_assays(m)
if TARGET not in assays:
    raise SystemExit(f"{TARGET} is not in the final assay list")
features = [a for a in assays if a != TARGET]
S = load_fingerprint_summary()
d = m[[TARGET] + features].join(S, how="inner")
d = d[d[TARGET].notna()]                       # only chemicals with a known target
y_all = d[TARGET].astype(int)
X_all = d[features + list(S.columns)]
print(f"Target {TARGET}: {len(d)} chemicals, {int(y_all.sum())} active "
      f"({y_all.mean() * 100:.1f}%).")
print(f"Features: {len(features)} assays + {S.shape[1]} structure columns")


def best_cut(y, p):
    """Cut-off chosen on TRAINING predictions only: the value that gives the best balanced
    accuracy. Same rule for every arm, so the arms are compared fairly."""
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


rows = []
for split in range(N_SPLITS):
    tr, te = train_test_split(d.index, test_size=0.2, stratify=y_all, random_state=split)
    Xtr, Xte, ytr, yte = X_all.loc[tr], X_all.loc[te], y_all.loc[tr], y_all.loc[te]
    complete_tr = Xtr.notna().all(axis=1)
    complete_te = Xte.notna().all(axis=1)

    # 1. Deletion
    ok_tr = complete_tr.sum() >= MIN_TRAIN_ROWS and ytr[complete_tr].nunique() == 2
    ok_te = complete_te.sum() > 10 and yte[complete_te].nunique() == 2
    if ok_tr and ok_te:
        cut_del = best_cut(ytr[complete_tr], cv_probs(rf(), Xtr[complete_tr], ytr[complete_tr]))
        p = rf().fit(Xtr[complete_tr], ytr[complete_tr]).predict_proba(Xte[complete_te])[:, 1]
        rows.append({"split": split, "method": "deletion", "test_set": "complete rows only",
                     "n_train": int(complete_tr.sum()), "n_test": int(complete_te.sum()),
                     **metrics(yte[complete_te], p, cut_del)})
    else:
        rows.append({"split": split, "method": "deletion", "test_set": "not possible",
                     "n_train": int(complete_tr.sum()), "n_test": int(complete_te.sum())})

    # 2. MICE (imputer fitted on training chemicals only)
    ps_all, ps_cv = [], []
    for k in range(N_IMPUTE):
        seed = 100 * split + k
        imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=seed)
        Xtr_f = pd.DataFrame(imp.fit_transform(Xtr), index=Xtr.index, columns=Xtr.columns)
        Xte_f = pd.DataFrame(imp.transform(Xte), index=Xte.index, columns=Xte.columns)
        model = rf().fit(Xtr_f, ytr)
        ps_all.append(model.predict_proba(Xte_f)[:, 1])
        ps_cv.append(cv_probs(rf(), Xtr_f, ytr))
    p_all = pd.Series(np.mean(ps_all, axis=0), index=Xte.index)
    cut_mice = best_cut(ytr, np.mean(ps_cv, axis=0))
    rows.append({"split": split, "method": "MICE", "test_set": "all test rows",
                 "n_train": len(Xtr), "n_test": len(Xte), **metrics(yte, p_all, cut_mice)})
    if complete_te.sum() > 10 and yte[complete_te].nunique() == 2:
        rows.append({"split": split, "method": "MICE", "test_set": "complete rows only",
                     "n_train": len(Xtr), "n_test": int(complete_te.sum()),
                     **metrics(yte[complete_te], p_all[complete_te], cut_mice)})

    # 3. Native handling of blanks
    hgb = HistGradientBoostingClassifier(class_weight="balanced", random_state=SEED).fit(Xtr, ytr)
    p_nat = pd.Series(hgb.predict_proba(Xte)[:, 1], index=Xte.index)
    hgb2 = HistGradientBoostingClassifier(class_weight="balanced", random_state=SEED)
    cut_nat = best_cut(ytr, cv_probs(hgb2, Xtr, ytr))
    rows.append({"split": split, "method": "native NaN", "test_set": "all test rows",
                 "n_train": len(Xtr), "n_test": len(Xte), **metrics(yte, p_nat, cut_nat)})
    print("split", split + 1, "of", N_SPLITS, "done")

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "deletion_vs_imputation.csv", index=False)
cols = ["avg_precision", "recall", "f1", "balanced_acc"]
grp = res.groupby(["method", "test_set"])[["n_train", "n_test"] + cols]
summary = grp.agg(["mean", "std"]).round(3)
print(summary.to_string())
print("\nCut-offs were chosen on training data only, the same way for each method.")
print("Fair head-to-head: compare 'deletion / complete rows only' with "
      "'MICE / complete rows only'. Improvement counts only if it is larger than the std.")