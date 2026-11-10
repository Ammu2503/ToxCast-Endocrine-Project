"""Final comparison: feature set x way of handling blanks x algorithm, on identical splits.
Feature sets : fingerprint only | assays only | assays + fingerprint
Blanks       : none needed (fingerprint only) | deletion | MICE | native (HGB only)
Algorithms   : logistic regression, random forest, gradient boosting, XGBoost (if installed),
               k-nearest neighbours, SVM, neural network
Scores       : MCC (main), average precision, ROC AUC, F1, balanced accuracy.
               MCC, F1 and balanced accuracy use a cut that calls the top share of chemicals
               active, where the share equals the TRAINING active rate (no test labels used).
Test set stays untouched. Imputer is fitted on training chemicals only.
Writes results/algorithm_comparison.csv and results/algorithm_comparison_summary.csv."""
import sys
from pathlib import Path
for _p in Path(__file__).resolve().parents:   # find er_common.py even if the file is nested
    for _c in (_p, _p / "Scripts"):
        if (_c / "er_common.py").exists():
            sys.path.insert(0, str(_c))
import warnings
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (average_precision_score, roc_auc_score, matthews_corrcoef,
                             f1_score, balanced_accuracy_score)
from er_common import ROOT, load_master, load_fingerprint_summary, final_assays, SEED

warnings.filterwarnings("ignore")
TARGET = "aeid_786"   # TO BE AGREED with Dr. Kattoum
N_SPLITS = 10         # repeated 80/20 splits (Ronia suggested about 20 if time allows)
SLOW_SPLITS = 5       # SVM and neural network run only on the first SLOW_SPLITS splits
N_IMPUTE = 2
MIN_TRAIN_ROWS = 100

m = load_master()
assays = [a for a in final_assays(m) if a != TARGET]
S = load_fingerprint_summary()
fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")
d = m[[TARGET] + assays].join(S, how="inner").join(fp, how="inner")
d = d[d[TARGET].notna()]
y_all = d[TARGET].astype(int)
A = d[assays]                    # assay results (with blanks)
Z = d[list(S.columns)]           # 20 structure summary columns, only used as MICE clues
B = d[list(fp.columns)]          # raw fingerprint bits
print(f"Target {TARGET}: {len(d)} chemicals, {int(y_all.sum())} active "
      f"({100 * y_all.mean():.1f}%). {len(assays)} assay features, {B.shape[1]} bits.")


try:
    from xgboost import XGBClassifier
except ImportError:
    XGBClassifier = None
    print("xgboost not installed: skipping it. Install with: pip install xgboost")


def models(split=0):
    out = {
        "logistic regression": make_pipeline(
            StandardScaler(), LogisticRegression(max_iter=2000, C=0.1,
                                                 class_weight="balanced")),
        "random forest": RandomForestClassifier(
            300, class_weight="balanced_subsample", min_samples_leaf=3,
            n_jobs=-1, random_state=SEED),
        "gradient boosting": HistGradientBoostingClassifier(
            class_weight="balanced", random_state=SEED),
        "k-nearest neighbours": make_pipeline(
            StandardScaler(), KNeighborsClassifier(25, weights="distance")),
    }
    if XGBClassifier is not None:
        out["XGBoost"] = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1,
                                       subsample=0.8, colsample_bytree=0.5, n_jobs=-1,
                                       random_state=SEED, eval_metric="logloss")
    if split < SLOW_SPLITS:
        out["SVM"] = make_pipeline(StandardScaler(), SVC(class_weight="balanced"))
        out["neural network"] = make_pipeline(
            StandardScaler(), MLPClassifier((64,), early_stopping=True,
                                            max_iter=300, random_state=SEED))
    return out


def score(model, X):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    return model.decision_function(X)


def get(fs, parts):
    """Build the input table for a feature set from the assay part and the fingerprint."""
    if fs == "fingerprint only":
        return parts["B"]
    if fs == "assays only":
        return parts["A"]
    return pd.concat([parts["A"], parts["B"]], axis=1)


rows = []


def add(split, fs, handling, algo, test_set, ntr, yte, p, prev):
    if len(set(yte)) < 2:
        return
    pred = (p >= np.quantile(p, 1 - prev)).astype(int)
    rows.append({"split": split, "feature_set": fs, "handling": handling,
                 "algorithm": algo, "test_set": test_set, "n_train": ntr,
                 "n_test": len(yte),
                 "mcc": matthews_corrcoef(yte, pred), "f1": f1_score(yte, pred),
                 "balanced_acc": balanced_accuracy_score(yte, pred),
                 "avg_precision": average_precision_score(yte, p),
                 "roc_auc": roc_auc_score(yte, p)})


for split in range(N_SPLITS):
    tr, te = train_test_split(d.index, test_size=0.2, stratify=y_all, random_state=split)
    ytr, yte = y_all.loc[tr], y_all.loc[te]
    comp_tr = A.loc[tr].notna().all(axis=1)
    comp_te = A.loc[te].notna().all(axis=1)
    yc = yte[comp_te]
    prev = ytr.mean()

    # 1. fingerprint only (no blanks to handle)
    for algo, mod in models(split).items():
        mod.fit(B.loc[tr], ytr)
        p = pd.Series(score(mod, B.loc[te]), index=te)
        add(split, "fingerprint only", "none", algo, "all test chemicals", len(tr), yte, p,
            prev)
        add(split, "fingerprint only", "none", algo, "complete-assay chemicals", len(tr),
            yc, p[comp_te], prev)

    # 2. deletion
    ok = comp_tr.sum() >= MIN_TRAIN_ROWS and ytr[comp_tr].nunique() == 2
    if ok:
        tr_c = comp_tr[comp_tr].index
        for fs in ("assays only", "assays + fingerprint"):
            parts_tr = {"A": A.loc[tr_c], "B": B.loc[tr_c]}
            parts_te = {"A": A.loc[te][comp_te], "B": B.loc[te][comp_te]}
            for algo, mod in models(split).items():
                mod.fit(get(fs, parts_tr), ytr.loc[tr_c])
                p = pd.Series(score(mod, get(fs, parts_te)), index=yc.index)
                add(split, fs, "deletion", algo, "complete-assay chemicals",
                    len(tr_c), yc, p, ytr.loc[tr_c].mean())

    # 3. MICE: imputer fitted on training chemicals only, results averaged over completions
    acc = {}
    for k in range(N_IMPUTE):
        seed = 100 * split + k
        imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=seed)
        Xtr0 = pd.concat([A.loc[tr], Z.loc[tr]], axis=1)
        Xte0 = pd.concat([A.loc[te], Z.loc[te]], axis=1)
        Xtr_f = pd.DataFrame(imp.fit_transform(Xtr0), index=tr, columns=Xtr0.columns)
        Xte_f = pd.DataFrame(imp.transform(Xte0), index=te, columns=Xte0.columns)
        for fs in ("assays only", "assays + fingerprint"):
            parts_tr = {"A": Xtr_f[assays], "B": B.loc[tr]}
            parts_te = {"A": Xte_f[assays], "B": B.loc[te]}
            for algo, mod in models(split).items():
                mod.fit(get(fs, parts_tr), ytr)
                acc.setdefault((fs, algo), []).append(score(mod, get(fs, parts_te)))
    for (fs, algo), ps in acc.items():
        p = pd.Series(np.mean([pd.Series(x).rank().values for x in ps], axis=0), index=te)
        add(split, fs, "MICE", algo, "all test chemicals", len(tr), yte, p, prev)
        add(split, fs, "MICE", algo, "complete-assay chemicals", len(tr), yc, p[comp_te],
            prev)

    # 4. native handling of blanks (gradient boosting only)
    for fs in ("assays only", "assays + fingerprint"):
        mod = HistGradientBoostingClassifier(class_weight="balanced", random_state=SEED)
        mod.fit(get(fs, {"A": A.loc[tr], "B": B.loc[tr]}), ytr)
        p = pd.Series(score(mod, get(fs, {"A": A.loc[te], "B": B.loc[te]})), index=te)
        add(split, fs, "native blanks", "gradient boosting", "all test chemicals",
            len(tr), yte, p, prev)
        add(split, fs, "native blanks", "gradient boosting", "complete-assay chemicals",
            len(tr), yc, p[comp_te], prev)
    print("split", split + 1, "of", N_SPLITS, "done", flush=True)

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "algorithm_comparison.csv", index=False)
keys = ["feature_set", "handling", "algorithm", "test_set"]
vals = ["n_train", "mcc", "avg_precision", "roc_auc", "f1", "balanced_acc"]
summ = res.groupby(keys)[vals].agg(["mean", "std"]).round(3)
summ.to_csv(ROOT / "results" / "algorithm_comparison_summary.csv")
pd.set_option("display.width", 200)
print(summ.to_string())
print("\nBest MCC on complete-assay chemicals (the fair head-to-head):")
if "complete-assay chemicals" in summ.index.get_level_values("test_set"):
    fair = summ.xs("complete-assay chemicals", level="test_set")
    top = fair[("mcc", "mean")].sort_values(ascending=False)
    print(top.head(10).to_string())
else:
    print("No complete-assay test chemicals in these splits.")