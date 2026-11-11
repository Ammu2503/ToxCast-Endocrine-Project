"""80/20 hold-out with 5-fold cross-validation tuning and confidence intervals, for every usable ER assay.
Design (the usual one, as in Walter et al.):
  1. Split the chemicals 80/20 (stratified). The 20% hold-out is never used for fitting or tuning.
  2. On the 80% training part, blanks in the predictor assays are handled by
       deletion | MICE (fitted on the training part only) | native blanks (gradient boosting).
  3. Each model's settings are chosen by 5-fold cross-validation on the training part (scored by
     average precision). The 5-fold validation score is saved next to the hold-out score.
  4. The tuned model is scored once on the hold-out chemicals that have complete predictors
     (fair head-to-head with deletion); MICE and native blanks are also scored on ALL hold-out chemicals.
  5. 95% confidence intervals: bootstrap over hold-out chemicals (per assay) and a t interval over
     assays for means and for MICE minus deletion.
Note: MICE is fitted on the training part before the 5-fold tuning, so the tuning folds share one
imputer (it sees predictor values of the validation fold, never the labels). The hold-out is clean.
Writes results/holdout_cv_comparison.csv, holdout_cv_missingness.csv, holdout_cv_summary.csv,
holdout_cv_paired.csv.   Run: python 37_holdout_cv_ci.py  (optional: add a number = hold-out splits, default 1)"""
import os
os.environ["PYTHONWARNINGS"] = "ignore"
import sys
from pathlib import Path
for _p in Path(__file__).resolve().parents:   # find er_common.py even if the file is nested
    for _c in (_p, _p / "Scripts"):
        if (_c / "er_common.py").exists():
            sys.path.insert(0, str(_c))
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV
from sklearn.metrics import average_precision_score, matthews_corrcoef
from er_common import ROOT, load_master, load_summary, load_fingerprint_summary, final_assays, SEED

warnings.filterwarnings("ignore")
N_SPLITS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
N_BOOT = 1000
N_IMPUTE = 2
MIN_CLASS = 50
MIN_TRAIN_ROWS = 100
N_ITER = 6
SIBLINGS = {"aeid_789": ["aeid_2053"], "aeid_2053": ["aeid_789"]}
try:
    from xgboost import XGBClassifier
except ImportError:
    XGBClassifier = None
    print("xgboost not installed: skipping it. Install with: pip install xgboost")

m = load_master()
summ = load_summary()
S = load_fingerprint_summary()
fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")
fp.columns = [f"bit{i}" for i in range(fp.shape[1])]
core = final_assays(m)
is_via = summ["endpoint"].str.lower().str.contains("viability")
has_fp = m.index.intersection(fp.index)
targets = []
for c in m.columns:
    if is_via.get(c, False):
        continue
    y = m.loc[has_fp, c].dropna()
    if (y == 1).sum() >= MIN_CLASS and (y == 0).sum() >= MIN_CLASS:
        targets.append(c)
print(f"{len(targets)} target assays: {targets}", flush=True)

CV = StratifiedKFold(5, shuffle=True, random_state=SEED)


def searchers(native=False):
    """Untuned estimator and its settings grid, for each algorithm."""
    if native:
        return {"gradient boosting": (HistGradientBoostingClassifier(
            class_weight="balanced", random_state=SEED),
            {"learning_rate": [0.03, 0.1], "max_depth": [None, 4, 8],
             "min_samples_leaf": [20, 50], "l2_regularization": [0, 1]})}
    out = {
        "logistic regression": (make_pipeline(StandardScaler(), LogisticRegression(
            max_iter=2000, class_weight="balanced")),
            {"logisticregression__C": [0.01, 0.1, 1, 10]}),
        "random forest": (RandomForestClassifier(
            class_weight="balanced_subsample", n_jobs=1, random_state=SEED),
            {"n_estimators": [200, 400], "max_depth": [None, 10, 20],
             "min_samples_leaf": [1, 3, 5], "max_features": ["sqrt", 0.2]}),
    }
    if XGBClassifier is not None:
        out["XGBoost"] = (XGBClassifier(n_jobs=1, random_state=SEED, eval_metric="logloss"),
                          {"n_estimators": [200, 400], "max_depth": [3, 4, 6],
                           "learning_rate": [0.03, 0.1], "subsample": [0.7, 0.9],
                           "colsample_bytree": [0.3, 0.5]})
    return out


def tune(est, grid, X, y):
    s = RandomizedSearchCV(est, grid, n_iter=N_ITER, cv=CV, scoring="average_precision",
                           random_state=SEED, n_jobs=-1, refit=True)
    s.fit(X, y)
    sd = s.cv_results_["std_test_score"][s.best_index_]
    return s, s.best_score_, sd


def get(fs, A_, B_):
    return A_ if fs == "assays only" else pd.concat([A_, B_], axis=1)


rng = np.random.default_rng(SEED)


def boot(y, pred, p):
    y = np.asarray(y)
    mc, ap = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max():
            continue
        mc.append(matthews_corrcoef(y[i], pred[i]))
        ap.append(average_precision_score(y[i], p[i]))
    return np.percentile(mc, [2.5, 97.5]), np.percentile(ap, [2.5, 97.5])


rows, miss = [], []


def add(tgt, split, fs, handling, algo, test_set, ntr, yte, p, prev, cv_ap, cv_sd, params):
    if len(set(yte)) < 2:
        return
    p = np.asarray(p)
    pred = (p >= np.quantile(p, 1 - prev)).astype(int)
    (m_lo, m_hi), (a_lo, a_hi) = boot(yte, pred, p)
    rows.append({"target": tgt, "split": split, "feature_set": fs, "handling": handling,
                 "algorithm": algo, "test_set": test_set, "n_train": ntr, "n_test": len(yte),
                 "mcc": matthews_corrcoef(yte, pred), "mcc_lo": m_lo, "mcc_hi": m_hi,
                 "avg_precision": average_precision_score(yte, p), "ap_lo": a_lo, "ap_hi": a_hi,
                 "cv5_avg_precision": cv_ap, "cv5_sd": cv_sd, "settings": str(params)})


for ti, tgt in enumerate(targets):
    assays = [a for a in core if a != tgt and a not in SIBLINGS.get(tgt, [])]
    assays = [a for a in assays if m[a].notna().any()]
    d = m[[tgt] + assays].join(S, how="inner").join(fp, how="inner")
    d = d[d[tgt].notna()]
    y_all = d[tgt].astype(int)
    A, Z, B = d[assays], d[list(S.columns)], d[list(fp.columns)]
    print(f"[{ti + 1}/{len(targets)}] {tgt}: {len(d)} chemicals, {int(y_all.sum())} active, "
          f"{len(assays)} predictor assays", flush=True)
    for split in range(N_SPLITS):
        tr, te = train_test_split(d.index, test_size=0.2, stratify=y_all, random_state=split)
        ytr, yte = y_all.loc[tr], y_all.loc[te]
        comp_tr = A.loc[tr].notna().all(axis=1)
        comp_te = A.loc[te].notna().all(axis=1)
        yc = yte[comp_te]
        miss.append({"target": tgt, "split": split, "chemicals": len(d), "active": int(y_all.sum()),
                     "predictor_assays": len(assays),
                     "pct_predictor_cells_blank": round(100 * A.isna().values.mean(), 1),
                     "pct_chemicals_with_any_blank": round(100 * A.isna().any(axis=1).mean(), 1),
                     "train_chemicals": len(tr), "train_complete_deletion": int(comp_tr.sum()),
                     "test_chemicals": len(te), "test_complete": int(comp_te.sum())})
        if yc.nunique() < 2:
            continue
        prev = ytr.mean()
        # deletion
        if comp_tr.sum() >= MIN_TRAIN_ROWS and ytr[comp_tr].sum() >= 25 \
                and (ytr[comp_tr] == 0).sum() >= 25:
            tr_c = comp_tr[comp_tr].index
            for fs in ("assays only", "assays + fingerprint"):
                for algo, (est, grid) in searchers().items():
                    s, cvs, cvsd = tune(est, grid, get(fs, A.loc[tr_c], B.loc[tr_c]), ytr.loc[tr_c])
                    p = s.predict_proba(get(fs, A.loc[te][comp_te], B.loc[te][comp_te]))[:, 1]
                    add(tgt, split, fs, "deletion", algo, "complete-assay chemicals", len(tr_c), yc,
                        p, ytr.loc[tr_c].mean(), cvs, cvsd, s.best_params_)
        # MICE (fitted on the training part only; tuned on the first imputation)
        acc = {}
        for k in range(N_IMPUTE):
            imp = IterativeImputer(max_iter=10, sample_posterior=True, random_state=100 * split + k)
            Xtr0 = pd.concat([A.loc[tr], Z.loc[tr]], axis=1)
            Xte0 = pd.concat([A.loc[te], Z.loc[te]], axis=1)
            Xtr_f = pd.DataFrame(imp.fit_transform(Xtr0), index=tr, columns=Xtr0.columns)
            Xte_f = pd.DataFrame(imp.transform(Xte0), index=te, columns=Xte0.columns)
            for fs in ("assays only", "assays + fingerprint"):
                for algo, (est, grid) in searchers().items():
                    Xt, Xe = get(fs, Xtr_f[assays], B.loc[tr]), get(fs, Xte_f[assays], B.loc[te])
                    if k == 0:
                        s, cvs, cvsd = tune(est, grid, Xt, ytr)
                        acc[(fs, algo)] = {"ps": [], "cv": (cvs, cvsd), "best": s.best_estimator_,
                                           "params": s.best_params_}
                        acc[(fs, algo)]["ps"].append(s.predict_proba(Xe)[:, 1])
                    else:
                        mod = clone(acc[(fs, algo)]["best"]).fit(Xt, ytr)
                        acc[(fs, algo)]["ps"].append(mod.predict_proba(Xe)[:, 1])
        for (fs, algo), v in acc.items():
            p = pd.Series(np.mean([pd.Series(x).rank().values for x in v["ps"]], axis=0), index=te)
            add(tgt, split, fs, "MICE", algo, "all test chemicals", len(tr), yte, p, prev,
                v["cv"][0], v["cv"][1], v["params"])
            add(tgt, split, fs, "MICE", algo, "complete-assay chemicals", len(tr), yc, p[comp_te],
                prev, v["cv"][0], v["cv"][1], v["params"])
        # native blanks
        for fs in ("assays only", "assays + fingerprint"):
            for algo, (est, grid) in searchers(native=True).items():
                s, cvs, cvsd = tune(est, grid, get(fs, A.loc[tr], B.loc[tr]), ytr)
                p = pd.Series(s.predict_proba(get(fs, A.loc[te], B.loc[te]))[:, 1], index=te)
                add(tgt, split, fs, "native blanks", algo, "all test chemicals", len(tr), yte, p,
                    prev, cvs, cvsd, s.best_params_)
                add(tgt, split, fs, "native blanks", algo, "complete-assay chemicals", len(tr), yc,
                    p[comp_te], prev, cvs, cvsd, s.best_params_)
        print(f"   split {split + 1} done", flush=True)

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "holdout_cv_comparison.csv", index=False)
pd.DataFrame(miss).to_csv(ROOT / "results" / "holdout_cv_missingness.csv", index=False)
c = res[res["test_set"] == "complete-assay chemicals"]
per_t = c.groupby(["feature_set", "handling", "algorithm", "target"])[
    ["mcc", "avg_precision", "cv5_avg_precision", "n_train"]].mean().reset_index()


def tci(x):
    x = np.asarray(x, float)
    se = x.std(ddof=1) / np.sqrt(len(x))
    t = stats.t.ppf(0.975, len(x) - 1)
    return x.mean() - t * se, x.mean() + t * se


srows = []
for (fs, h, algo), g in per_t.groupby(["feature_set", "handling", "algorithm"]):
    lo, hi = tci(g["mcc"]) if len(g) > 2 else (np.nan, np.nan)
    srows.append({"feature_set": fs, "handling": h, "algorithm": algo, "assays": len(g),
                  "mcc_mean": g["mcc"].mean(), "mcc_ci_low": lo, "mcc_ci_high": hi,
                  "holdout_ap_mean": g["avg_precision"].mean(),
                  "cv5_ap_mean": g["cv5_avg_precision"].mean(), "n_train_mean": g["n_train"].mean()})
summary = pd.DataFrame(srows)
summary = summary.round(3)
summary.to_csv(ROOT / "results" / "holdout_cv_summary.csv", index=False)
out = []
for (fs, algo), g in per_t.groupby(["feature_set", "algorithm"]):
    dele = g[g["handling"] == "deletion"].set_index("target")
    mi = g[g["handling"] == "MICE"].set_index("target")
    s = mi.index.intersection(dele.index)
    if len(s) < 3:
        continue
    for met in ("mcc", "avg_precision"):
        diff = mi.loc[s, met] - dele.loc[s, met]
        lo, hi = tci(diff)
        out.append({"feature_set": fs, "algorithm": algo, "metric": met, "assays": len(diff),
                    "deletion_mean": round(dele.loc[s, met].mean(), 3),
                    "mean_diff": round(diff.mean(), 4), "ci_low": round(lo, 4),
                    "ci_high": round(hi, 4), "assays_mice_better": int((diff > 0).sum())})
paired = pd.DataFrame(out)
paired.to_csv(ROOT / "results" / "holdout_cv_paired.csv", index=False)
pd.set_option("display.width", 230)
pd.set_option("display.max_rows", 300)
print("\nMissingness per assay:")
print(pd.DataFrame(miss).drop(columns="split").to_string(index=False))
print("\nMean over assays, hold-out test chemicals with complete predictors (95% t interval over assays):")
print(summary.to_string(index=False))
if len(paired):
    print("\nMICE minus deletion, paired by assay:")
    print(paired.to_string(index=False))