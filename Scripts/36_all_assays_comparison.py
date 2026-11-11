"""Deletion versus MICE versus native blanks, with EVERY usable ER assay taken as the target in turn
(no single chosen target). Results are averaged over assays, as in Walter et al. and Kim et al.
Targets : every ER activity endpoint with at least MIN_CLASS active and MIN_CLASS not active
          chemicals that also have a fingerprint. Viability endpoints are not targets.
Inputs  : the core assay list (results/final_assays.txt) minus the target (and minus its
          same-assay sibling at another dose, SIBLINGS), with or without the Morgan fingerprint.
Handling: deletion | MICE | native blanks (gradient boosting).  Models: logistic regression,
          random forest, gradient boosting, XGBoost (if installed).
Scores  : MCC (main), average precision, ROC AUC, F1, balanced accuracy on test chemicals with
          complete inputs (fair head-to-head).  Splits: N_SPLITS random 80/20 splits per assay.
Writes results/all_assays_comparison.csv (per assay and split), all_assays_summary.csv (mean over
assays) and all_assays_paired.csv (MICE minus deletion, paired by assay).
Run: python 36_all_assays_comparison.py   (optional: add a number to set the splits, default 5)"""
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
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (average_precision_score, roc_auc_score, matthews_corrcoef,
                             f1_score, balanced_accuracy_score)
from er_common import ROOT, load_master, load_summary, load_fingerprint_summary, final_assays, SEED

warnings.filterwarnings("ignore")
N_SPLITS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
N_IMPUTE = 2
MIN_CLASS = 50
MIN_TRAIN_ROWS = 100
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


def models():
    out = {
        "logistic regression": make_pipeline(StandardScaler(), LogisticRegression(
            max_iter=2000, C=0.1, class_weight="balanced")),
        "random forest": RandomForestClassifier(
            300, class_weight="balanced_subsample", min_samples_leaf=3, n_jobs=-1,
            random_state=SEED),
        "gradient boosting": HistGradientBoostingClassifier(
            class_weight="balanced", random_state=SEED),
    }
    if XGBClassifier is not None:
        out["XGBoost"] = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1,
                                       subsample=0.8, colsample_bytree=0.5, n_jobs=-1,
                                       random_state=SEED, eval_metric="logloss")
    return out


def get(fs, A_, B_):
    return A_ if fs == "assays only" else pd.concat([A_, B_], axis=1)


rows = []


def add(tgt, split, fs, handling, algo, ntr, yte, p, prev):
    if len(set(yte)) < 2:
        return
    pred = (p >= np.quantile(p, 1 - prev)).astype(int)
    rows.append({"target": tgt, "split": split, "feature_set": fs, "handling": handling,
                 "algorithm": algo, "n_train": ntr, "n_test": len(yte),
                 "mcc": matthews_corrcoef(yte, pred), "f1": f1_score(yte, pred),
                 "balanced_acc": balanced_accuracy_score(yte, pred),
                 "avg_precision": average_precision_score(yte, p),
                 "roc_auc": roc_auc_score(yte, p)})


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
        if yc.nunique() < 2:
            continue
        prev = ytr.mean()
        if comp_tr.sum() >= MIN_TRAIN_ROWS and ytr[comp_tr].nunique() == 2:
            tr_c = comp_tr[comp_tr].index
            for fs in ("assays only", "assays + fingerprint"):
                for algo, mod in models().items():
                    mod.fit(get(fs, A.loc[tr_c], B.loc[tr_c]), ytr.loc[tr_c])
                    p = pd.Series(mod.predict_proba(
                        get(fs, A.loc[te][comp_te], B.loc[te][comp_te]))[:, 1], index=yc.index)
                    add(tgt, split, fs, "deletion", algo, len(tr_c), yc, p,
                        ytr.loc[tr_c].mean())
        acc = {}
        for k in range(N_IMPUTE):
            imp = IterativeImputer(max_iter=10, sample_posterior=True,
                                   random_state=100 * split + k)
            Xtr0 = pd.concat([A.loc[tr], Z.loc[tr]], axis=1)
            Xte0 = pd.concat([A.loc[te], Z.loc[te]], axis=1)
            Xtr_f = pd.DataFrame(imp.fit_transform(Xtr0), index=tr, columns=Xtr0.columns)
            Xte_f = pd.DataFrame(imp.transform(Xte0), index=te, columns=Xte0.columns)
            for fs in ("assays only", "assays + fingerprint"):
                for algo, mod in models().items():
                    mod.fit(get(fs, Xtr_f[assays], B.loc[tr]), ytr)
                    acc.setdefault((fs, algo), []).append(
                        mod.predict_proba(get(fs, Xte_f[assays], B.loc[te]))[:, 1])
        for (fs, algo), ps in acc.items():
            p = pd.Series(np.mean([pd.Series(x).rank().values for x in ps], axis=0), index=te)
            add(tgt, split, fs, "MICE", algo, len(tr), yc, p[comp_te], prev)
        for fs in ("assays only", "assays + fingerprint"):
            mod = HistGradientBoostingClassifier(class_weight="balanced", random_state=SEED)
            mod.fit(get(fs, A.loc[tr], B.loc[tr]), ytr)
            p = pd.Series(mod.predict_proba(get(fs, A.loc[te], B.loc[te]))[:, 1], index=te)
            add(tgt, split, fs, "native blanks", "gradient boosting", len(tr), yc,
                p[comp_te], prev)

res = pd.DataFrame(rows)
res.to_csv(ROOT / "results" / "all_assays_comparison.csv", index=False)
per_t = res.groupby(["feature_set", "handling", "algorithm", "target"])[
    ["mcc", "avg_precision", "n_train"]].mean().reset_index()
summary = per_t.groupby(["feature_set", "handling", "algorithm"]).agg(
    assays=("target", "nunique"), mcc_mean=("mcc", "mean"), mcc_sd_over_assays=("mcc", "std"),
    ap_mean=("avg_precision", "mean"), n_train_mean=("n_train", "mean")).round(3)
summary.to_csv(ROOT / "results" / "all_assays_summary.csv")

out = []
for (fs, algo), g in per_t.groupby(["feature_set", "algorithm"]):
    dele = g[g["handling"] == "deletion"].set_index("target")
    mi = g[g["handling"] == "MICE"].set_index("target")
    s = mi.index.intersection(dele.index)
    if len(s) < 3:
        continue
    for met in ("mcc", "avg_precision"):
        diff = mi.loc[s, met] - dele.loc[s, met]
        se = diff.std(ddof=1) / np.sqrt(len(diff))
        t = stats.t.ppf(0.975, len(diff) - 1)
        out.append({"feature_set": fs, "algorithm": algo, "metric": met, "assays": len(diff),
                    "deletion_mean": round(dele.loc[s, met].mean(), 3),
                    "mean_diff": round(diff.mean(), 4),
                    "ci_low": round(diff.mean() - t * se, 4),
                    "ci_high": round(diff.mean() + t * se, 4),
                    "assays_mice_better": int((diff > 0).sum())})
paired = pd.DataFrame(out)
paired.to_csv(ROOT / "results" / "all_assays_paired.csv", index=False)
pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 300)
print("\nMean over assays (test chemicals with complete inputs):")
print(summary.to_string())
if len(paired):
    print("\nMICE minus deletion, paired by assay:")
    print(paired.to_string(index=False))