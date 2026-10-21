"""Evidence for the MICE conditions in the article (section 4.5), for the final assay list.
1. Overlap: how many chemicals were tested in both assays of every pair.
2. Informative clues: how strongly the assays agree with each other (phi correlation).
3. Class balance: how rare are active results.
4. Convergence: does the imputed table stop changing as rounds increase.
"""
import warnings
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.exceptions import ConvergenceWarning
from er_common import ROOT, load_master, load_fingerprint_summary, final_assays, SEED

MIN_OVERLAP = 200          # pairs sharing fewer chemicals cannot teach MICE a link
m = load_master()
assays = final_assays(m)
S = load_fingerprint_summary()
d = m[assays].join(S, how="inner")
print(f"{len(assays)} assays | {len(d)} chemicals with a structure")

# 1. Overlap
obs = d[assays].notna().astype(int)
overlap = obs.T @ obs
overlap.to_csv(ROOT / "results" / "mice_check_overlap.csv")
low = [(a, b, int(overlap.loc[a, b])) for i, a in enumerate(assays) for b in assays[i + 1:]
       if overlap.loc[a, b] < MIN_OVERLAP]
print(f"\n1. OVERLAP: {len(low)} of {len(assays) * (len(assays) - 1) // 2} pairs share fewer than "
      f"{MIN_OVERLAP} chemicals")
for a, b, n in low[:40]:
    print(f"   {a} & {b}: {n}")
print("   Pairs with 0 overlap cannot inform each other at all.")

# 2. Agreement between assays
phi = d[assays].corr(min_periods=MIN_OVERLAP)
phi.to_csv(ROOT / "results" / "mice_check_correlation.csv")
vals = phi.where(~np.eye(len(assays), dtype=bool)).stack()
print("\n2. CLUES: correlation between assays (observed pairs only)")
print(f"   median {vals.median():.2f}, strongest {vals.max():.2f}, "
      f"pairs above 0.3: {(vals.abs() > 0.3).sum() // 2}")

# 3. Class balance
bal = pd.DataFrame({"tested": d[assays].notna().sum(),
                    "active": d[assays].sum().astype(int)})
bal["active_pct"] = (bal["active"] / bal["tested"] * 100).round(1)
bal.to_csv(ROOT / "results" / "mice_check_class_balance.csv")
print("\n3. CLASS BALANCE")
print(bal.to_string())

# 4. Convergence: change in imputed values between rounds
print("\n4. CONVERGENCE (mean absolute change in the filled-in cells vs the previous setting)")
blank = d.isna().values
prev, rows = None, []
for it in [1, 2, 3, 5, 8, 12, 20]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", ConvergenceWarning)
        imp = IterativeImputer(max_iter=it, tol=1e-3, random_state=SEED)
        out = imp.fit_transform(d)
    filled = out[blank]
    change = np.nan if prev is None else float(np.mean(np.abs(filled - prev)))
    rows.append({"max_iter": it, "mean_change_vs_previous": change})
    prev = filled
conv = pd.DataFrame(rows)
conv.to_csv(ROOT / "results" / "mice_check_convergence.csv", index=False)
print(conv.to_string(index=False))
print("   If the change is near 0 by the last rows, the loop has settled.")