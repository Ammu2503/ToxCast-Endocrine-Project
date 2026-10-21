import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
m = pd.read_csv(ROOT / "Data" / "Processed" / "er_master_all.csv", index_col=0)
summ = pd.read_csv(ROOT / "results" / "er_master_all_missingness.csv", index_col=0)
fp_ids = set(pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv",
                         usecols=["dtxsid"])["dtxsid"])

# Short labels, ordered by program so related endpoints sit together
summ["label"] = summ["endpoint"].str.replace("TOX21_", "T21_").str.replace("_ratio", "") \
    .str.replace("_Antagonist", "_Ant").str.replace("_Agonist", "_Ag") \
    .str.replace("_viability", "_via").str.replace("_VM7", "")
order = summ.sort_values(["category", "cell_line", "endpoint"]).index.tolist()
miss = m[order].isna()
labels = [summ.loc[c, "label"] for c in order]
miss.columns = labels

# 1. Bar chart: % missing per endpoint
fig, ax = plt.subplots(figsize=(10, 8))
pct = miss.mean() * 100
ax.barh(range(len(pct)), pct.values)
ax.set_yticks(range(len(pct))); ax.set_yticklabels(pct.index, fontsize=7)
ax.invert_yaxis(); ax.set_xlabel("% of chemicals never tested"); ax.set_xlim(0, 100)
ax.set_title("Missing data per endpoint")
fig.tight_layout(); fig.savefig(ROOT / "results" / "missing_bar.png", dpi=150); plt.close(fig)

# 2. Pattern heatmap: chemicals grouped by identical missing pattern
pat = miss.astype(int)
key = pat.astype(str).agg("".join, axis=1)
counts = key.value_counts()
print("Distinct missing patterns:", len(counts))
print("Top 10 patterns cover", round(counts.head(10).sum() / len(key) * 100, 1), "% of chemicals")
pd.DataFrame({"n_chemicals": counts.values, "pattern": counts.index}).head(30) \
    .to_csv(ROOT / "results" / "missing_patterns_top30.csv", index=False)
sorted_idx = key.map(counts).sort_values(ascending=False, kind="stable").index
fig, ax = plt.subplots(figsize=(10, 7))
ax.imshow(pat.loc[sorted_idx].values, aspect="auto", cmap="Greys", interpolation="nearest")
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=6)
ax.set_ylabel("Chemicals (grouped by pattern)"); ax.set_title("Missing (black) vs tested (white)")
fig.tight_layout(); fig.savefig(ROOT / "results" / "missing_heatmap.png", dpi=150); plt.close(fig)

# 3. Correlation of missingness indicators (plus: chemical has no structure)
ind = pat.copy()
ind["no_structure"] = (~ind.index.isin(fp_ids)).astype(int)
ind = ind.loc[:, ind.std() > 0]
corr = ind.corr()
corr.to_csv(ROOT / "results" / "missing_correlation.csv")
fig, ax = plt.subplots(figsize=(10, 9))
im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="coolwarm")
ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.columns, fontsize=6)
fig.colorbar(im, ax=ax, fraction=0.04); ax.set_title("Correlation of 'missing' between endpoints")
fig.tight_layout(); fig.savefig(ROOT / "results" / "missing_correlation.png", dpi=150); plt.close(fig)

print("Chemicals with no structure:", int(ind["no_structure"].sum()))
print("Saved: missing_bar.png, missing_heatmap.png, missing_correlation.png,")
print("       missing_patterns_top30.csv, missing_correlation.csv  (in results/)")