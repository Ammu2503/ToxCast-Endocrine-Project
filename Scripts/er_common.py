"""Shared helpers for scripts 25 to 28. Keep this file in the Scripts folder."""
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import TruncatedSVD

ROOT = Path(__file__).resolve().parent.parent
N_COMP = 20      # fingerprint summary columns
SEED = 0


def load_master():
    return pd.read_csv(ROOT / "Data" / "Processed" / "er_master_all.csv", index_col=0)


def load_summary():
    return pd.read_csv(ROOT / "results" / "er_master_all_missingness.csv", index_col=0)


def load_fingerprint_summary():
    """Returns a table of N_COMP numbers per chemical (compressed Morgan fingerprint)."""
    fp = pd.read_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv").set_index("dtxsid")
    svd = TruncatedSVD(n_components=N_COMP, random_state=SEED)
    S = pd.DataFrame(svd.fit_transform(fp.values), index=fp.index,
                     columns=[f"fp{i}" for i in range(N_COMP)])
    return S


def final_assays(master, min_tested=500):
    """Reads results/final_assays.txt (one aeid_NNN per line, # for comments).
    If the file does not exist yet, falls back to every endpoint tested on at least
    min_tested chemicals, and says so. Replace this with your decided list."""
    path = ROOT / "results" / "final_assays.txt"
    if path.exists():
        names = [ln.split("#")[0].strip() for ln in path.read_text().splitlines()]
        names = [n for n in names if n]
        missing = [n for n in names if n not in master.columns]
        if missing:
            raise SystemExit(f"Not in master table: {missing}")
        return names
    print("WARNING: results/final_assays.txt not found. Using every endpoint with at least "
          f"{min_tested} chemicals tested. Create the file once you decide the assay list.")
    return [c for c in master.columns if master[c].notna().sum() >= min_tested]