"""Writes down which assays each analysis used, so the runs can be reproduced.
- results/assays_masking_test.txt : endpoints with at least 500 chemicals tested
  (the rule script 27 falls back to when final_assays.txt does not exist)
- results/assays_prediction_core.txt : the list in final_assays.txt (scripts 28 and 32)
Also checks the first list against results/mice_masking_test.csv."""
import sys
from pathlib import Path
for _p in Path(__file__).resolve().parents:   # find er_common.py even if the file is nested
    for _c in (_p, _p / "Scripts"):
        if (_c / "er_common.py").exists():
            sys.path.insert(0, str(_c))
import pandas as pd
from er_common import ROOT, load_master, final_assays

m = load_master()
rule = [c for c in m.columns if m[c].notna().sum() >= 500]
(ROOT / "results" / "assays_masking_test.txt").write_text(
    "# endpoints with at least 500 chemicals tested\n" + "\n".join(rule) + "\n")
core = final_assays(m)
(ROOT / "results" / "assays_prediction_core.txt").write_text(
    "# list used by scripts 28 and 32 (from final_assays.txt)\n" + "\n".join(core) + "\n")
print(len(rule), "assays meet the 500-chemical rule;", len(core), "in the core list")

f = ROOT / "results" / "mice_masking_test.csv"
if f.exists():
    used = set(pd.read_csv(f)["assay"])
    print("Assays in the masking result file:", len(used))
    print("In rule list but not in result file:", sorted(set(rule) - used))
    print("In result file but not in rule list:", sorted(used - set(rule)))