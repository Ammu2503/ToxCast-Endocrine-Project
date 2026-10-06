import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
f = ROOT / "Data" / "Raw" / "CCD-Batch-Search_2026-10-06_11_48_06.xlsx"

xl = pd.ExcelFile(f)
print(xl.sheet_names)
for s in xl.sheet_names:
    d = xl.parse(s)
    print("\nSHEET:", s, d.shape)
    print(d.columns.tolist())
    print(d.head(3).to_string())