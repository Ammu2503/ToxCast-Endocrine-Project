import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
files = sorted((ROOT / "Data" / "Raw" / "aeid_files").rglob("*.xlsx"))

for f in files:
    xl = pd.ExcelFile(f)
    print("=====", f.name, xl.sheet_names)
    for s in xl.sheet_names:
        d = xl.parse(s)
        print(s, d.shape)
        print(d.columns.tolist())
        print(d.head(3).to_string())
    print()