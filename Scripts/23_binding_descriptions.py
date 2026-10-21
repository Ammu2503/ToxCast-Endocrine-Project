"""Prints the General sheet (assay description) for the binding endpoints and any others
you list. Use it to check species, method and cell system against my guesses."""
import re
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FOLDER = ROOT / "Data" / "Raw" / "aeid_files_all_er"
AEIDS = [708, 714, 725]      # add more aeids here to check them too

for f in sorted(FOLDER.glob("*.xlsx")):
    aeid = int(re.search(r"aeid(\d+)_", f.name).group(1))
    if aeid not in AEIDS:
        continue
    print("=" * 80)
    print("aeid", aeid, "|", f.name)
    xl = pd.ExcelFile(f)
    print("sheets:", xl.sheet_names)
    sheet = "General" if "General" in xl.sheet_names else xl.sheet_names[0]
    g = pd.read_excel(f, sheet_name=sheet, header=None)
    with pd.option_context("display.max_colwidth", 200, "display.width", 250):
        print(g.to_string(index=False, header=False))