import re
import zipfile
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ZIP = Path(r"C:\Users\harik\Downloads\pubchem_invitrodb_v4_3_19AUG2026.zip")

s = pd.read_csv(ROOT / "results" / "er_assay_survey.csv")
name = "assay_component_endpoint_name"

# Rev-ErbA is a different receptor (not estrogen receptor beta). The raw channel
# readouts (ch1, ch2) are not separate assays.
s = s[~s[name].str.contains("REV_ERB", case=False)]
s = s[~s[name].str.contains(r"_ch[12]$", case=False)]

in_zip = {}
with zipfile.ZipFile(ZIP) as z:
    for info in z.infolist():
        m = re.search(r"aeid(\d+)_", info.filename)
        if m:
            in_zip[int(m.group(1))] = info.file_size

s["file_in_zip"] = s["aeid"].isin(in_zip.keys())
s["file_size_kb"] = s["aeid"].map(lambda a: round(in_zip.get(a, 0) / 1024))
out = s[["aeid", name, "mechanism_guess", "direction_from_name", "signal_direction",
         "cell_short_name", "file_in_zip", "file_size_kb", "in_current_master"]]
out = out.sort_values(["mechanism_guess", name])
out.to_csv(ROOT / "results" / "er_assay_availability.csv", index=False)
print("Endpoints after removing Rev-ErbA and ch1/ch2:", len(out))
print("Found in the zip:", int(out["file_in_zip"].sum()))
print(out.to_string(index=False))