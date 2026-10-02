import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
path = ROOT / "Data" / "Raw" / "assay_annotations_invitrodb_v4_3_AUG2025.xlsx"

df = pd.read_excel(path, sheet_name="annotations_combined")

name = df["assay_component_endpoint_name"].astype(str)
mask = (name.str.upper().str.startswith("TOX21_ER")
        & ~name.str.upper().str.contains("ERR")
        & name.str.upper().str.contains("ANTAGONIST"))
sel = df[mask].copy()

cols = ["aeid", "assay_component_endpoint_name", "signal_direction",
        "cell_short_name", "detection_technology_type",
        "assay_function_type", "cell_viability_assay",
        "burst_assay", "timepoint_hr", "data_usability"]

print(len(sel), "endpoints")
print(sel[cols].to_string())

out = ROOT / "results"
out.mkdir(exist_ok=True)
sel[cols].to_csv(out / "er_antagonist_endpoints.csv", index=False)