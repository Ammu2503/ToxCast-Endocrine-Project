import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XLSX = ROOT / "Data" / "Raw" / "assay_annotations_invitrodb_v4_3_AUG2025.xlsx"

ER_TOKENS = {"er", "era", "erb", "her", "ber", "mera", "ere", "eregfp",
             "eraera", "eraerb", "erberb", "eralpha", "erbeta"}


def is_er(name):
    """True for estrogen receptor endpoints. ERR (a different receptor) is excluded."""
    tokens = str(name).lower().split("_")
    if any(t.startswith("err") for t in tokens):
        return False
    if any(t in ER_TOKENS for t in tokens):
        return True
    return "t47d" in tokens and "positive" in tokens  # ACEA cell growth endpoint


def mechanism_guess(name):
    """Rough guess from the endpoint name only. Must be checked against the EPA descriptions."""
    n = str(name).lower()
    if n.startswith("nvs_nr"):
        return "binding (guess)"
    if "eraera" in n or "eraerb" in n or "erberb" in n:
        return "dimerization (guess)"
    if "t47d" in n:
        return "cell growth (guess)"
    if "bla" in n or "luc" in n or n.startswith("atg_") or "eregfp" in n:
        return "gene transcription (guess)"
    return "unknown, check"


def direction(name):
    n = str(name).lower()
    if "antagonist" in n:
        return "antagonist"
    if "agonist" in n:
        return "agonist"
    return "unclear from name"


if __name__ == "__main__":
    df = pd.read_excel(XLSX, sheet_name="annotations_combined")
    name_col = "assay_component_endpoint_name"
    hit = df[df[name_col].apply(is_er)].copy()
    hit["mechanism_guess"] = hit[name_col].apply(mechanism_guess)
    hit["direction_from_name"] = hit[name_col].apply(direction)
    hit["in_current_master"] = hit["aeid"].isin([786, 787, 789, 790, 2053, 2054, 2119, 2120])
    want = ["aeid", name_col, "mechanism_guess", "direction_from_name", "signal_direction",
            "cell_short_name", "detection_technology_type", "assay_function_type",
            "assay_design_type", "biological_process_target", "cell_viability_assay",
            "in_current_master"]
    want = [c for c in want if c in hit.columns]
    hit = hit[want].sort_values(["mechanism_guess", name_col])
    hit.to_csv(ROOT / "results" / "er_assay_survey.csv", index=False)
    print("ER related endpoints found:", len(hit))
    print(hit["mechanism_guess"].value_counts().to_string())
    print(hit[[name_col, "mechanism_guess", "direction_from_name"]].to_string(index=False))