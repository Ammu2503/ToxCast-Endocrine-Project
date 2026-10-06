import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
df = pd.read_csv(ROOT / "results" / "er_antagonist_endpoints.csv")

decisions = {
    786:  ("include", "Primary ERa antagonist, BLA, loss of signal"),
    789:  ("include", "Primary ERa antagonist, LUC VM7, 0.5 nM E2"),
    2053: ("include", "Primary ERa antagonist, LUC VM7, 0.1 nM E2"),
    787:  ("include", "Viability partner of 786"),
    790:  ("include", "Viability partner of 789"),
    2054: ("include", "Viability partner of 2053"),
    2119: ("pending", "ER beta, depends on whether beta is in scope"),
    2120: ("pending", "Viability partner of 2119"),
    1189: ("exclude", "Raw BLA channel, feeds the ratio, not a separate assay"),
    1190: ("exclude", "Raw BLA channel, feeds the ratio, not a separate assay"),
    2117: ("exclude", "Raw BLA channel, feeds the ratio, not a separate assay"),
    2118: ("exclude", "Raw BLA channel, feeds the ratio, not a separate assay"),
}

df["include"] = df["aeid"].map(lambda a: decisions[a][0])
df["reason"] = df["aeid"].map(lambda a: decisions[a][1])
df.to_csv(ROOT / "results" / "er_antagonist_endpoint_decisions.csv", index=False)
print(df[["aeid", "assay_component_endpoint_name", "include", "reason"]].to_string())