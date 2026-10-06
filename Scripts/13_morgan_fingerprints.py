import pandas as pd
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
f = ROOT / "Data" / "Raw" / "CCD-Batch-Search_2026-10-06_11_48_06.xlsx"

d = pd.read_excel(f, sheet_name="Main Data")
print("rows:", len(d), "| unique DTXSID:", d["DTXSID"].nunique())
print("missing SMILES:", int(d["SMILES"].isna().sum()),
      "| missing QSAR-ready SMILES:", int(d["QSAR_READY_SMILES"].isna().sum()))

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rows, bad, used = [], [], {"qsar_ready": 0, "plain": 0}
for _, r in d.iterrows():
    mol, src = None, None
    for col, name in [("QSAR_READY_SMILES", "qsar_ready"), ("SMILES", "plain")]:
        smi = r[col]
        if pd.notna(smi):
            mol = Chem.MolFromSmiles(str(smi))
            if mol is not None:
                src = name
                break
    if mol is None:
        bad.append({"DTXSID": r["DTXSID"], "PREFERRED_NAME": r["PREFERRED_NAME"],
                    "SMILES": r["SMILES"], "QSAR_READY_SMILES": r["QSAR_READY_SMILES"]})
        continue
    used[src] += 1
    rows.append([r["DTXSID"]] + [int(x) for x in gen.GetFingerprintAsNumPy(mol)])

fp = pd.DataFrame(rows, columns=["dtxsid"] + [f"b{i}" for i in range(2048)])
fp.to_csv(ROOT / "Data" / "Processed" / "morgan_fingerprints.csv", index=False)
pd.DataFrame(bad).to_csv(ROOT / "results" / "smiles_failed.csv", index=False)

master = pd.read_csv(ROOT / "Data" / "Processed" / "er_antagonist_master.csv", index_col=0)
with_fp = master.index.isin(fp["dtxsid"])

summary = pd.DataFrame([{
    "chemicals_in_export": len(d),
    "fingerprints_made": len(fp),
    "failed_or_no_smiles": len(bad),
    "used_qsar_ready": used["qsar_ready"],
    "used_plain_smiles": used["plain"],
    "master_chemicals": len(master),
    "master_chemicals_with_fingerprint": int(with_fp.sum()),
}])
print(summary.T.to_string(header=False))
summary.to_csv(ROOT / "results" / "structure_summary.csv", index=False)