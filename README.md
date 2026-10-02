# ToxCast Endocrine Project

Goal: build a master dataset for estrogen antagonist assays (ToxCast/Tox21) and study missing assay data.

## Data sources (not stored in this repo)
- EPA assay annotations: assay_annotations_invitrodb_v4_3_AUG2025.xlsx, from the EPA ToxCast download page (invitrodb v4.3), downloaded 2026-10-01
- Assay endpoint list: assay-lists-2026-10-01.csv, from the CompTox Dashboard, downloaded 2026-10-01

## Steps done
- Scripts/02_er_antagonist_endpoints.py selects the ER antagonist endpoints, their viability partners and the raw channel readouts from the annotations workbook.
- Output: results/er_antagonist_endpoints.csv (12 endpoints)

## Endpoint decisions so far
- Include: 786, 789, 2053 (primary ERa antagonist endpoints) and 787, 790, 2054 (their viability partners)
- Pending: 2119, 2120 (ER beta), until scope is confirmed
- Exclude from master sheet: 1189, 1190, 2117, 2118 (raw BLA channel readouts that feed the ratio)

## Next
- Get chemical-level hit calls, build the master spreadsheet, describe missing data, add structures (SMILES, RDKit, Morgan).