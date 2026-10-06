# ToxCast Endocrine Project

## Goal
Build a master dataset for estrogen antagonist assays from EPA ToxCast/Tox21 data, describe how much assay data is missing, and add chemical structure information. This prepares the ground for later missing-data work (for example imputation with MICE) and comparison with an EPA consensus approach.

Scope for now: estrogen antagonist assays only. Agonist, androgen and other pathways are left out until this one is clean.

## Key terms
- Assay: one lab test run on cells, for example TOX21_ERa_BLA_Antagonist_ratio.
- Endpoint: one analysis of an assay readout. Each endpoint has an ID called aeid.
- Antagonist: a chemical that blocks the receptor, so the signal goes down (signal_direction = loss).
- Viability endpoint: a cell health test on the same cells. It helps show whether a signal drop is only dead cells.
- BLA and LUC: two ways of measuring the signal (fluorescence from beta-lactamase, or light from luciferase).

## Data sources 
- EPA assay annotations: assay_annotations_invitrodb_v4_3_AUG2025.xlsx, from the EPA ToxCast download page (invitrodb v4.3, released August 2025). Downloaded 2026-10-01. Place it in Data/Raw.
- Assay endpoint list: assay-lists-2026-10-01.csv, from the CompTox Dashboard assay endpoints page. Downloaded 2026-10-01. Place it in Data/Raw. Used for early exploration only.

The small annotations workbook is stored in Data/Raw so the scripts can be rerun. Other raw and processed data stay on the local computer (the Data folder is listed in .gitignore, and the workbook was added with git add -f).d data stay on the local computer.

## Folder layout
- Data/Raw: original downloads, never edited by hand
- Data/Processed: local working files (not in Git)
- Scripts: one numbered script per step
- results: small output tables that are safe to share
- docs: notes and documentation

## Scripts
- Scripts/02_er_antagonist_endpoints.py
  Reads the annotations_combined sheet of the EPA workbook. Selects endpoints whose names start with TOX21_ER, contain ANTAGONIST and do not contain ERR (ERR is a different receptor). Output: results/er_antagonist_endpoints.csv (12 endpoints).
- Scripts/03_endpoint_decisions.py
  Reads the file above and adds an include column (include, pending or exclude) and a reason for each endpoint. Output: results/er_antagonist_endpoint_decisions.csv.

To run: put the workbook in Data/Raw, then run script 02, then script 03. Python packages used: pandas, openpyxl.

## Endpoint decisions so far
Included (3 primary ERa antagonist endpoints and their viability partners):
- 786 TOX21_ERa_BLA_Antagonist_ratio (HEK293T cells, 24 hr, fluorescence)
- 787 TOX21_ERa_BLA_Antagonist_viability (viability partner of 786)
- 789 TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2 (VM7 cells, 22 hr, luminescence)
- 790 TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2_viability (viability partner of 789)
- 2053 TOX21_ERa_LUC_VM7_Antagonist_0.1nM_E2 (same assay as 789 with less estradiol)
- 2054 TOX21_ERa_LUC_VM7_Antagonist_0.1nM_E2_viability (viability partner of 2053)

Pending (ER beta, waiting on scope confirmation):
- 2119 TOX21_ERb_BLA_Antagonist_ratio
- 2120 TOX21_ERb_BLA_Antagonist_viability

Excluded from the master sheet (raw BLA channel readouts that feed the ratio, not separate assays):
- 1189, 1190 (ERa ch1 and ch2)
- 2117, 2118 (ERb ch1 and ch2)

Also left out: ERR assays (estrogen-related receptor, a different target), agonist assays, binding-only assays and steroidogenesis assays.

## Notes and limits
- The 0.1 nM and 0.5 nM E2 endpoints are the same assay run at two estradiol levels. They are two conditions, not two independent experiments.
- Endpoint reasons so far come from the annotation workbook columns (signal_direction, assay_function_type, cell_viability_assay). They still need checking against the EPA assay description documents.
- Many apparent antagonist hits can be caused by cell death, which is why each primary endpoint is paired with a viability endpoint. A rule for handling cytotoxic hits is still to be decided and documented.

## Open issues
- Chemical-level results (active or inactive call per chemical per endpoint) have not been obtained yet. The CompTox Dashboard Batch Search export is being tested to see whether it separates tested-and-inactive from never-tested chemicals. If it cannot, other routes will be checked (CTX Bioactivity API, invitrodb database package).
- Scope question: whether estrogen antagonist includes ER beta.

## Next steps
- Get chemical-level hit calls for endpoints 786, 787, 789, 790, 2053 and 2054.
- Build the master spreadsheet with one row per chemical and one column per endpoint. Document the rule for duplicate results.
- Describe the data: number of unique chemicals, number of endpoints, percent missing overall and per endpoint, active and inactive counts.
- Add SMILES for each chemical and generate Morgan fingerprints with RDKit. Report how many structures fail.