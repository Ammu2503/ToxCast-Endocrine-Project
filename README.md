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
- DTXSID: the EPA chemical ID used to join all tables.

## Data sources
- EPA assay annotations: assay_annotations_invitrodb_v4_3_AUG2025.xlsx, from the EPA ToxCast download page (invitrodb v4.3, released August 2025). Downloaded 2026-10-01. Stored in Data/Raw in this repo (small file) so scripts 02 and 03 can be rerun.
- Assay endpoint list: assay-lists-2026-10-01.csv, from the CompTox Dashboard assay endpoints page. Downloaded 2026-10-01. Early exploration only. Kept local.
- Chemical level results: invitrodb v4.3 summary files, one Excel file per endpoint (aeid), from pubchem_invitrodb_v4_3_19AUG2026.zip on the EPA Clowder site. Downloaded 2026-10-06. Only 8 files are used (aeid 786, 787, 789, 790, 2053, 2054, 2119, 2120). Kept local.
- Chemical structures: CompTox Dashboard Batch Search export (CCD-Batch-Search_2026-10-06_11_48_06.xlsx) for the 8,947 DTXSIDs in the master table. Kept local.

Everything in the Data folder except the annotations workbook stays on the local computer (Data/ is listed in .gitignore). To repeat the work, download the files above and place them in Data/Raw.

## Folder layout
- Data/Raw: original downloads, never edited by hand
- Data/Processed: local working files (not in Git)
- Scripts: one numbered script per step
- results: small output tables that are safe to share
- docs: notes and documentation

## Scripts
Run in this order. Python packages used: pandas, openpyxl, rdkit.
- 02_er_antagonist_endpoints.py: selects endpoints whose names start with TOX21_ER, contain ANTAGONIST and do not contain ERR (ERR is a different receptor) from the annotations workbook. Output: results/er_antagonist_endpoints.csv (12 endpoints).
- 03_endpoint_decisions.py: adds an include column (include, pending or exclude) and a reason for each endpoint. Output: results/er_antagonist_endpoint_decisions.csv.
- 05_extract_endpoint_files.py: pulls the 8 endpoint files out of the EPA zip into Data/Raw/aeid_files.
- 06_inspect_aeid_files.py and 07_check_results.py: inspection of the file layout and outcome codes.
- 08_build_long_table.py: stacks the 8 files into one long table (Data/Processed/er_antagonist_long.csv).
- 09_build_master.py: builds the master table and the missingness summary. Output: results/er_antagonist_missingness.csv.
- 10_cytotox_and_agreement.py: overlap of antagonist hits with viability hits, and agreement between endpoints. Output: results/er_antagonist_cytotox_check.csv.
- 11_cytotox_rule.py: AC50 based cytotoxicity rule. Output: results/er_antagonist_cytotox_rule.csv.
- 12_cytotox_sensitivity.py: repeats the rule at several cutoffs. Output: results/er_antagonist_cytotox_sensitivity.csv.
- 13a_inspect_structures.py and 13_morgan_fingerprints.py: structure export check, then Morgan fingerprints. Output: results/structure_summary.csv and results/smiles_failed.csv.

## Endpoint decisions
Included (3 primary ERa antagonist endpoints and their viability partners):
- 786 TOX21_ERa_BLA_Antagonist_ratio (HEK293T cells, 24 hr, fluorescence)
- 787 TOX21_ERa_BLA_Antagonist_viability (viability partner of 786)
- 789 TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2 (VM7 cells, 22 hr, luminescence)
- 790 TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2_viability (viability partner of 789)
- 2053 TOX21_ERa_LUC_VM7_Antagonist_0.1nM_E2 (same assay as 789 with less estradiol)
- 2054 TOX21_ERa_LUC_VM7_Antagonist_0.1nM_E2_viability (viability partner of 2053)

Pending (ER beta, waiting on scope confirmation). Both are in the master table:
- 2119 TOX21_ERb_BLA_Antagonist_ratio
- 2120 TOX21_ERb_BLA_Antagonist_viability

Excluded from the master table (raw BLA channel readouts that feed the ratio, not separate assays):
- 1189, 1190 (ERa ch1 and ch2)
- 2117, 2118 (ERb ch1 and ch2)

Also left out: ERR assays (estrogen-related receptor, a different target), agonist assays, binding only assays and steroidogenesis assays.

## Rules used
- Outcome codes in the EPA files: 2 = active (HITC 0.90 or higher), 1 = inactive, 3 = wrong direction (negative HITC).
- Master table values: 1 = active, 0 = not active (codes 1 and 3 together), blank = chemical never tested in that endpoint.
- If a chemical has several samples for one endpoint, it is called active if any sample is active. Samples disagree for 3.5% to 8.3% of the chemicals that have repeats, so a stricter rule still needs to be tested.

## Results so far
Master table (Data/Processed/er_antagonist_master.csv, local):
- 8,947 unique chemicals and 8 endpoints. 9.6% of cells are missing overall. 7,229 chemicals have results in all 8 endpoints.
- Missing by endpoint: 7.2% for 786, 787, 789, 790 (8,305 chemicals tested) and 12.0% for 2053, 2054, 2119, 2120 (7,871 tested).
- A chemical was tested in all four endpoints of a group or in none, so the missing data comes from the screening design and is not random.
- Active chemicals (share of those tested): 786 13.1%, 789 10.2%, 2053 12.7%, 2119 18.7%. Viability partners: 787 3.9%, 790 7.3%, 2054 7.3%, 2120 13.6%.

Cytotoxicity:
- Share of antagonist hits that are also active in the paired viability endpoint: 26.0% (786), 47.6% (789), 39.5% (2053), 57.6% (2119).
- Rule tested: keep an antagonist hit if the viability endpoint is inactive, or if the antagonist AC50 is at least 3 times lower than the viability AC50. Hits kept: 82.2% (786), 64.2% (789), 68.8% (2053), 51.8% (2119).
- The 3 fold cutoff is a working choice, not an EPA rule. It is still to be checked against the EPA method and discussed.

Agreement between primary endpoints (correlation of active calls): 786 vs 789 0.53, 786 vs 2053 0.55, 789 vs 2053 0.76, 786 vs 2119 0.66. 419 chemicals are active in all four primary endpoints and 2,036 in at least one. The two LUC conditions are the same assay at two estradiol levels, so their agreement is expected.

Structures:
- Morgan fingerprints, radius 2, 2048 bits, from QSAR ready SMILES with plain SMILES as a fallback.
- 8,619 of 8,947 chemicals have a fingerprint (8,212 from QSAR ready SMILES, 407 from plain SMILES). 328 chemicals (3.7%) have no usable structure, listed in results/smiles_failed.csv.

## Notes and limits
- The 0.1 nM and 0.5 nM E2 endpoints are two conditions of the same assay, not two independent experiments.
- Endpoint reasons come from the annotation workbook columns (signal_direction, assay_function_type, cell_viability_assay). They still need checking against the EPA assay description documents. For 2053 the workbook says 22 hr but the EPA description says 48 hours.
- Viability activity is a flag, not proof, that an antagonist hit is caused by cell death.
- Fingerprint settings (radius 2, 2048 bits) are defaults and have not been tested on this data.
- Why the 328 structures failed has not been checked yet.

## Open issues
- Scope question: whether estrogen antagonist includes ER beta.
- Which cytotoxicity rule and cutoff to use.
- Whether treating inactive and wrong direction results together as not active is acceptable.
- Check each endpoint against the EPA assay description documents.

## Next steps
- Check whether chemicals without a fingerprint differ in activity or missing status.
- Study the missing data pattern against chemical structure.
- Imputation test: hide known results, fill them in (for example with MICE), and compare against a simple baseline.