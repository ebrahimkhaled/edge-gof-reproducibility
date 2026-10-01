# Index of the code: what produced each reported result

This index lists, for every table, figure and declared block of the paper, the script that produces it
and the files it reads. Paths are relative to the archive root. Every script resolves its paths from
the archive root as described in `README.md`, "How to run". The study and its analysis ran in a working
folder called `battery/`. In the archive that folder is split in two: `results/analysis/` holds what was
in `battery/analysis/`, and `results/blocks/` holds everything else.

The scripts in `code/simulations/` are of four kinds:

* **Reported.** These scripts produce something the paper reports. They are listed in sections 1–5.
* **Shared code.** These are libraries that the reported scripts source (section 6).
* **Checks made before a block ran.** These are self-tests, independent reviews and verifications of
  the code that the reported scripts use (section 7). They support the results but do not produce any.
* **Not part of the reported study.** This covers exploratory work, earlier papers, package development
  and audits (section 8). These scripts are deposited so that the record is complete. The paper does
  not rely on them.

## 1. Tables

All thirteen generated tables are written by **`code/paper/make_tables_paper3.R`**. It reads deposited
analysis files only and recomputes nothing from raw p-values.

| table | file in `manuscript/tables/` | what the builder reads | script that wrote those files |
|---|---|---|---|
| 1 | `tab_size.tex` | `results/analysis/rule_A_gate.csv` | `analyse_M_battery.R` |
| 2 | `tab_families.tex` | `results/analysis/rule_A_families.csv` | `analyse_M_battery.R` |
| 3 | `tab_rivals.tex` | `results/blocks/8/_paired.csv`, `results/blocks/8L/_paired.csv`, `results/analysis/rule_A_membership.csv`, `results/blocks/T/timing.csv`, `results/blocks/T/timing_lecessie_n3.csv` | `run_M_rivals.R --summary` (8), `run_M_block8L.R` (8L), `analyse_M_battery.R`, `bench_time_T.R`, `bench_time_T_lecessie_n3.R` |
| 4 | `tab_severity.tex` | `results/analysis/corruption_partitions.csv` | `analyse_corruption_partitions.R` (reads blocks 9 and C1b) |
| S1 | `tab_hypotheses.tex` | `results/analysis/hypotheses.csv`, `results/analysis/mcnemar.csv` | `analyse_M_battery.R` |
| S2 | `tab_rivals_contam.tex` | `results/blocks/9b/analysis/rates_9b.csv`, `seconds_9b.csv`, `results/blocks/9bL/_summary.csv` | `analyse_block9b.R`, `run_M_block9bL.R` |
| S3 | `tab_robust.tex` | `results/blocks/9d/analysis_9d.csv` | `analyse_block9d.R` |
| S4 | `tab_edgefr.tex` | `results/analysis/edgefr.csv`, `results/analysis/edgefr_paired.csv` | `analyse_edgefr.R` (reads blocks 9R, 9R2, 9R2n, 9R3) |
| S5 | `tab_contamination.tex` | `results/blocks/9/analysis/cell_test.csv`, `results/blocks/8L/_summary.csv`, `results/blocks/T/*.csv` | `analyse_block9.R`, `run_M_block8L.R`, `bench_time_T*.R` |
| S6 | `tab_external.tex` | `results/blocks/EXT/_summary.csv` | `analyse_blockEXT.R` |
| S7 | `tab_pairedfam.tex` | `results/analysis/paired_families.csv` | `analyse_paired_families.R` |
| S8 | `tab_tw.tex` | `results/blocks/TW/_summary.csv` | `run_M_blockTW.R`, `run_M_blockTW_HH.R`, `analyse_blockTW.R` |
| S9 | `tab_twext.tex` | `results/blocks/TW_ext/_summary.csv` | `run_M_blockTW_ext.R` |

Two tables are typed in the LaTeX source rather than generated:

* **Table 5 (`tab:cohortcorrupt`, Section 7).** This gives the rejection rates at the 5% level from the
  per-replicate files in `results/blocks/R2_p3/`, which were produced by
  `run_M_blockR_realdata.R` with `RR_COHORT=p3 RR_VAR=number_inpatient`. Each entry is the proportion of
  replicates with p ≤ 0.05. No separate analysis script exists for this table.
* **Table 6 (`tab:whennot`, Section 8).** This table contains no numbers.

## 2. Figures

All the figures are written to `manuscript/figures/`. Each figure script sources `_ek_theme.R`.

| figure | file | script | reads |
|---|---|---|---|
| 1 | `Fig13_scenarios.pdf` | `_fig13_p3.R` | `results/analysis/census_extended.csv`, `results/analysis/_census_power_paper2.csv`, `results/analysis/rule_A_membership.csv`, `results/blocks/{8,8L,8R,8S}/_summary.csv` |
| 2 | `Fig14_census.pdf` | `_fig13_p3.R` (the same run) | as Figure 1 |
| 3 | `Fig15_cost.pdf` | `_fig15_p3.R` | `results/blocks/T/*.csv`, `results/blocks/{8,8L}/_summary.csv`, `results/analysis/census_extended.csv`, `results/analysis/corruption_partitions.csv`, `results/blocks/9/analysis/cell_test.csv` |
| 4 | `Fig12_partition.pdf` | `_fig12_partition_p3.R` | `results/blocks/9c/analysis/cell_G_test.csv` (from `analyse_block9c.R`) |
| 5 | `Fig10_diabetes.pdf` | `_fig10_diabetes_p3.R` | `results/cohort/runI_p3_{calcurve,dev,summary}.csv` |
| 6 | `Fig11_belt.pdf` | `_fig11_belt_p3.R` | `results/cohort/runI_p3_calcurve.csv` |

`manuscript/figures/Fig16_rough_omitted.pdf` was drawn by `_fig16_8R.R` for the paper-2 version. The
submitted paper does not include it.

The file `results/analysis/_census_power_paper2.csv` is the census cache that `_fig13_scenarios.R`
writes from the block summaries. `analyse_census_extended.R` recomputes the tests in that file and
stops if any of them disagrees with it.

## 3. Numbers quoted in the text

| where | number | script | output |
|---|---|---|---|
| Section 3 | θ ≈ 0.3 (x → 4x) and ≈ 0.7 (x → −4x) | `check_theta.R` | `results/analysis/theta.csv` |
| Section 3 and Appendix (bounded-influence proposition) | the numerical check of the bound | `robust_influence_check.R`, `robust_influence_check2.R` | written to `declarations/` (the original wrote them to the theory folder) |
| Section 5 (census, the HL_w, PH, Tsiatis and Xie columns) | lead-or-tie counts | `analyse_census_extended.R` | `results/analysis/census_extended*.csv` |
| Section 5 (cost frontier, asymmetric links) | cubic LR against EDGE | `_cubic_frontier.R` | console |
| Section 6 (sign error) | slope moves by 14% | `check_slope_signerror.R` | `results/analysis/slope_signerror.csv` |
| Section 6 (EDGE-FR) | quoted values | `_edgefr_numbers.R` | console |
| Section 6 (the cost in power) | power at the contaminated null | `analyse_contaminated_power.R` | `results/analysis/contaminated_power.csv` |
| Section 6 / SI S4 | the robust fit at n = 5000 (block 9dX) | `run_M_block9dX.R` | `results/blocks/9dX/` (rates are read from the per-replicate files) |
| Section 7 | the cohort, the partition and validation | `run_I_bigdata_p3.R` (with `_cohort_p3.R`, `_pstar_giviti_harness.R`) | `results/cohort/runI_p3_{dev,val,calcurve}.csv` |
| Section 7 | calibration slope, CITL, AUC, Cox test | `run_I_p3_calslope.R` | `results/cohort/runI_p3_summary.csv` |
| Section 7 | the external-mode statistic | `edge_external.R` (when run as a script it checks itself against `runI_p3_val.csv`) | console |
| Section 8 | the basis at each usage boundary | `analyse_basis_boundaries.R` | `results/analysis/basis_at_boundaries.csv` |
| SI S0, criterion B9.4 | the reformulated criterion | `analyse_block9_B9_4_corrected.R` | `results/blocks/9/analysis/B9_4_fit_corrected.csv` |
| SI S8 | streaming monitor: level, repeated looks, power, cost | `run_M_blockSTREAM.R` | `results/blocks/STREAM/_summary.csv` |

## 4. The declared blocks

Each block was declared and hashed before it ran. The declaration is in `declarations/`. The runner
regenerates every data set from its stored seed. Where an analysis script exists, it reads the
runner's per-replicate output.

| block | declaration | runner | analysis | deposit |
|---|---|---|---|---|
| main battery (0, 1a, 1b, 2–7) | `PREDECLARATION_restructure_battery.md` | `run_M_all.R` → `run_M_battery.R` (cells: `_battery_cells.R`; sample sizes: `battery_nplan.R`) | `analyse_M_battery.R` | summaries in `results/blocks/<0–7>/`; analysis in `results/analysis/`; per-replicate files not deposited |
| 8 | same, section E9 | `run_M_rivals.R` | `run_M_rivals.R --summary` | `results/blocks/8/` |
| 8L | `PREDECLARATION_block8L_lecessie.md` | `run_M_block8L.R` | `analyse_block8L.R` | `results/blocks/8L/` |
| 8R | `PREDECLARATION_block8R_rivals_rough_omitted.md` | `run_M_block8R.R` | its own `--summary` | `results/blocks/8R/` |
| 8S | `PREDECLARATION_block8S_smalln_rivals.md` | `run_M_block8S.R` | its own `--summary` | `results/blocks/8S/` |
| 9 | `PREDECLARATION_block9_contamination.md` | `run_M_block9.R` (`_block9_contam.R`) | `analyse_block9.R`, `analyse_block9_B9_4_corrected.R` | `results/blocks/9/` |
| 9b | `PREDECLARATION_block9b_slow_rivals_robustness.md` | `run_M_block9b.R` (`_block9b_rivals.R`) | `analyse_block9b.R` | `results/blocks/9b/` |
| 9bL | `PREDECLARATION_block9bL_lecessie.md` | `run_M_block9bL.R` | in the runner | `results/blocks/9bL/` |
| 9c | `PREDECLARATION_block9c_groupsize.md` | `run_M_block9c.R` | `analyse_block9c.R` | `results/blocks/9c/` |
| 9d | `PREDECLARATION_block9d_robust.md` | `run_M_block9d.R` (`_block9d_robust.R`) | `analyse_block9d.R` | `results/blocks/9d/` |
| 9R | `PREDECLARATION_block9R_remedies.md` | `run_M_block9R.R` (`_block9R_remedies.R`) | `analyse_edgefr.R` | `results/blocks/9R/` |
| 9R2, 9R2n | `PREDECLARATION_block9R2_hybrid.md` | `run_M_block9R2.R` (`_block9R2_hybrid.R`), `run_M_block9R2n.R` | `analyse_block9R2.R`, `analyse_edgefr.R` | `results/blocks/9R2/`, `results/blocks/9R2n/` |
| 9R3 | `PREDECLARATION_block9R3_tolerance.md` | `run_M_block9R3.R` | `analyse_edgefr.R` | `results/blocks/9R3/` |
| C1b, 9dX | `PREDECLARATION_blockZ_signerror_and_robust5000.md` | `run_M_blockC1b.R` (`_blockC1b_contam.R`), `run_M_block9dX.R` | `analyse_blockC1b.R`, `analyse_corruption_partitions.R` | `results/blocks/C1b/`, `results/blocks/9dX/` |
| R, R2 | `PREDECLARATION_blockR_realdata_corruption.md`, `PREDECLARATION_blockR2_realdata_leverage.md` | `run_M_blockR_realdata.R` (`RR_VAR` selects the covariate) | rates read from the per-replicate files | `results/blocks/R/`, `results/blocks/R2/` (paper-2 cohort, superseded by the next row) |
| R_p3, R2_p3 | as R and R2, rerun on the paper-3 cohort | `run_M_blockR_realdata.R` with `RR_COHORT=p3` | rates read from the per-replicate files (Table 5) | `results/blocks/R_p3/`, `results/blocks/R2_p3/` |
| RQ | `PREDECLARATION_blockRQ_robust_tests.md` | `run_M_blockRQ.R` | `analyse_blockRQ.R` | `results/blocks/RQ/` |
| EXT | `PREDECLARATION_blockEXT_external_validation.md` | `run_M_blockEXT.R` (`edge_external.R`) | `analyse_blockEXT.R`, `analyse_blockEXT_paired.R` | `results/blocks/EXT/` |
| STREAM | `PREDECLARATION_blockSTREAM.md` | `run_M_blockSTREAM.R` | in the runner | `results/blocks/STREAM/` |
| T | `PREDECLARATION_blockT_computation.md` | `bench_time_T.R`, `bench_time_T_lecessie_n3.R` | in the runner | `results/blocks/T/` |

The design files that the runners read are at the top of `results/blocks/`. These are `cells.csv`
(the cell table, with every seed base), `nplan.csv` (the sample sizes), `launch_plan.csv`,
`RULE_weighting.txt`, `seed_overlaps.csv` and `old_id_check.csv`.

The runners refuse to overwrite a finished block. To regenerate a block, pass `--out` with a new folder
(for example `--out battery/9_rerun`). A relative `--out` is created under `code/simulations/`.

## 5. Real data

The cohort analysis reads `code/data_large/diabetic_data.csv`. This file is not redistributed. See
`README.md`, "The clinical data". The UIS data used by the real-data cells of the main battery are in
`code/simulations/uis_data.rds`.

## 6. Shared code (sourced by the reported scripts)

| file | role |
|---|---|
| `_battery_tests.R` | one replicate of every test, the harness's own implementations (`bt_*`), and the real-data loaders |
| `_battery_cells.R` | the cell table and the seed rules |
| `_dgp_library.R` | data-generating processes |
| `_harness.R` | parallel scaffolding of the earlier grids, also used by some runners |
| `_proj_test.R` | the Liu et al. projection test |
| `_block9_contam.R`, `_block9b_rivals.R`, `_block9d_robust.R`, `_block9R_remedies.R`, `_block9R2_hybrid.R`, `_blockC1b_contam.R` | per-block replicate functions |
| `_cohort_p3.R`, `_pstar_giviti_harness.R`, `edge_external.R` | the Section 7 cohort, GiViTI and the external mode |
| `_ek_theme.R` | the figure style |

## 7. Checks made before a block ran (not result producers)

These scripts are `battery_selftest.R`, `battery_dryrun.R`, `battery_review_*.R`, `battery_verify2_*.R`,
`battery_verify_sym1_raw.R`, `rivals_selftest.R`, `block8L_selftest.R`, `block8R_selftest.R`,
`block9_selftest.R`, `block9d_selftest.R`, `block9R2_selftest.R`, `analyse_selftest.R`,
`analysis_review_*.R` and `analysis_verify2_*.R`.

Several of these scripts build or read scratch roots that are not deposited:
`battery/dryrun`, `battery/_review`, `battery/_test`. They also test the path guards of the runners
against the author's folder layout. They document how the code was checked. They are not expected to
run unchanged from the archive.

## 8. Not part of the reported study

The following scripts are deposited for completeness. The paper does not rely on them.

* **Earlier papers and early exploration.** These are `Paper_Simu_Inshallah.R`, `_addvar.R`,
  `_bench_*.R`, `_combo.R`, `_comprehensive.R`, `_enh_*.R`, `_ensemble.R`, `_fig_sims*.R`,
  `_figures_plot*.R`, `_final_*.R`, `_funcs_only.R`, `_harness_selfcheck.R`, `_imhof_check.R`,
  `_master*.R`, `_moredef*.R`, `_null_validation.R`, `_realdata_*.R`, `_separation*.R`, `_stackvote.R`,
  `_worked_example.R` and `wp_b*.R`. Some of them source the author's function files from before the
  package existed. Those files are not deposited, and their path maps to `code/legacy_not_deposited/`.
* **Figures and tables of earlier versions of this paper.** These are `_fig1_*`–`_fig9_*`,
  `_fig10_diabetes.R`, `_fig11_belt.R`, `_fig12_partition.R`, `_fig13_final.R`, `_fig13_scenarios.R`,
  `_fig13_versions.R`, `_fig15_cost.R`, `_fig15_options.R`, `_fig16_8R.R` and `emit_tables_*.R`. They
  write to `output/`.
* **The earlier experiment grids and probes.** These are `grid_*.R`, `bench_compute.R`,
  `bench_slow_timing.R`, `null_calibration_checks.R`, `run_all.R`, `run_A_*` to `run_L*`,
  `run_I_bigdata.R` (the paper-2 cohort), `scout_envelope.R`, `probe_probit.R`, `clinical_visibility.R`,
  `giviti_mechanism_check.R`, `grouped_*_check.R`, `grouped_weight_poly3.R`, `analyse_A.R`,
  `analyse_paired.R`, `analyse_runL.R`, `uis_*.R`, `recount*.R`, `make_headline_recount.R` and
  `proof_check_all.R`.
* **Audits of the earlier runs.** These are `audit_*.R`.
* **Package development.** These are `pkg280_*.R` and `map_pkg_*.R`, plus `map_manuscript_*.R` and
  `battery_verify2_srcident.R`. They check `ebrahim.gof` against its own source tree. They need a
  source checkout of the package, which is passed as `EBRAHIM_GOF_SRC`. Their outputs map to
  `declarations/` or to `tempdir()`.
* **Numbers quoted in the text from single analyses.** `analyse_PR_reported.R` (the Pulkstenis--Robinson
  paragraph of S10, output `results/analysis/PR_*.csv`), `run_I_p3_ici.R` (the calibration index of
  Section 7, output `results/cohort/runI_p3_ici.csv`) and `analyse_corruption_ici.R` (the calibration
  index of the clean records under corruption, Section 6.1, output `results/analysis/corruption_ici.csv`).
* **Robustness exploration.** This is `robust_influence_rivals.R`, which its header describes as
  exploratory.
