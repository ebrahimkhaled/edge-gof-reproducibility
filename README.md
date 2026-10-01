# Reproducibility archive — *A grouped calibration test for logistic regression that tolerates a few corrupted records*

Ebrahim Khaled Ebrahim, Osama Abd El-Aziz Hussein and Ahmed El-Kotory, Department of Applied Statistics, Alexandria University.

This archive holds the pre-declarations, the simulation code, the deposited results and the manuscript
source of the study reported in the paper. It is the archive Section 4 of the paper refers to, and it is
what allows every number in the paper to be checked without re-running anything.

## What is here

```
declarations/      every pre-declaration, with the sha256 that froze it
code/INDEX.md      which script produced each table, figure and block, and which scripts are exploratory
code/simulations/  the study's R code: generators, tests, block runners, self-tests, analyses
code/paper/        the table builder and this archive's assembler
results/analysis/  the declared analysis files the paper's tables are built from
results/blocks/    the per-replicate deposits of the blocks the paper reports, the block summaries of
                   the main battery, and the battery's design files (cell table, sample sizes)
results/cohort/    the outputs of the clinical-cohort analysis of Section 7
manuscript/        the LaTeX source, the generated tables, the figures and the compiled PDF
MANIFEST.sha256    sha256 of every file above
```

## How to run

**Software.** The study ran under R 4.4.1 on Windows. The reported scripts use these packages:

* `data.table` 1.15.4 and `parallel`, everywhere;
* `ebrahim.gof` (from CRAN, version 2.9.0 or later), for the comparators and the identity checks;
* `BAGofT` 1.0.0 and `givitiR` 1.3, for the rival tests;
* `robustbase` 0.99.4, for block 9d;
* `ggplot2` 4.0.3, `scales`, `patchwork` and `ggrepel`, for the figures;
* `RhpcBLASctl`, optionally, to pin BLAS to one thread on the workers.

Install them from CRAN. The scripts were written against the development copy of `ebrahim.gof`. In
this archive, every `devtools::load_all()` of that copy has been replaced by
`library(ebrahim.gof)`, and a comment marks each place.

**The root.** Every script finds its files from the archive root. It takes the root from the first of
these that is available:

1. the environment variable `EDGE_ARCHIVE_ROOT`;
2. when the script is run with `Rscript`, the folder two levels above the script (`code/<folder>/x.R`
   sits two levels below the root);
3. the working directory.

The script then sets `EDGE_ARCHIVE_ROOT` itself, so that the worker processes it starts use the same
root. Setting the variable is the safest choice:

```
cd <archive root>
export EDGE_ARCHIVE_ROOT="$PWD"          # Windows: set EDGE_ARCHIVE_ROOT=%CD%
Rscript code/paper/make_tables_paper3.R
```

**How the paths were made portable.** The study was run in the authors' working folders. When this
archive is assembled, `code/paper/make_archive.py` rewrites every absolute path in the R code by one
fixed table. A short block at the top of each rewritten script defines `edge_path()`, `edge_battery()`
and `edge_out()`.

| in the authors' folders | in this archive |
|---|---|
| `simulations/` (the code) | `code/simulations/` |
| `simulations/battery/analysis/` | `results/analysis/` |
| `simulations/battery/<block or file>` | `results/blocks/<block or file>` |
| `simulations/runI_p3_*.csv` | `results/cohort/` |
| `paper_EDGE/paper3/` | `manuscript/` |
| `paper_EDGE/theory/` | `declarations/` |
| output folders of earlier paper versions | `output/<name>/`, created when a script writes there |
| temporary folders | `tempdir()` |
| the development copy of `ebrahim.gof` | the installed package; scripts that inspect its source read `EBRAHIM_GOF_SRC` |
| the authors' function files from before the package | `code/legacy_not_deposited/` (not deposited; only exploratory scripts use them) |

The assembler refuses to finish while any deposited script still names a folder on the authors'
machine. Nothing else in the code was changed, apart from one comment that misnamed a test (see the
end of this section). Each rewritten line can be compared with the original by searching for
`edge_path`, `edge_battery`, `edge_out` or `archive:`.

**Which script rebuilds what.**

* `Rscript code/paper/make_tables_paper3.R` writes the thirteen generated tables into `manuscript/tables/`.
  Set `EDGE_TABLES_OUT` to write them to another folder and compare them with the deposited ones. The
  output is identical, byte for byte.
* The six figures are written by `_fig13_p3.R` (Figures 1 and 2), `_fig15_p3.R` (Figure 3),
  `_fig12_partition_p3.R` (Figure 4), `_fig10_diabetes_p3.R` (Figure 5) and `_fig11_belt_p3.R`
  (Figure 6), all in `code/simulations/`. Each reads deposited files only.
* `code/INDEX.md` lists, for every table, figure, quoted number and declared block, the script and the
  input files. It also lists the scripts that are not part of the reported study.

A runner writes into `results/blocks/<block>/` by default. It refuses to overwrite a finished block,
so to regenerate one, pass `--out` with a new folder name.

**A correction made in the deposited copy.** The header comments of `analyse_census_extended.R` and
`analyse_corruption_partitions.R` called the HL_w column the "Hosmer–Hjort weighted decile test". The
code behind that column (`bt_hlw()` in `_battery_tests.R`) is the Hosmer–Lemeshow statistic on
equal-width risk groups. That is also what the paper calls it. The comment has been corrected, and no
computation changed.

## The discipline this archive records

Every block of the study was declared before it ran. A declaration fixes the scenarios, the decision
rules, the thresholds and what would count against the claim, and it was hashed before any scenario of
that block was computed; `declarations/` holds one document per block with its hash. Most carry a
single hash. Two carry a short log instead, because the document was revised more than once *before* its
block ran and each revision was recorded with its time and its reason; in those files the live document
matches the **last** hash listed. One document, `PREDECLARATION_restructure_battery_FROZEN_AD.md`, is a
frozen copy kept beside its hashed original and carries no hash of its own.
`code/paper/make_archive.py` re-checks every hash on each build and refuses to assemble the archive if
one has moved.

Declared claims that went against the test are reported as failures in the paper, in the same place and
type as the ones that held: H1, H3 and H4 of the main study, criterion 2 of block 9R2, and criterion 1 of
block 9R3.

## How to check a number in the paper

1. **From the analysis files.** Every table is generated by `code/paper/make_tables_paper3.R` from
   `results/analysis/` and `results/blocks/`. Running it rewrites `manuscript/tables/` from the deposited
   results; the files it writes should be identical to the ones here.
2. **From the per-replicate p-values.** `results/blocks/<block>/` holds one file per scenario with one row
   per replicate and one column per test, so any rejection rate in the paper can be recomputed directly.
3. **From scratch.** Each block has a runner in `code/simulations/` named `run_M_block<block>.R`. Every
   data set is regenerated from a stored seed base plus the replicate index, so a re-run reproduces the
   same data sets. The blocks that were added after the main study also carry an identity gate: they
   recompute the published test on the regenerated data and refuse the result unless it matches the
   stored p-value to 1e-8.

Entry points, by block:

| block | runner | what it reports |
|---|---|---|
| main battery | `run_M_all.R`, which calls `run_M_battery.R` block by block | size and power over the six departure families |
| 8 | `run_M_rivals.R` | the Liu projection test and BAGofT on the link and tail scenarios |
| 8L | `run_M_block8L.R` | le Cessie's smoothed test on the same data sets |
| 8R | `run_M_block8R.R` | the three rivals on the omitted-term and rough-misfit scenarios |
| 9 | `run_M_block9.R` | the contamination study: four corruptions, 84 scenarios |
| 9b, 9bL | `run_M_block9b.R`, `run_M_block9bL.R` | the resampling tests and le Cessie under corruption |
| 9c | `run_M_block9c.R` | the partition surface, data held fixed |
| 9d | `run_M_block9d.R` | maximum likelihood against a robust fit |
| 9R, 9R2, 9R3 | `run_M_block9R.R`, `run_M_block9R2.R`, `run_M_block9R3.R` | the remedies, the hybrid partition, and the tolerance question behind EDGE-FR |
| 9R2n | `run_M_block9R2n.R` | the matched null those variants are size-adjusted against |
| 8S | `run_M_block8S.R` | the three rivals on the link and tail family at n = 200 |
| C1b | `run_M_blockC1b.R` | the sign error: the severe member of the corruption family |
| 9dX | `run_M_block9dX.R` | the robust fit at n = 5000, the paper's own boundary |
| R, R2; R_p3, R2_p3 | `run_M_blockR_realdata.R` (`RR_VAR` picks the covariate, `RR_COHORT=p3` the paper's cohort) | corrupted predictors on the real cohort, outcome drawn from the frozen model |
| RQ | `run_M_blockRQ.R` | robust quasi-deviance versions of Stukel's test under corrupted covariates |
| EXT | `run_M_blockEXT.R` | the external mode against the usual external-validation tests |
| STREAM | `run_M_blockSTREAM.R` | `edge.stream()` under drift and repeated looks |
| T | `bench_time_T.R`, `bench_time_T_lecessie_n3.R` | the one-core timing benchmark |

Two analyses are deposited for reasons that belong in the record rather than in a footnote.
`analyse_block9_B9_4_corrected.R` documents a declared criterion that had to be reformulated after it
was computed: B9.4 was declared as a relative change in every fitted coefficient, and the design's
true intercept is zero, so the relative change in it is not interpretable. The script states the
reformulation, regenerates the corrected file, and reproduces the slope figures the paper quotes.
`analyse_contaminated_power.R` reads power under contamination at the critical value of the
contaminated null rather than the nominal level, which is what the paper's own reading rule requires
and what the first draft of Section 6.5 did not do.

## The clinical data

The application of Section 7 uses **Diabetes 130-US Hospitals for Years 1999–2008** from the UCI Machine
Learning Repository (DOI 10.24432/C5230J). The file is not redistributed here. To reproduce Section 7,
download `diabetic_data.csv` from the repository and place it at `code/data_large/diabetic_data.csv`
(one level above the simulation directory). That is where `bt_real_data0()` in
`code/simulations/_battery_tests.R` and `cohort_p3()` in `code/simulations/_cohort_p3.R` read it. The file the study used is 19,159,383 bytes with

```
sha256  0689e7ec031237dc63031b938805c48377748761a3b26acab621567afa24df97
```

so a reader can confirm they hold the same version. The preparation — the outcome definition, the
predictors and the development/validation split — is in that same function and in `run_I_bigdata.R`,
not done by hand.

## Software

The tests are those of the R package `ebrahim.gof`. The archive's code calls the package for the
comparators and recomputes the directed test itself, operation for operation, in
`code/simulations/_battery_tests.R`; the two agree to 1e-8 by the identity gates named above. The package
version used for the study is stated in the paper.

## What is not here

The per-replicate deposits of the main power battery (blocks 0, 1a, 1b and 2 to 7) come to about
500 MB. They are left out so that the archive stays small enough to download and check. Their block
summaries (`results/blocks/<block>/_summary.csv`) are included, and so is everything the paper reports
from them, in `results/analysis/`. The per-replicate files can be regenerated exactly from the runners
and the stored seeds. The paired McNemar tests of `analyse_M_battery.R` read those per-replicate
files, so that script needs the regenerated blocks before it can be rerun.

## Licence

Code: GPL-3, matching the `ebrahim.gof` package it accompanies. Declarations, results and manuscript
text: CC BY 4.0. Cite the archive by its DOI and the paper by its own.
