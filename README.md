# Balanced subsampling + subject-level permutation K-S analysis

`balanced_subsample_perm_ks_analysis_v2.m` compares the **distribution** of a measurement between experimental groups. It is designed for data in which each group has several subjects (for example, animals) and each subject contributes many observations (for example, trials).

It uses the two-sample Kolmogorov-Smirnov (K-S) statistic and adds two safeguards that a plain K-S test on pooled data lacks.

1. **Every subject gets equal weight.** Pooling all observations would let subjects with many trials dominate their group. Instead, the script draws the same number of observations from every subject before comparing groups.
2. **Subjects, not observations, are the unit of inference.** A K-S test on thousands of pooled trials treats every trial as independent and gives misleadingly small p-values. Here, significance comes from shuffling whole **subjects** between groups, so the p-value reflects variability between subjects.

When more than two groups are compared, the script reports p-values with three multiple-comparison corrections.

## Contents

- [Requirements](#requirements)
- [Preparing the input Excel file](#preparing-the-input-excel-file)
- [Running the script](#running-the-script)
- [What the analysis does](#what-the-analysis-does)
- [Output files](#output-files)
- [Reproducibility](#reproducibility)
- [Interpreting the results](#interpreting-the-results)

## Requirements

- **MATLAB R2020a or later.** The script uses `exportgraphics` (R2020a), `sheetnames` (R2019b) and `readcell`/`writecell` (R2019a).
- **Statistics and Machine Learning Toolbox** may be needed for `quantile`, depending on your MATLAB release. Running `matlab.codetools.requiredFilesAndProducts('balanced_subsample_perm_ks_analysis_v2.m')` lists exactly what your installation needs.
- **Desktop MATLAB.** The script uses file-selection, checkbox and input dialogs, so it cannot run in `-batch` mode or without a display.
- **Excel is not required** to run the script, but output files must be closed in Excel before re-running. The script overwrites them.

No other files or toolboxes are needed. The K-S statistic and all corrections are implemented as local functions at the end of the script.

## Preparing the input Excel file

The input is a single `.xlsx` or `.xls` workbook.

| Excel element | Meaning |
|---|---|
| **Sheet** | One experimental **group**. The sheet name is used as the group name throughout the output. |
| **Column** | One **subject** in that group |
| **First non-empty cell of a column** | The subject's **name** |
| **Cells below the name** | That subject's numeric **observations**, one per cell |

Example: a sheet named `ACC`:

|   | A | B | C |
|---|---|---|---|
| 1 | Mouse01 | Mouse02 | Mouse03 |
| 2 | 0.412 | 0.398 | 0.455 |
| 3 | 0.387 | 0.421 | 0.430 |
| 4 | 0.441 | 0.405 |  |
| … | … | … | … |

with further sheets (for example `SAT1`, `PSE2`....) laid out the same way.

**Rules and things to watch for:**

- **Every column must have a name in its first non-empty cell.** The script always treats that cell as the name, even if it is a number. A column without a header would silently lose its first observation and be named after it.
- Columns can have **different lengths**. Subjects do not need the same number of observations.
- **Blank, text, `NaN` and `Inf` cells below the name are ignored.** Columns with no valid numbers are skipped.
- Each group needs **at least 2 subjects** with valid data. Groups with fewer are excluded automatically with a warning, because subject-level permutation is not possible with a single subject.
- Sheet names become part of output file names, so keep them short and descriptive.

## Running the script

Open the script in MATLAB and press **Run** (or type `balanced_subsample_perm_ks_analysis_v2` with its folder on the path). Then follow the dialogs:

1. **Select the input Excel file.**
2. **Choose the groups to compare.** Every sheet with at least 2 valid subjects is listed, along with its number of subjects and the smallest number of observations in any of its subjects. All are ticked by default; untick any you want to leave out. At least 2 must stay selected.
3. **Set the analysis parameters:**

   | Parameter | Meaning | Default |
   |---|---|---|
   | `N_SAMP` | Observations drawn from **each** subject in each iteration, without replacement. Must be between 1 and the smallest per-subject count among the selected groups. | That smallest count |
   | `N_ITER` | Number of balanced subsampling iterations on the real data | 10000 |
   | `N_PERM` | Number of subject permutations used to build the null distribution | 10000 |

4. **Set the CDF grid** (start, end and bin width) used for the mean-CDF curves and table. The defaults cover the full data range with roughly 100 bins. This grid affects **only** the CDF plots and the `Mean_CDF_by_Pair` sheet, not the K-S statistics or p-values.

The console then prints per-subject sample counts, progress bars appear during the three resampling stages, and a summary table is printed at the end.

Three settings are fixed at the top of the script and can be edited there:

| Setting | Meaning | Default |
|---|---|---|
| `ALPHA` | Significance level for all tests and corrections, and the width of the reported intervals (`1 − ALPHA`) | `0.05` |
| `PLOT_CORRECTION` | Which corrected p-value decides the `*` marks on the summary figure: `'BH'`, `'Holm'` or `'WY'` | `'BH'` |
| `SEED` | Random seed for reproducible results. `0` gives a different random stream each run. | `42` |

## What the analysis does

### 1. Balanced subsampling of the real data

This is repeated `N_ITER` times:

- From every subject in every selected group, draw `N_SAMP` observations at random **without replacement**.
- Pool the draws within each group. Each group now has `N_SAMP × number of subjects` values, and every subject contributes equally.
- For each pair of groups, compute the two-sample K-S statistic **D**, the largest vertical gap between the two groups' empirical cumulative distribution functions (CDFs). D ranges from 0 (identical distributions) to 1 (no overlap).
- Record each group's empirical CDF on the chosen grid.

The result is a distribution of `N_ITER` observed D values for each pair, whose mean is used as the test statistic. Averaging over many subsamples removes the dependence on any single random draw.

### 2. Subject-level permutation (null distribution)

This is repeated `N_PERM` times for each pair of groups:

- Pool the **subjects** of the two groups and randomly reassign them to two pseudo-groups of the original sizes.
- Draw `N_SAMP` observations per subject exactly as in step 1 and compute D.

This shows how large D gets when group labels carry no information, while keeping each subject's observations together.

### 3. p-values and multiple-comparison corrections

**Uncorrected p-value** for each pair:

```
p = (1 + number of null D values ≥ mean observed D) / (1 + N_PERM)
```

Adding 1 to the numerator and denominator keeps p from being exactly 0. The smallest attainable p is therefore `1 / (1 + N_PERM)`.

When 3 or more groups are selected (so there are several pairs), three corrections are applied across all pairs:

| Correction | Controls | How |
|---|---|---|
| **Benjamini-Hochberg (BH)** | False discovery rate | Standard step-up adjustment of the uncorrected p-values |
| **Holm-Bonferroni** | Family-wise error rate | Standard step-down adjustment of the uncorrected p-values |
| **Westfall-Young** (single-step min-P) | Family-wise error rate | A separate set of `N_PERM` **global** permutations shuffles subjects across *all* selected groups at once and computes D for every pair. Each D is converted to a p-value within its own pair, and the smallest p across pairs is kept for each permutation. A pair's corrected p is the fraction of permutations whose smallest p is at least as small as that pair's observed p. This is the nonparametric counterpart of Tukey's HSD and accounts for the correlation between comparisons that share a group. |

With only 2 groups there is a single comparison. All corrected p-values then equal the uncorrected one, and no global permutation is run.

## Output files

All outputs are saved in the **same folder as the input file**, with names built from the input file's name (`<name>` below). Re-running the script on the same input overwrites them.

### Excel workbook: `<name>_BalSubPerm_KS_results.xlsx`

| Sheet | Contents |
|---|---|
| `Parameters` | Everything needed to reproduce the run: input file, groups analysed, `N_SAMP`, `N_ITER`, `N_PERM`, permutation scheme, p-value formula, corrections used, `ALPHA`, random seed, CDF grid, and run date |
| `Subject_Statistics` | One row per subject (see below) |
| `KS_Summary` | One row per pair of groups, with the main results (see below) |
| `Observed_D_Iterations` | The full `N_ITER × pairs` matrix of observed D values. Column 1 is the iteration number, and each further column is one pair. |
| `Permutation_D_Null` | The full `N_PERM × pairs` matrix of null D values from the pairwise permutations, in the same layout |
| `Global_Permutation_D_Null` | *(3+ groups only)* The `N_PERM × pairs` matrix of D values from the global permutations used for the Westfall-Young correction |
| `Mean_CDF_by_Pair` | The mean CDF of each group for every pair (see below) |

**`Subject_Statistics` columns:**

| Column | Meaning |
|---|---|
| `Group`, `Subject` | Sheet and subject name |
| `N_Samples` | Number of valid observations for the subject |
| `Fraction_Drawn_Per_Iter` | `N_SAMP / N_Samples`: the share of this subject's data used in each iteration. Subjects at `1` contribute all of their data in every iteration. |
| `Min`, `Max`, `Mean`, `Median` | Summary of the subject's observations |
| `Std_Dev`, `SEM` | Standard deviation and standard error of the mean |
| `CV_pct` | Coefficient of variation, `100 × Std_Dev / Mean`. This is meaningless if the mean is near 0 or the data can be negative. |

**`KS_Summary` columns:**

| Column | Meaning |
|---|---|
| `Comparison`, `Group_1`, `Group_2` | The pair compared |
| `N_Subjects_G1`, `N_Subjects_G2` | Number of subjects in each group |
| `N_Distinct_Subject_Assignments` | Number of distinct ways to split the pair's subjects into two groups of these sizes. Below about 1000 (also flagged in the console), the permutation null has few distinct values and the smallest attainable p-value is limited by the number of subjects, not by `N_PERM`. |
| `Obs_Mean_D`, `Obs_Median_D`, `Obs_SD_D` | Mean, median and SD of the observed D across iterations |
| `Obs_D_Lo_95Pct`, `Obs_D_Hi_95Pct` | Central 95 % interval of the observed D across iterations (the percentage follows `ALPHA`). This shows how much D depends on which observations were drawn, not a confidence interval for the population D. |
| `Null_Mean_D`, `Null_Median_D` | Mean and median of the permutation null D |
| `Null_Critical_D_95thPct` | The 95th percentile of the null D. An observed mean D above this value corresponds to an uncorrected p below about 0.05. |
| `p_uncorrected`, `Sig_uncorrected` | Uncorrected permutation p-value, and `*` if `p ≤ ALPHA` |
| `p_BH_corrected`, `Sig_BH_corrected` | Benjamini-Hochberg adjusted p-value, and `*` if significant |
| `p_Holm_corrected`, `Sig_Holm_corrected` | Holm-Bonferroni adjusted p-value, and `*` if significant |
| `p_WestfallYoung_corrected`, `Sig_WestfallYoung_corrected` | Westfall-Young adjusted p-value, and `*` if significant |

**`Mean_CDF_by_Pair` layout:** each pair occupies a block of 10 columns, and the blocks sit side by side, separated by one blank column. Within each block:

- **Row 1** holds the pair's name.
- **Row 2** holds the column headers.
- **Row 3 onward** has one row per CDF grid point.

| Column in block | Meaning |
|---|---|
| `Bin_End` | Grid value *x* |
| `<G1>_MeanCDF`, `<G1>_SD` | Group 1's mean and SD of P(value ≤ *x*) across iterations |
| `<G1>_Lo_95Pct`, `<G1>_Hi_95Pct` | Central 95 % interval of group 1's CDF across iterations |
| `<G2>_…` | The same four columns for group 2 |
| `MeanCDF_Diff_G1_minus_G2` | Difference between the two mean CDFs at each grid point |

### Figures (PNG, 300 dpi)

Each pair name is turned into a file-name tag by replacing spaces and punctuation with `_`. For example, `Control vs Treatment_A` becomes `Control_vs_Treatment_A`.

| File | Contents |
|---|---|
| `<name>_BalSubPerm_KS_ObsVsNull_<pair>.png` | One per pair. Histograms of the observed D values and the permutation null D values. A red line marks the mean observed D and a dashed line marks the null 95th percentile. The title shows all four p-values. |
| `<name>_BalSubPerm_KS_MeanCDF_<pair>.png` | One per pair. The mean CDF of each group, with shaded 95 % intervals across iterations. |
| `<name>_BalSubPerm_KS_Summary.png` | All pairs on one figure: bars show the mean observed D with 95 % interval whiskers, and diamonds show the null 95th percentile. `*` marks pairs that are significant under `PLOT_CORRECTION`. |

## Reproducibility

- With `SEED` set to a non-zero value, re-running the script on the **same input, group selection and parameters** in the same MATLAB release gives identical results. Changing any of these, including the order of the sheets in the workbook, changes the random draws.
- The `Parameters` sheet records every setting used. Keep it alongside any reported result.
- Results depend on `N_SAMP`. Using the default (the largest allowed value) uses as much data as possible from the smallest subject, but subjects with exactly `N_SAMP` observations then contribute all of their data in every iteration. The console notes how many subjects this applies to.

## Interpreting the results

- **D measures a difference in the whole distribution,** not just in the mean. A significant result can reflect a shift, a change in spread, or a change in shape. Use the mean-CDF plots to see which.
- **The p-value tests whether groups differ by more than reassigning subjects at random would produce.** With few subjects per group, even large differences may not reach significance, and `N_Distinct_Subject_Assignments` shows how fine-grained the null can be.
- **Which correction to report:** BH (the default for the plot) controls the expected proportion of false positives among significant results. Holm and Westfall-Young control the chance of *any* false positive, and are more conservative. Westfall-Young is usually more powerful than Holm, because it accounts for correlations between comparisons that share a group. Decide which to report before looking at the results.

## In this study

[USED FOR ALL THE SYNAPSE CDH ANALYSES IN FIGURES 2, 3, 5, 7, S1, and S3]

## Version history

- **v1:** balanced within-subject subsampling, subject-level permutation K-S test, uncorrected and Benjamini-Hochberg p-values.
- **v2:** adds Holm-Bonferroni and Westfall-Young (min-P, global subject permutation) corrections, and the `PLOT_CORRECTION` option. 
