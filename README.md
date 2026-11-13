# Analysis code

Code for the analyses reported in the manuscript on identifying (same-day) and
forecasting (next-day) clinically meaningful pain increases among veterans from
daily diary and wearable data (MAVERICK study). Scripts are numbered in the order
they are run.

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `01_imputation.R` | Sets physiologically implausible wearable values to missing, applies the missingness exclusions, and produces 50 imputed datasets by multilevel multiple imputation by chained equations | `<DATA_DIR>/imputed_datasets/` |
| 2 | `02_gpboost_pipeline.R` (run by `run_models.sh`) | Fits one of the four models: builds the analytic dataset in every imputation, splits each participant's days 80/20, tunes and fits a GPBoost model per imputation, saves test-day predictions, pools gain importance, selects features by scree-plot elbow, computes SHAP values | `results/<model>/` |
| 3 | `03_performance_tables.R` | Pooled performance (AUC; tree-only and intercept-only AUC; PR-AUC; sensitivity, specificity and precision) with 95% CIs whose within-imputation variance comes from a participant-level bootstrap (participants resampled with replacement, 2,000 times per imputation), within-person AUC, and sample sizes (Tables 2 and 3) | `results/tables/` |
| 4 | `04_shap_direction.R` | Spearman correlation between feature value and SHAP value for each retained feature (Supplementary Tables 2 and 3) | `results/<model>/`, `results/tables/` |
| 5 | `05_figures.R` | SHAP beeswarm plots (Figure 1, Supplementary Figure 1) | `results/figures/` |
| 6 | `06_table1.R` | Descriptive statistics of the analytic sample from the observed data (Table 1, Supplementary Table 1) | `results/tables/` |

## Running

From this folder:

```bash
Rscript 01_imputation.R      # raw data -> 50 imputed datasets
./run_all.sh                 # steps 2 to 6
```

`run_all.sh` fits the four models one after another (`run_models.sh`) and then
runs steps 3 to 6. A single model can be fitted with, for example,
`HORIZON=nextday FEATURES=wearable Rscript 02_gpboost_pipeline.R`.

Environment variables (all optional):

- `DATA_DIR` (default `data`): folder holding `mavdata.csv` and `imputed_datasets/`.
- `RESULTS_DIR` (default `results`): where outputs are written.
- `N_BOOT` (default `2000`): participant resamples per imputation for the confidence intervals of step 3.

## The four models

Models are named `<horizon>_<features>`:

- `sameday_full`, `sameday_wearable`: same-day identification (a day's features and that day's pain outcome).
- `nextday_full`, `nextday_wearable`: next-day forecasting (the prior day's features and today's pain outcome).
- `full` = daily diary and wearable features; `wearable` = wearable features only. Both include the baseline survey characteristics.

## Output files per model (`results/<model>/`)

| File | Contents |
|---|---|
| `<model>_analytic_sample_per_imputation.csv` | Participants, person-days and pain-increase days in each imputed dataset |
| `<model>_tuned_hyperparameters.csv` | Selected hyperparameters, boosting rounds, CV log-loss and Youden threshold per imputation |
| `<model>_model_imp<i>.json` | The fitted GPBoost model for imputation i |
| `<model>_test_predictions.csv` | Held-out test days: outcome and predicted probability from the full model, the tree component alone and the random intercept alone, per imputation |
| `<model>_gain_importance.csv`, `<model>_variable_selection.csv`, `<model>_scree.pdf` | Pooled gain importance, scree-plot elbow and the retained features |
| `<model>_shap_per_imputation.rds` | SHAP values, feature values and person-day keys for the test days of every imputation |
| `<model>_shap_importance.csv` | Mean absolute SHAP value per feature after averaging across imputations within person-day |
| `<model>_within_person_auc.csv`, `<model>_shap_direction.csv` | Written by steps 3 and 4 |

Feature names carry the daily variable's name plus a suffix: `_absdelta` is the
daily deviation from the participant's reference-period mean, `_baseline_mean`
and `_baseline_sd` are the reference-period (days 1 to 14) mean and SD. Baseline
survey characteristics end in `_0`. Wearable variables start with `d.`.

## Data

The data are not included in this folder. `mavdata.csv` holds one row per
participant-day with the daily diary, wearable and baseline survey variables.

## Software

R 4.5.2 with gpboost 1.6.6, pROC 1.19.0.1, PRROC 1.4, dplyr 1.2.0, tidyr 1.3.2,
readr 2.2.0, ggplot2 4.0.2, ggbeeswarm 0.7.3 and patchwork 1.3.2. The imputation
(November 2025) used the mice and miceadds packages.

## Run time

Each model takes about 4 to 6 hours on an 8-core laptop with 16 GB of memory
(hyperparameter tuning across the 50 imputations accounts for nearly all of it;
the pipeline uses four parallel R processes with two threads each). Steps 3 to 6
take a few minutes. All random seeds are fixed.

## License and citation

The code is released under the MIT License (see `LICENSE`). Please cite the
article and the archived release of this repository (https://github.com/sks9496/maverick-pain-who-vs-when;
Zenodo DOI: https://doi.org/10.5281/zenodo.24184005). Citation metadata is in `CITATION.cff`.

## Data availability

The MAVERICK data are not distributed with this code; see the article's Data
Availability Statement for the access route.
