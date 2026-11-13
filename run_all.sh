#!/usr/bin/env bash
# Reproduce the analyses from the imputed datasets: fit the four models, then
# compute the performance tables, SHAP directions, figures and Table 1.
# Run from this folder:  ./run_all.sh
# (01_imputation.R, which produces the imputed datasets from the raw data, is
# run separately: Rscript 01_imputation.R)
set -e
cd "$(dirname "$0")"
./run_models.sh
Rscript 03_performance_tables.R
Rscript 04_shap_direction.R
Rscript 05_figures.R
Rscript 06_table1.R
