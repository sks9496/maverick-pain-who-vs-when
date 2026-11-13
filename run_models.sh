#!/usr/bin/env bash
# Fit the four reported models one after another. Each model's console output
# goes to <RESULTS_DIR>/<model>/run_log.txt. Run from this folder:
#   ./run_models.sh
# DATA_DIR and RESULTS_DIR are passed through to 02_gpboost_pipeline.R.
cd "$(dirname "$0")"
RESULTS_DIR="${RESULTS_DIR:-results}"

for horizon in sameday nextday; do
  for features in full wearable; do
    model="${horizon}_${features}"
    mkdir -p "$RESULTS_DIR/$model"
    echo "=== $model  start $(date)"
    HORIZON=$horizon FEATURES=$features RESULTS_DIR="$RESULTS_DIR" \
      Rscript 02_gpboost_pipeline.R > "$RESULTS_DIR/$model/run_log.txt" 2>&1
    echo "=== $model  done  $(date)"
  done
done
