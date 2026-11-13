# ===================================================================
# GPBoost pipeline for same-day identification and next-day forecasting of
# clinically meaningful pain increases from daily diary and wearable data.
#
# One run fits ONE of the four reported models, chosen by two environment
# variables:
#   HORIZON  = sameday | nextday    same-day identification or next-day forecast
#   FEATURES = full    | wearable   daily diary + wearable features, or wearable
#                                   features only (both feature sets also include
#                                   the baseline survey characteristics)
#
# Steps:
#   1. Load the 50 imputed datasets; compute each participant's reference-period
#      (days 1-14) mean and SD of every daily variable.
#   2. Build the analytic dataset for each imputation: daily deviation features,
#      the pain-increase outcome and, for next-day models, prior-day features.
#   3. For each imputation: random 80/20 within-person train/test split,
#      hyperparameter tuning by 10-fold cross-validation on the training days,
#      and a Youden's J classification threshold from the out-of-fold predictions.
#   4. For each imputation: fit the tuned model on the training days, save it,
#      and save its predictions for the held-out test days.
#   5. Pool gain-based feature importance across imputations and retain the
#      features above the scree-plot elbow.
#   6. Compute SHAP values on the held-out test days, per imputation and pooled.
#
# Performance metrics and the manuscript tables are computed from the saved
# predictions by 03_performance_tables.R; figures by 05_figures.R.
#
# Usage (from this folder):
#   HORIZON=sameday FEATURES=full Rscript 02_gpboost_pipeline.R
#   ./run_models.sh            # all four models in sequence
# Optional environment variables: DATA_DIR (folder that holds imputed_datasets/;
# default "data") and RESULTS_DIR (default "results").
# ===================================================================

Sys.setenv("R_MAX_VSIZE" = "16Gb")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(pROC)
  library(gpboost)
  library(parallel)
})

# -------------------------------------------------------------------
# 1. Configuration
# -------------------------------------------------------------------
HORIZON  <- Sys.getenv("HORIZON",  unset = "sameday")
FEATURES <- Sys.getenv("FEATURES", unset = "full")
stopifnot(HORIZON %in% c("sameday", "nextday"), FEATURES %in% c("full", "wearable"))

model_name  <- paste0(HORIZON, "_", FEATURES)
data_dir    <- file.path(Sys.getenv("DATA_DIR", unset = "data"), "imputed_datasets")
results_dir <- file.path(Sys.getenv("RESULTS_DIR", unset = "results"), model_name)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
results_dir <- normalizePath(results_dir)
out_file <- function(suffix) file.path(results_dir, paste0(model_name, "_", suffix))

M                  <- 50      # number of imputed datasets
reference_period   <- 1:14    # days used for each participant's reference mean and SD
min_reference_pain <- 0.5     # participants below this reference-period mean pain are excluded
train_prop         <- 0.80    # share of each participant's days used for training
tune_n_folds       <- 10      # cross-validation folds for tuning
tune_n_combos      <- 30      # hyperparameter combinations tried per imputation
n_workers          <- 4       # parallel R processes (each limited to 2 OpenMP threads)
seed               <- 123

cat(sprintf("\nModel: %s (M = %d imputations)\nResults: %s\n", model_name, M, results_dir))

# -------------------------------------------------------------------
# 2. Variables and feature sets
# -------------------------------------------------------------------
outcome_var <- "pain_daily"

wearable_daily_vars <- c(
  "d.distance_sum_value", "d.steps_sum_value", "d.stairs_sum_value",
  "d.sleep.efficiency_mean_value",
  "d.sleep.minutes.tosleep_mean_value", "d.sleep.minutes.afterwake_sum_value",
  "d.sleep.minutes.asleep_sum_value", "d.sleep.timein.bed_sum_value",
  "d.sleep.minutes.awake_sum_value", "d.sleep.rem.minutes_sum_value"
)

diary_daily_vars <- c(
  "sleephours_daily", "sleepquality_daily", "stress_daily",
  "dep_daily", "anx_daily", "pcl4tot_daily",
  "panaspos_daily", "panasneg_daily", "cigquant_daily_sum",
  "vapequant_daily_sum", "canhours_daily_sum", "alcdrinks_daily_sum", "SSS_daily"
)

baseline_survey_vars <- c(
  "minority_0", "phqcase_0", "gadcase_0", "bmacq_0", "tbiscreen_0",
  "age_0", "deploy_0", "lectot_0", "isi_0", "pcl_0", "ucla_0",
  "pain_0", "uppstot_0", "mmq_0", "passtot_0", "pctitot_0",
  "alcbinge_0", "alcalone_0", "alcbingealone_0", "alcdays_0",
  "mjdays_0", "mjalone_0", "woman_0"
)

# Each daily variable yields three features, named with these suffixes:
#   _absdelta       daily deviation from the participant's reference-period mean
#   _baseline_mean  reference-period mean (the participant's typical level)
#   _baseline_sd    reference-period SD (the participant's typical variability)
# Next-day models also use the prior-day pain deviation as a feature, so pain
# joins the list of variables that get a deviation feature there.
deviation_vars <- c(wearable_daily_vars, diary_daily_vars)
if (HORIZON == "nextday") deviation_vars <- c(deviation_vars, outcome_var)
reference_vars <- unique(c(deviation_vars, outcome_var))

feature_block <- function(vars) {
  c(paste0(vars, "_absdelta"), paste0(vars, "_baseline_mean"), paste0(vars, "_baseline_sd"))
}
features <- if (FEATURES == "full") {
  c(feature_block(wearable_daily_vars),
    feature_block(diary_daily_vars),
    if (HORIZON == "nextday") "pain_daily_absdelta",
    "pain_daily_baseline_mean", "pain_daily_baseline_sd",
    baseline_survey_vars)
} else {
  c(feature_block(wearable_daily_vars), baseline_survey_vars)
}
cat(sprintf("Features: %d\n\n", length(features)))

# -------------------------------------------------------------------
# 3. Helper functions
# -------------------------------------------------------------------

# Reference-period mean and SD of every daily variable, per participant.
# A participant with a single reference-period value gets SD 0.
reference_stats <- function(imputed) {
  imputed %>%
    filter(day_number %in% reference_period) %>%
    group_by(uid) %>%
    summarise(across(all_of(reference_vars),
                     list(baseline_mean = ~ mean(.x, na.rm = TRUE),
                          baseline_sd   = ~ sd(.x, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
    mutate(across(ends_with("_baseline_sd"), ~ replace_na(.x, 0)))
}

# Daily deviation = daily value minus the participant's reference-period mean.
add_deviation_features <- function(df, ref_stats) {
  df <- df %>% left_join(ref_stats, by = "uid")
  for (v in deviation_vars) {
    df[[paste0(v, "_absdelta")]] <- df[[v]] - df[[paste0(v, "_baseline_mean")]]
  }
  df %>% mutate(across(where(is.numeric), ~ ifelse(is.finite(.x), .x, NA_real_)))
}

# Next-day models: shift every deviation feature to the prior day within
# participant. The first day per participant and any day that follows a gap in
# consecutive days are dropped, so "prior day" always means yesterday.
shift_to_prior_day <- function(df) {
  deviation_cols <- paste0(deviation_vars, "_absdelta")
  df %>%
    arrange(uid, day_number) %>%
    group_by(uid) %>%
    mutate(across(all_of(deviation_cols), ~ dplyr::lag(.x, n = 1)),
           .day_gap = day_number - dplyr::lag(day_number)) %>%
    filter(.day_gap == 1) %>%
    select(-.day_gap) %>%
    ungroup()
}

# Random 80/20 split of each participant's days into training and test days.
split_within_person <- function(df, split_seed) {
  set.seed(split_seed)
  df <- df %>%
    group_by(uid) %>%
    mutate(.n_days = n(),
           .is_train = sample(rep(c(TRUE, FALSE),
                                  times = c(floor(train_prop * .n_days[1]),
                                            .n_days[1] - floor(train_prop * .n_days[1]))))) %>%
    ungroup()
  list(train = df %>% filter(.is_train)  %>% select(-starts_with(".")),
       test  = df %>% filter(!.is_train) %>% select(-starts_with(".")))
}

# Person-by-class observation weights: each outcome class carries equal total
# weight, and within a class every participant's days carry equal total weight.
# Weights are normalised to a mean of 1 and passed to GPModel(), so the tree
# component and the random intercept are both estimated under the weighted loss.
make_personclass_weights <- function(uid, y) {
  df <- data.frame(uid = uid, y = y, .row = seq_along(y))
  cnt <- df %>% group_by(uid, y) %>% summarise(n = n(), .groups = "drop")
  per_class_persons <- cnt %>% group_by(y) %>% summarise(np = n(), .groups = "drop")
  cnt <- cnt %>% left_join(per_class_persons, by = "y") %>%
    mutate(w_cell = 1 / (n * np))
  df <- df %>% left_join(cnt %>% select(uid, y, w_cell), by = c("uid", "y"))
  df$w_cell <- df$w_cell / mean(df$w_cell, na.rm = TRUE)
  df$w_cell[order(df$.row)]
}

# Run a per-imputation function on a cluster of separate R processes. The
# objects named in export_names are copied from this session to every worker.
run_parallel <- function(imps, worker_fn, export_names) {
  cl <- makePSOCKcluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "2")
    suppressPackageStartupMessages({ library(dplyr); library(pROC); library(gpboost) })
  })
  clusterExport(cl, export_names, envir = globalenv())
  clusterSetRNGStream(cl, iseed = seed)
  parLapply(cl, imps, worker_fn)
}

# -------------------------------------------------------------------
# 4. Load the imputed datasets and build the analytic dataset per imputation
# -------------------------------------------------------------------
cat("Building analytic datasets...\n")
model_data  <- vector("list", M)
sample_rows <- vector("list", M)
for (i in 1:M) {
  imputed <- read.csv(file.path(data_dir, sprintf("imputed_dataset_%d.csv", i)))
  d <- imputed %>%
    filter(day_number > max(reference_period)) %>%         # days after the reference period
    add_deviation_features(reference_stats(imputed)) %>%
    filter(pain_daily_baseline_mean >= min_reference_pain) %>%
    # Outcome: daily pain at least 2 points above the participant's reference-period
    # mean, with the mean rounded to the nearest 0.5 (the scale's resolution).
    mutate(pain_reference_rounded = round(pain_daily_baseline_mean * 2) / 2,
           Y = as.integer(pain_daily >= pain_reference_rounded + 2)) %>%
    filter(!is.na(Y))
  if (HORIZON == "nextday") d <- shift_to_prior_day(d)

  model_data[[i]]  <- d
  sample_rows[[i]] <- data.frame(imputation = i, n_participants = n_distinct(d$uid),
                                 n_person_days = nrow(d), n_pain_increase_days = sum(d$Y))
  if (i %% 10 == 0) cat(sprintf("  %d/%d\n", i, M))
}

missing_features <- setdiff(features, names(model_data[[1]]))
if (length(missing_features) > 0) {
  stop("Features missing from the data: ", paste(missing_features, collapse = ", "))
}
sample_sizes <- bind_rows(sample_rows)
write_csv(sample_sizes, out_file("analytic_sample_per_imputation.csv"))
cat(sprintf("Analytic person-days per imputation: %d-%d (%d-%d participants)\n\n",
            min(sample_sizes$n_person_days), max(sample_sizes$n_person_days),
            min(sample_sizes$n_participants), max(sample_sizes$n_participants)))

# -------------------------------------------------------------------
# 5. Hyperparameter tuning: per imputation, 10-fold CV on the training days
# -------------------------------------------------------------------
# Random search over 30 of the 180 combinations below (the same 30 for every
# imputation). Maximum tree depth is unconstrained and line-search step length
# is enabled in every fit.
param_grid <- expand.grid(
  learning_rate    = c(0.001, 0.01, 0.05, 0.1, 1),
  min_data_in_leaf = c(10, 20, 50),
  lambda_l2        = c(0, 1, 10),
  num_leaves       = c(7, 15, 31, 63),
  stringsAsFactors = FALSE
)
set.seed(seed)
grid_sample  <- param_grid[sample(nrow(param_grid), tune_n_combos), , drop = FALSE]
fixed_params <- list(objective = "binary", max_depth = -1, line_search_step_length = TRUE, verbose = -1)

# For one imputation: split off the training days, then for each candidate
# combination run 10-fold CV with early stopping (up to 1,000 boosting rounds,
# stopping after 10 rounds without improvement in validation log-loss). The
# combination with the lowest mean validation log-loss wins; its number of
# boosting rounds is the mean best iteration over the folds; the classification
# threshold is Youden's J applied to its out-of-fold predictions.
tune_one_imputation <- function(imp) {
  train_d   <- split_within_person(as.data.frame(model_data[[imp]]), split_seed = seed + imp * 1000)$train
  X_train   <- as.matrix(train_d[, features, drop = FALSE])
  y_train   <- as.numeric(train_d$Y)
  uid_train <- train_d$uid
  w_train   <- make_personclass_weights(uid_train, y_train)

  set.seed(seed + imp)
  fold_ids <- sample(rep(1:tune_n_folds, length.out = length(y_train)))

  n_try         <- nrow(grid_sample)
  combo_scores  <- numeric(n_try)
  combo_nrounds <- integer(n_try)
  combo_oof     <- vector("list", n_try)

  for (ci in 1:n_try) {
    combo        <- grid_sample[ci, ]
    fold_losses  <- numeric(tune_n_folds)
    fold_nrounds <- integer(tune_n_folds)
    oof_preds    <- rep(NA_real_, length(y_train))
    combo_ok     <- TRUE

    for (k in 1:tune_n_folds) {
      train_k <- fold_ids != k
      val_k   <- fold_ids == k

      gp_k <- GPModel(group_data = uid_train[train_k], likelihood = "bernoulli_logit",
                      weights = w_train[train_k])
      gp_k$set_prediction_data(group_data_pred = uid_train[val_k])

      dtrain_k <- gpb.Dataset(data = X_train[train_k, , drop = FALSE], label = y_train[train_k])
      dval_k   <- gpb.Dataset(data = X_train[val_k, , drop = FALSE],   label = y_train[val_k])

      fold_seed <- seed + imp * 100 + k
      set.seed(fold_seed)
      bst_k <- tryCatch(
        suppressWarnings(
          gpb.train(data = dtrain_k, gp_model = gp_k, nrounds = 1000,
                    params = c(fixed_params, as.list(combo), list(seed = fold_seed)),
                    valids = list(val = dval_k),
                    use_gp_model_for_validation = TRUE,
                    early_stopping_rounds = 10, verbose = 0)
        ), error = function(e) NULL)
      if (is.null(bst_k)) { combo_ok <- FALSE; break }

      fold_nrounds[k] <- bst_k$best_iter
      fold_losses[k]  <- bst_k$best_score

      pred_k <- tryCatch(
        suppressWarnings(
          predict(bst_k, data = X_train[val_k, , drop = FALSE],
                  group_data_pred = uid_train[val_k], pred_latent = FALSE)
        ), error = function(e) NULL)
      if (!is.null(pred_k)) oof_preds[val_k] <- as.numeric(pred_k$response_mean)

      rm(bst_k, gp_k, dtrain_k, dval_k); gc(verbose = FALSE)
    }

    if (combo_ok) {
      combo_scores[ci]  <- mean(fold_losses, na.rm = TRUE)
      combo_nrounds[ci] <- max(1, round(mean(fold_nrounds)))
    } else {
      combo_scores[ci]  <- Inf      # a combination that fails in any fold is not eligible
    }
    combo_oof[[ci]] <- oof_preds
  }

  if (!is.finite(min(combo_scores))) stop("Every hyperparameter combination failed for imputation ", imp)
  best_ci   <- which.min(combo_scores)
  oof_best  <- combo_oof[[best_ci]]
  valid_oof <- !is.na(oof_best)
  threshold <- 0.5
  if (sum(valid_oof) >= 10 && length(unique(y_train[valid_oof])) >= 2) {
    roc_obj   <- roc(response = y_train[valid_oof], predictor = oof_best[valid_oof],
                     levels = c(0, 1), direction = "<", quiet = TRUE)
    threshold <- coords(roc_obj, "best", best.method = "youden", ret = "all")$threshold[1]
  }

  list(imputation = imp, params = as.list(grid_sample[best_ci, ]), nrounds = combo_nrounds[best_ci],
       cv_logloss = combo_scores[best_ci], threshold = threshold)
}

cat(sprintf("Tuning: %d combinations x %d folds x %d imputations on %d workers...\n",
            tune_n_combos, tune_n_folds, M, n_workers))
tuned_list <- run_parallel(
  1:M, tune_one_imputation,
  export_names = c("model_data", "features", "grid_sample", "fixed_params", "seed",
                   "train_prop", "tune_n_folds", "split_within_person", "make_personclass_weights"))

hyperparameters <- bind_rows(lapply(tuned_list, function(t) {
  data.frame(imputation = t$imputation, learning_rate = t$params$learning_rate,
             num_leaves = t$params$num_leaves, min_data_in_leaf = t$params$min_data_in_leaf,
             lambda_l2 = t$params$lambda_l2, nrounds = t$nrounds,
             cv_logloss = t$cv_logloss, threshold = t$threshold)
}))
write_csv(hyperparameters, out_file("tuned_hyperparameters.csv"))
cat(sprintf("  Thresholds: %.3f-%.3f (mean %.3f)\n\n", min(hyperparameters$threshold),
            max(hyperparameters$threshold), mean(hyperparameters$threshold)))

# -------------------------------------------------------------------
# 6. Fit the tuned model per imputation and predict the held-out test days
# -------------------------------------------------------------------
fit_one_imputation <- function(imp) {
  sp      <- split_within_person(as.data.frame(model_data[[imp]]), split_seed = seed + imp * 1000)
  train_d <- sp$train
  test_d  <- sp$test
  X_train <- as.matrix(train_d[, features, drop = FALSE])
  y_train <- as.numeric(train_d$Y)
  X_test  <- as.matrix(test_d[, features, drop = FALSE])
  w_train <- make_personclass_weights(train_d$uid, y_train)
  tuned   <- tuned_list[[imp]]

  gp_model <- GPModel(group_data = train_d$uid, likelihood = "bernoulli_logit", weights = w_train)
  set.seed(seed + imp)
  bst <- suppressWarnings(
    gpboost(data = X_train, label = y_train, gp_model = gp_model, nrounds = tuned$nrounds,
            params = c(fixed_params, tuned$params, list(seed = seed + imp)), verbose = 0))

  # Predicted probability of a pain increase from the full model, from the tree
  # component alone (random intercept set to zero) and from the random intercept alone.
  pred_full   <- suppressWarnings(predict(bst, data = X_test, group_data_pred = test_d$uid, pred_latent = FALSE))
  pred_latent <- suppressWarnings(predict(bst, data = X_test, group_data_pred = test_d$uid, pred_latent = TRUE))
  predictions <- data.frame(
    uid = test_d$uid, day_number = test_d$day_number, Y = as.numeric(test_d$Y),
    full_prob = as.numeric(pred_full$response_mean),
    tree_prob = plogis(as.numeric(pred_latent$fixed_effect)),
    re_prob   = plogis(as.numeric(pred_latent$random_effect_mean)),
    imputation = imp)

  gain <- gpb.importance(bst, percentage = TRUE)
  importance <- data.frame(variable = gain$Feature, importance = gain$Gain, imputation = imp)

  suppressWarnings(gpb.save(bst, filename = out_file(sprintf("model_imp%d.json", imp))))
  list(predictions = predictions, importance = importance)
}

cat(sprintf("Fitting %d models on %d workers...\n", M, n_workers))
fit_results <- run_parallel(
  1:M, fit_one_imputation,
  export_names = c("model_data", "features", "tuned_list", "fixed_params", "seed", "train_prop",
                   "split_within_person", "make_personclass_weights", "out_file",
                   "results_dir", "model_name"))

predictions <- bind_rows(lapply(fit_results, `[[`, "predictions"))
write_csv(predictions, out_file("test_predictions.csv"))
cat(sprintf("  Saved %d test-day predictions (%d distinct person-days)\n\n",
            nrow(predictions), nrow(distinct(predictions, uid, day_number))))

# -------------------------------------------------------------------
# 7. Pooled gain importance and scree-plot feature selection
# -------------------------------------------------------------------
# Gain importance (percentage of the training-loss reduction attributable to
# splits on the feature) averaged over the 50 imputations; a feature never split
# on in an imputation contributes 0 there.
importance_by_imp <- bind_rows(lapply(fit_results, `[[`, "importance"))
pooled_importance <- expand_grid(variable = unique(importance_by_imp$variable), imputation = 1:M) %>%
  left_join(importance_by_imp, by = c("variable", "imputation")) %>%
  mutate(importance = replace_na(importance, 0)) %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance), sd_importance = sd(importance), .groups = "drop") %>%
  arrange(desc(mean_importance))
write_csv(pooled_importance, out_file("gain_importance.csv"))

# Elbow of the ranked-importance curve: the rank with the greatest perpendicular
# distance from the straight line joining its first and last points (both axes
# rescaled to 0-1). Features at or above the elbow are retained.
n_vars <- nrow(pooled_importance)
x <- 1:n_vars
y <- pooled_importance$mean_importance
x_n <- (x - min(x)) / (max(x) - min(x))
y_n <- (y - min(y)) / (max(y) - min(y))
line_vec  <- c(x_n[n_vars] - x_n[1], y_n[n_vars] - y_n[1])
perp_dist <- abs(line_vec[2] * (x_n - x_n[1]) - line_vec[1] * (y_n - y_n[1])) / sqrt(sum(line_vec^2))
n_retain  <- which.max(perp_dist)

write_csv(data.frame(rank = x, variable = pooled_importance$variable, mean_importance = y,
                     selected = x <= n_retain),
          out_file("variable_selection.csv"))
pdf(out_file("scree.pdf"), width = 10, height = 6)
plot(x, y, type = "b", pch = 19, cex = 0.5, xlab = "Feature rank", ylab = "Mean gain importance (%)",
     main = sprintf("%s: %d of %d features retained", model_name, n_retain, n_vars))
abline(v = n_retain, col = "red", lty = 2)
invisible(dev.off())
cat(sprintf("Scree elbow: %d of %d features retained\n\n", n_retain, n_vars))

# -------------------------------------------------------------------
# 8. SHAP values on the held-out test days
# -------------------------------------------------------------------
# For each imputation, reload the saved model and compute SHAP values (LightGBM
# predcontrib) for that imputation's test days, taken from the saved
# predictions. A row's SHAP values sum to the tree component's log-odds, so
# every imputation is checked against its saved tree probability. The values are
# then averaged across imputations within person-day, and each feature's SHAP
# importance is the mean absolute pooled value.
cat("SHAP values...\n")
shap_list <- vector("list", M)
feat_list <- vector("list", M)
key_list  <- vector("list", M)
for (imp in 1:M) {
  bst    <- suppressWarnings(gpb.load(filename = out_file(sprintf("model_imp%d.json", imp))))
  keys   <- predictions %>% filter(imputation == imp) %>% select(uid, day_number, tree_prob)
  test_d <- merge(keys, as.data.frame(model_data[[imp]]), by = c("uid", "day_number"), sort = FALSE)
  stopifnot(nrow(test_d) == nrow(keys))
  X_test <- as.matrix(test_d[, features, drop = FALSE])

  shap_raw  <- suppressWarnings(predict(bst, data = X_test, predcontrib = TRUE))
  prob_diff <- max(abs(plogis(rowSums(shap_raw)) - test_d$tree_prob))
  if (!is.finite(prob_diff) || prob_diff > 1e-6) {
    stop(sprintf("SHAP check failed for imputation %d (max probability difference %.2e)", imp, prob_diff))
  }
  shap_mat <- shap_raw[, 1:(ncol(shap_raw) - 1), drop = FALSE]   # last column is the bias term
  colnames(shap_mat) <- features

  shap_list[[imp]] <- shap_mat
  feat_list[[imp]] <- X_test
  key_list[[imp]]  <- test_d[, c("uid", "day_number")]
  rm(bst); gc(verbose = FALSE)
}
saveRDS(list(shap_list = shap_list, feat_list = feat_list, key_list = key_list,
             features = features, model_name = model_name),
        out_file("shap_per_imputation.rds"))

pooled_shap <- bind_rows(lapply(1:M, function(i) {
  as.data.frame(shap_list[[i]]) %>%
    mutate(uid = key_list[[i]]$uid, day_number = key_list[[i]]$day_number)
})) %>%
  group_by(uid, day_number) %>%
  summarise(across(all_of(features), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

shap_importance <- data.frame(
  variable      = features,
  mean_abs_shap = colMeans(abs(as.matrix(pooled_shap[, features]))),
  n_imputations = M
) %>% arrange(desc(mean_abs_shap))
write_csv(shap_importance, out_file("shap_importance.csv"))
cat(sprintf("  Pooled over %d imputations: %d distinct person-days\n", M, nrow(pooled_shap)))

cat(sprintf("\n%s complete. Outputs in %s\n", model_name, results_dir))
