# ===================================================================
# Model performance (Tables 2 and 3) from the saved test-day predictions.
# No model is refit.
#
# For each of the four models and each imputation, using that imputation's
# tuned classification threshold:
#   - AUC of the full model, of the tree component alone (random intercept set
#     to zero) and of the random intercept alone
#   - PR-AUC, and the outcome base rate on the test days (its chance level)
#   - sensitivity, specificity and precision at the threshold
# Each metric's within-imputation variance comes from a participant-level
# bootstrap: participants are resampled with replacement N_BOOT times (a
# participant drawn k times contributes all of their test days k times), the
# metric is recomputed in every resample with the fitted predictions and the
# threshold held fixed, and the variance of the resampled values is taken. This
# accounts for the clustering of test days within participants, which the
# DeLong (AUC) and binomial (proportion) variances for independent observations
# do not. Participant draws are seeded by imputation number, so models that
# share a test set (same horizon) share the same draws.
#
# Metrics are pooled across the 50 imputations by Rubin's rules: total variance =
# mean within-imputation variance + (1 + 1/M) x between-imputation variance, with
# Rubin (1987) degrees of freedom for the t-based 95% CI.
#
# Within-person AUC: for each participant, the AUC over that participant's own
# test days is computed in every imputation in which those days include both
# outcome classes, then averaged across those imputations. The distribution
# across participants is summarised by its median, interquartile range and the
# proportion above .50 (a description of the distribution, so no CI).
#
# Also writes the sample sizes behind each model: participants and person-days
# per imputation, pooled test-set person-days and the test-day base rate.
#
# Usage: Rscript 03_performance_tables.R
#   RESULTS_DIR optional (default "results"); N_BOOT optional (default 2000).
# Outputs, in <RESULTS_DIR>/tables/: pooled_performance.csv,
#   within_person_auc_summary.csv, sample_sizes.csv; and per model
#   <model>_within_person_auc.csv in the model's folder.
# ===================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pROC)
  library(PRROC)
})

results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
n_boot      <- as.integer(Sys.getenv("N_BOOT", unset = "2000"))
seed_base   <- 20261004L
tables_dir  <- file.path(results_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

models <- c("sameday_full", "sameday_wearable", "nextday_full", "nextday_wearable")
model_file <- function(model, suffix) file.path(results_dir, model, paste0(model, "_", suffix))

# ---- Helpers ------------------------------------------------------

# AUC with a fixed orientation (class 0 = controls, class 1 = cases), so that a
# participant whose pain-increase days score below their other days gets an AUC
# below .50 rather than being flipped above it.
roc_fixed <- function(y, p) roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_fixed <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(roc_fixed(y, p)$auc)
}
# The same AUC as a Mann-Whitney statistic (average ranks for ties); used inside
# the bootstrap because it is much faster than building a ROC object.
auc_rank <- function(y, p) {
  n1 <- sum(y == 1L); n0 <- length(y) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(p)[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
pr_auc <- function(y, p) {
  if (sum(y == 1) == 0 || sum(y == 0) == 0) return(NA_real_)
  pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE)$auc.integral
}
threshold_metrics <- function(y, pred_class) {
  tp <- sum(pred_class == 1L & y == 1L); fn <- sum(pred_class == 0L & y == 1L)
  tn <- sum(pred_class == 0L & y == 0L); fp <- sum(pred_class == 1L & y == 0L)
  c(sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
    specificity = if (tn + fp > 0) tn / (tn + fp) else NA_real_,
    precision   = if (tp + fp > 0) tp / (tp + fp) else NA_real_)
}

# Participant-level bootstrap variance of every metric for one imputation.
bootstrap_variances <- function(g, pred_class, n_boot, seed) {
  y <- as.integer(g$Y)
  uids <- sort(unique(g$uid)); n_uid <- length(uids)
  rows_by_uid <- split(seq_len(nrow(g)), factor(g$uid, levels = uids))
  set.seed(seed)
  draws <- matrix(sample.int(n_uid, n_uid * n_boot, replace = TRUE), nrow = n_boot)
  stats <- matrix(NA_real_, n_boot, 7,
                  dimnames = list(NULL, c("auc", "auc_tree_only", "auc_intercept_only", "pr_auc",
                                          "sensitivity", "specificity", "precision")))
  for (b in seq_len(n_boot)) {
    rows <- unlist(rows_by_uid[draws[b, ]], use.names = FALSE)
    yb <- y[rows]
    stats[b, ] <- c(auc_rank(yb, g$full_prob[rows]), auc_rank(yb, g$tree_prob[rows]),
                    auc_rank(yb, g$re_prob[rows]), pr_auc(yb, g$full_prob[rows]),
                    threshold_metrics(yb, pred_class[rows]))
  }
  apply(stats, 2, var, na.rm = TRUE)
}

# Rubin's rules. Q = per-imputation estimates; U = per-imputation within-imputation
# variances.
pool_rubin <- function(Q, U) {
  keep <- !is.na(Q)
  Q <- Q[keep]
  m <- length(Q)
  Q_bar <- mean(Q)
  B     <- var(Q)                                   # between-imputation variance
  U_bar <- mean(U[keep], na.rm = TRUE)              # within-imputation variance
  total_var <- U_bar + (1 + 1 / m) * B
  if (U_bar > 0) {
    r      <- (1 + 1 / m) * B / U_bar
    df     <- (m - 1) * (1 + 1 / r)^2
    t_crit <- qt(0.975, df)
  } else {
    df     <- Inf
    t_crit <- 1.96
  }
  data.frame(estimate = Q_bar,
             ci_lower = Q_bar - t_crit * sqrt(total_var),
             ci_upper = Q_bar + t_crit * sqrt(total_var),
             df = df, n_imputations = m,
             within_imputation_var = U_bar, between_imputation_var = B)
}

# Per-imputation metrics on the test days, with bootstrap variances.
per_imputation_metrics <- function(preds, thresholds) {
  preds %>%
    group_by(imputation) %>%
    group_modify(function(g, key) {
      thr <- thresholds$threshold[thresholds$imputation == key$imputation]
      pred_class <- as.integer(g$full_prob > thr)
      tm <- threshold_metrics(as.integer(g$Y), pred_class)
      v  <- bootstrap_variances(g, pred_class, n_boot, seed_base + key$imputation)
      data.frame(
        auc                = auc_fixed(g$Y, g$full_prob),  auc_var                = v[["auc"]],
        auc_tree_only      = auc_fixed(g$Y, g$tree_prob),  auc_tree_only_var      = v[["auc_tree_only"]],
        auc_intercept_only = auc_fixed(g$Y, g$re_prob),    auc_intercept_only_var = v[["auc_intercept_only"]],
        pr_auc             = pr_auc(g$Y, g$full_prob),     pr_auc_var             = v[["pr_auc"]],
        sensitivity        = tm[["sensitivity"]],          sensitivity_var        = v[["sensitivity"]],
        specificity        = tm[["specificity"]],          specificity_var        = v[["specificity"]],
        precision          = tm[["precision"]],            precision_var          = v[["precision"]])
    }) %>%
    ungroup()
}

fmt <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))

# ---- Per model ----------------------------------------------------
cat(sprintf("Participant-level bootstrap with %d resamples per imputation\n", n_boot))
pooled_rows <- list(); within_rows <- list(); size_rows <- list()
for (m in models) {
  preds      <- read_csv(model_file(m, "test_predictions.csv"), show_col_types = FALSE)
  thresholds <- read_csv(model_file(m, "tuned_hyperparameters.csv"), show_col_types = FALSE)
  sizes      <- read_csv(model_file(m, "analytic_sample_per_imputation.csv"), show_col_types = FALSE)

  # Pooled performance
  pim <- per_imputation_metrics(preds, thresholds)
  pooled <- bind_rows(
    auc                = pool_rubin(pim$auc,                pim$auc_var),
    auc_tree_only      = pool_rubin(pim$auc_tree_only,      pim$auc_tree_only_var),
    auc_intercept_only = pool_rubin(pim$auc_intercept_only, pim$auc_intercept_only_var),
    pr_auc             = pool_rubin(pim$pr_auc,             pim$pr_auc_var),
    sensitivity        = pool_rubin(pim$sensitivity,        pim$sensitivity_var),
    specificity        = pool_rubin(pim$specificity,        pim$specificity_var),
    precision          = pool_rubin(pim$precision,          pim$precision_var),
    .id = "metric")
  base_rate <- data.frame(metric = "base_rate", estimate = mean(preds$Y),
                          ci_lower = NA, ci_upper = NA, df = NA, n_imputations = NA,
                          within_imputation_var = NA, between_imputation_var = NA)
  pooled_rows[[m]] <- bind_rows(pooled, base_rate) %>% mutate(model = m, .before = 1)

  # Within-person AUC
  per_person <- preds %>%
    group_by(imputation, uid) %>%
    summarise(auc = if (sum(Y) > 0 && sum(1 - Y) > 0) auc_fixed(Y, full_prob) else NA_real_,
              .groups = "drop") %>%
    filter(!is.na(auc)) %>%
    group_by(uid) %>%
    summarise(within_person_auc = mean(auc), n_imputations_evaluable = n(), .groups = "drop")
  write_csv(per_person, model_file(m, "within_person_auc.csv"))
  within_rows[[m]] <- data.frame(
    model = m,
    median_auc = median(per_person$within_person_auc),
    iqr_lower  = as.numeric(quantile(per_person$within_person_auc, 0.25)),
    iqr_upper  = as.numeric(quantile(per_person$within_person_auc, 0.75)),
    n_above_50 = sum(per_person$within_person_auc > 0.50),
    n_evaluable = nrow(per_person),
    n_participants = n_distinct(preds$uid))
  within_rows[[m]]$proportion_above_50 <- within_rows[[m]]$n_above_50 / within_rows[[m]]$n_evaluable

  # Sample sizes
  size_rows[[m]] <- data.frame(
    model = m,
    participants_per_imputation_min = min(sizes$n_participants),
    participants_per_imputation_max = max(sizes$n_participants),
    person_days_per_imputation_min  = min(sizes$n_person_days),
    person_days_per_imputation_max  = max(sizes$n_person_days),
    pooled_test_person_days = nrow(distinct(preds, uid, day_number)),
    test_day_base_rate = mean(preds$Y))

  p <- function(metric) { r <- pooled[pooled$metric == metric, ]; sprintf("%s [%s, %s]", fmt(r$estimate), fmt(r$ci_lower), fmt(r$ci_upper)) }
  w <- within_rows[[m]]
  cat(sprintf("\n%s\n  AUC %s | tree-only %s | intercept-only %s\n  PR-AUC %s (base rate %s) | sensitivity %s | specificity %s | precision %s\n  within-person AUC median %s (IQR %s-%s), %d/%d (%.0f%%) above .50\n",
              m, p("auc"), p("auc_tree_only"), p("auc_intercept_only"),
              p("pr_auc"), fmt(mean(preds$Y)), p("sensitivity"), p("specificity"), p("precision"),
              fmt(w$median_auc), fmt(w$iqr_lower), fmt(w$iqr_upper), w$n_above_50, w$n_evaluable,
              100 * w$proportion_above_50))
}

write_csv(bind_rows(pooled_rows), file.path(tables_dir, "pooled_performance.csv"))
write_csv(bind_rows(within_rows), file.path(tables_dir, "within_person_auc_summary.csv"))
write_csv(bind_rows(size_rows),   file.path(tables_dir, "sample_sizes.csv"))
cat(sprintf("\nWrote pooled_performance.csv, within_person_auc_summary.csv and sample_sizes.csv to %s\n", tables_dir))
