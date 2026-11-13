# ===================================================================
# Direction of each retained feature's contribution (Supplementary Tables 2 and 3).
# No model is refit; reads the saved per-imputation SHAP values.
#
# For each model and each feature retained by the scree elbow: the Spearman
# correlation between feature value and SHAP value over the test days, computed
# within each of the 50 imputations and pooled by Fisher's z transformation
# (the mean z back-transformed, with a 95% CI from the between-imputation
# standard error of z). A positive correlation means that higher feature values
# shifted predictions toward a pain increase. Features are listed in descending
# order of pooled mean absolute SHAP value, the order used in Figure 1 and
# Supplementary Figure 1.
#
# Usage: Rscript 04_shap_direction.R      (RESULTS_DIR optional, default "results")
# Output: <model>_shap_direction.csv in each model's folder, and
#   <RESULTS_DIR>/tables/shap_direction_all_models.csv
# ===================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
tables_dir  <- file.path(results_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

models <- c("sameday_full", "sameday_wearable", "nextday_full", "nextday_wearable")
model_file <- function(model, suffix) file.path(results_dir, model, paste0(model, "_", suffix))

# Spearman correlation between a feature's values and its SHAP values in one
# imputation; NA when there are too few test days or too few distinct values.
spearman_one_imputation <- function(x, s) {
  ok <- !is.na(x) & !is.na(s)
  if (sum(ok) < 10 || length(unique(x[ok])) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], s[ok], method = "spearman"))
}

# Pool per-imputation correlations through Fisher's z.
pool_fisher <- function(r) {
  r <- r[!is.na(r) & abs(r) < 1]
  z <- atanh(r)
  z_mean <- mean(z)
  z_se   <- sd(z) / sqrt(length(z))
  data.frame(rho = tanh(z_mean), ci_lower = tanh(z_mean - 1.96 * z_se),
             ci_upper = tanh(z_mean + 1.96 * z_se), n_imputations = length(z))
}

all_rows <- list()
for (m in models) {
  shap       <- readRDS(model_file(m, "shap_per_imputation.rds"))
  selection  <- read_csv(model_file(m, "variable_selection.csv"), show_col_types = FALSE)
  importance <- read_csv(model_file(m, "shap_importance.csv"), show_col_types = FALSE)
  retained   <- importance %>%
    filter(variable %in% selection$variable[selection$selected]) %>%
    arrange(desc(mean_abs_shap))

  direction <- bind_rows(lapply(retained$variable, function(v) {
    r <- vapply(seq_along(shap$shap_list), function(i) {
      spearman_one_imputation(shap$feat_list[[i]][, v], shap$shap_list[[i]][, v])
    }, numeric(1))
    data.frame(variable = v, mean_abs_shap = retained$mean_abs_shap[retained$variable == v],
               pool_fisher(r))
  }))
  write_csv(direction, model_file(m, "shap_direction.csv"))
  all_rows[[m]] <- direction %>% mutate(model = m, .before = 1)

  cat(sprintf("\n%s (%d retained features)\n", m, nrow(direction)))
  for (i in seq_len(nrow(direction))) {
    cat(sprintf("  %-50s rho %+.2f [%+.2f, %+.2f]\n", direction$variable[i], direction$rho[i],
                direction$ci_lower[i], direction$ci_upper[i]))
  }
}
write_csv(bind_rows(all_rows), file.path(tables_dir, "shap_direction_all_models.csv"))
cat(sprintf("\nWrote shap_direction_all_models.csv to %s\n", tables_dir))
