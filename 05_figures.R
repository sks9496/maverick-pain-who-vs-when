# ===================================================================
# Figure 1 and Supplementary Figure 1: SHAP beeswarm plots of the retained
# features for the same-day and next-day models (full models in Figure 1,
# wearable-only models in Supplementary Figure 1).
# No model is refit; reads the saved per-imputation SHAP values.
#
# For each model, SHAP values and feature values are averaged across
# imputations within person-day (each point is one pooled test-set person-day).
# Features are ordered top to bottom by mean absolute pooled SHAP value. Point
# colour is the feature value, rescaled between its 5th and 95th percentiles
# (wider percentiles when those coincide); SHAP values beyond the 99.5th
# percentile of their absolute value are drawn at that limit.
#
# Usage: Rscript 05_figures.R      (RESULTS_DIR optional, default "results")
# Output: <RESULTS_DIR>/figures/Figure_1.pdf/.png, Supp_Figure_1.pdf/.png
# ===================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(patchwork)
})

results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
model_file <- function(model, suffix) file.path(results_dir, model, paste0(model, "_", suffix))

# ---- Feature labels -----------------------------------------------
# One name per daily variable, so daily and reference-period labels agree.
var_names <- c(
  dep_daily = "Depression", stress_daily = "Stress", anx_daily = "Anxiety", pcl4tot_daily = "PTSD Symptom",
  panaspos_daily = "Positive Affect", panasneg_daily = "Negative Affect", SSS_daily = "Social Interaction",
  sleephours_daily = "Sleep Hours", sleepquality_daily = "Sleep Quality", canhours_daily_sum = "Cannabis Hours",
  vapequant_daily_sum = "Vaping Quantity", alcdrinks_daily_sum = "Alcohol Drinks", cigquant_daily_sum = "Cigarettes Smoked",
  pain_daily = "Pain",
  d.steps_sum_value = "Step Count", d.distance_sum_value = "Walking Distance", d.stairs_sum_value = "Stairs Climbed",
  d.sleep.minutes.asleep_sum_value = "Sleep Duration", d.sleep.minutes.awake_sum_value = "Wake After Sleep Onset",
  d.sleep.minutes.afterwake_sum_value = "Post-Wake Time in Bed", d.sleep.minutes.tosleep_mean_value = "Sleep Onset Latency",
  d.sleep.timein.bed_sum_value = "Time in Bed", d.sleep.efficiency_mean_value = "Sleep Efficiency",
  d.sleep.rem.minutes_sum_value = "REM Minutes")

# Feature name -> plot label. Deviation features are "Daily ..." in the same-day
# models and "Prior-Day ..." in the next-day models.
label_var <- function(x, prior_day = FALSE) {
  vapply(x, function(v) {
    m <- regmatches(v, regexec("^(.*)_(absdelta|baseline_mean|baseline_sd)$", v))[[1]]
    if (length(m) == 3 && m[2] %in% names(var_names)) {
      src <- if (startsWith(m[2], "d.")) "Wearable" else "Diary"
      nm  <- var_names[[m[2]]]
      switch(m[3],
        absdelta      = sprintf("%s %s Deviation (%s)", if (prior_day) "Prior-Day" else "Daily", nm, src),
        baseline_mean = sprintf("Reference-Period %s Mean (%s)", nm, src),
        baseline_sd   = sprintf("Reference-Period %s SD (%s)", nm, src))
    } else {
      paste("Baseline Survey:", gsub("_0$", "", v))   # baseline survey characteristics
    }
  }, character(1), USE.NAMES = FALSE)
}

# ---- Pooled SHAP and feature values for the retained features ------------
load_pooled <- function(model) {
  shap      <- readRDS(model_file(model, "shap_per_imputation.rds"))
  selection <- read_csv(model_file(model, "variable_selection.csv"), show_col_types = FALSE)
  retained  <- selection$variable[selection$selected]
  si <- which(!sapply(shap$shap_list, is.null))
  common <- Reduce(intersect, lapply(shap$shap_list[si], colnames))
  retained <- intersect(retained, common)
  pool_mat <- function(lst) {
    bind_rows(lapply(si, function(i)
        as.data.frame(lst[[i]][, retained, drop = FALSE]) %>%
          mutate(uid = shap$key_list[[i]]$uid, day_number = shap$key_list[[i]]$day_number))) %>%
      group_by(uid, day_number) %>%
      summarise(across(all_of(retained), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
      select(-uid, -day_number) %>% as.matrix()
  }
  list(shap_mat = pool_mat(shap$shap_list), feat_mat = pool_mat(shap$feat_list), retained = retained)
}

# ---- Beeswarm panel -------------------------------------------------------
make_beeswarm <- function(shap_mat, feat_mat, panel_title, prior_day = FALSE) {
  mean_abs <- colMeans(abs(shap_mat)); feat_order <- names(sort(mean_abs))
  df <- bind_rows(lapply(feat_order, function(v) {
    fv <- feat_mat[, v]
    vmin <- quantile(fv, 0.05, na.rm = TRUE); vmax <- quantile(fv, 0.95, na.rm = TRUE)
    if (vmin == vmax) { vmin <- quantile(fv, 0.01, na.rm = TRUE); vmax <- quantile(fv, 0.99, na.rm = TRUE) }
    if (vmin == vmax) { vmin <- min(fv, na.rm = TRUE); vmax <- max(fv, na.rm = TRUE) }
    fvc <- pmin(pmax(fv, vmin), vmax)
    data.frame(feature_label = label_var(v, prior_day), shap_value = shap_mat[, v],
               feat_scaled = (fvc - vmin) / (vmax - vmin + 1e-10))
  })) %>% mutate(feature_label = factor(feature_label, levels = label_var(feat_order, prior_day)))
  q <- quantile(abs(df$shap_value), 0.995, na.rm = TRUE)
  df$shap_plot <- pmax(pmin(df$shap_value, q), -q)
  ggplot(df, aes(x = shap_plot, y = feature_label, color = feat_scaled)) +
    geom_quasirandom(size = 0.6, alpha = 0.65, width = 0.38, bandwidth = 0.8, groupOnX = FALSE) +
    scale_color_gradientn(colors = c("#1648A0","#6AACD6","#F7F7F7","#F4622A","#A50026"),
      values = c(0,0.25,0.5,0.75,1), name = "Feature\nValue",
      breaks = c(0.02,0.98), labels = c("Low","High"),
      guide = guide_colorbar(barheight = 8, barwidth = 0.9, ticks = FALSE, title.vjust = 1)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
    labs(title = panel_title, x = "SHAP value", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0),
          axis.text.y = element_text(size = 10, color = "black"), axis.text.x = element_text(size = 9),
          axis.title.x = element_text(size = 10, margin = ggplot2::margin(t = 6)),
          axis.line = element_line(color = "grey50", linewidth = 0.35),
          axis.ticks = element_line(color = "grey50", linewidth = 0.3), axis.ticks.y = element_blank(),
          legend.title = element_text(size = 10), legend.text = element_text(size = 9),
          legend.position = "right",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = ggplot2::margin(10, 12, 8, 6))
}

# PDF at the given size in inches, and a 300-dpi PNG at the given pixel height.
save_fig <- function(plot, name, w, h, png_px_h) {
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot, width = w, height = h, device = pdf)
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot, width = w, height = png_px_h / 300, dpi = 300, bg = "white")
  cat(sprintf("  wrote %s.pdf and %s.png\n", name, name))
}

cat("Loading pooled SHAP values...\n")
p_sf <- load_pooled("sameday_full");     p_nf <- load_pooled("nextday_full")
p_sw <- load_pooled("sameday_wearable"); p_nw <- load_pooled("nextday_wearable")

fig1 <- make_beeswarm(p_sf$shap_mat, p_sf$feat_mat, "(a) Same-day full model", FALSE) |
        make_beeswarm(p_nf$shap_mat, p_nf$feat_mat, "(b) Next-day full model", TRUE)
save_fig(fig1, "Figure_1", 18, max(6, max(ncol(p_sf$shap_mat), ncol(p_nf$shap_mat)) * 0.42), 3021)

sfig1 <- make_beeswarm(p_sw$shap_mat, p_sw$feat_mat, "(a) Same-day wearable-only model", FALSE) |
         make_beeswarm(p_nw$shap_mat, p_nw$feat_mat, "(b) Next-day wearable-only model", TRUE)
save_fig(sfig1, "Supp_Figure_1", 18, max(6, max(ncol(p_sw$shap_mat), ncol(p_nw$shap_mat)) * 0.42), 1800)
