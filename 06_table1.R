# ===================================================================
# Table 1 and Supplementary Table 1: descriptive statistics for the analytic
# sample (every participant who entered at least one of the 50 imputed models),
# computed from the observed, pre-imputation data, plus the observed rate of
# pain-increase days after the reference period.
#
# Usage: Rscript 06_table1.R
#   DATA_DIR (default "data") must hold mavdata.csv; RESULTS_DIR (default
#   "results") must hold the same-day full model's test predictions, which
#   identify the analytic sample.
# Output: <RESULTS_DIR>/tables/table1_descriptives.csv (one row per table cell)
# ===================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

data_dir    <- Sys.getenv("DATA_DIR", unset = "data")
results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")
tables_dir  <- file.path(results_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_csv(file.path(data_dir, "mavdata.csv"), show_col_types = FALSE)

# Analytic sample: every participant with test-day predictions in any imputation.
preds <- read_csv(file.path(results_dir, "sameday_full", "sameday_full_test_predictions.csv"),
                  show_col_types = FALSE)
analytic_uids <- sort(unique(preds$uid))

# Physiologically implausible wearable values set to missing (the same rules,
# in the same order, as 01_imputation.R).
dat <- dat %>% mutate(
  `d-sleep-duration_sum_value` = if_else(`d-sleep-duration_sum_value` > 24 & !is.na(`d-sleep-duration_sum_value`), NA_real_, `d-sleep-duration_sum_value`),
  `d-sedentary_minutes_activity_value_sum_value` = if_else(`d-sedentary_minutes_activity_value_sum_value` > 1440 & !is.na(`d-sedentary_minutes_activity_value_sum_value`), NA_real_, `d-sedentary_minutes_activity_value_sum_value`),
  `d-sleep-rem-minutes_sum_value` = if_else((`d-sleep-rem-minutes_sum_value` > `d-sleep-minutes-asleep_sum_value`) | (`d-sleep-rem-minutes_sum_value` > 1440), NA_real_, `d-sleep-rem-minutes_sum_value`, `d-sleep-rem-minutes_sum_value`),
  `d-distance_sum_value` = if_else(`d-distance_sum_value` > 100000 & !is.na(`d-distance_sum_value`), NA_real_, `d-distance_sum_value`),
  `d-sleep-minutes-asleep_sum_value` = if_else(`d-sleep-minutes-asleep_sum_value` > 1440 & !is.na(`d-sleep-minutes-asleep_sum_value`), NA_real_, `d-sleep-minutes-asleep_sum_value`),
  `d-sleep-timein-bed_sum_value` = if_else(`d-sleep-timein-bed_sum_value` > 1440 & !is.na(`d-sleep-timein-bed_sum_value`), NA_real_, `d-sleep-timein-bed_sum_value`))

fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_n_pct   <- function(x) { n <- sum(x == 1, na.rm = TRUE); sprintf("%d (%.1f%%)", n, 100 * n / sum(!is.na(x))) }
fmt_miss    <- function(x) sprintf("%.1f%%", 100 * mean(is.na(x)))

daily_cols <- c("d-steps_sum_value","d-distance_sum_value","d-stairs_sum_value","d-sleep-efficiency_mean_value",
  "d-sleep-minutes-tosleep_mean_value","d-sleep-minutes-awake_sum_value","d-sleep-minutes-afterwake_sum_value",
  "d-sleep-minutes-asleep_sum_value","d-sleep-timein-bed_sum_value","d-sleep-rem-minutes_sum_value",
  "sleephours_daily","sleepquality_daily","stress_daily","dep_daily","anx_daily","pcl4tot_daily","panaspos_daily",
  "panasneg_daily","cigquant_daily_sum","vapequant_daily_sum","canhours_daily_sum","alcdrinks_daily_sum","SSS_daily","pain_daily")
base_cols <- c("age_0","woman_0","minority_0","deploy_0","lectot_0","tbiscreen_0","pain_0","pcl_0","isi_0","phqcase_0",
  "gadcase_0","ucla_0","passtot_0","pctitot_0","uppstot_0","bmacq_0","mmq_0","alcbinge_0","alcalone_0",
  "alcbingealone_0","alcdays_0","mjdays_0","mjalone_0")
binary <- c("woman_0","minority_0","tbiscreen_0","phqcase_0","gadcase_0")

d <- dat %>% filter(uid %in% analytic_uids)
b <- d %>% distinct(uid, .keep_all = TRUE)      # one row per participant for baseline measures

rows <- bind_rows(
  tibble(level = "daily", variable = daily_cols,
         value   = sapply(daily_cols, function(cc) fmt_mean_sd(d[[cc]])),
         missing = sapply(daily_cols, function(cc) fmt_miss(d[[cc]]))),
  tibble(level = "baseline", variable = base_cols,
         value   = sapply(base_cols, function(cc) if (cc %in% binary) fmt_n_pct(b[[cc]]) else fmt_mean_sd(b[[cc]])),
         missing = sapply(base_cols, function(cc) fmt_miss(b[[cc]]))))

# Observed pain-increase rate after the reference period, using the outcome rule
# of the pipeline (daily pain at least 2 points above the rounded reference mean).
ref <- d %>% filter(day_number %in% 1:14) %>% group_by(uid) %>%
  summarise(reference_mean = mean(pain_daily, na.rm = TRUE), .groups = "drop") %>%
  mutate(threshold = round(reference_mean * 2) / 2 + 2)
post <- d %>% filter(day_number > 14, !is.na(pain_daily)) %>%
  inner_join(ref, by = "uid") %>% filter(!is.nan(reference_mean))
meta <- tibble(level = "sample",
               variable = c("n_participants", "n_person_days", "n_post_reference_days_with_observed_pain",
                            "pain_increase_rate_observed"),
               value = c(as.character(length(analytic_uids)), format(nrow(d), big.mark = ","),
                         format(nrow(post), big.mark = ","),
                         sprintf("%.1f%%", 100 * mean(post$pain_daily >= post$threshold))),
               missing = "")

table1 <- bind_rows(rows, meta)
write_csv(table1, file.path(tables_dir, "table1_descriptives.csv"))
print(meta)
cat(sprintf("Wrote table1_descriptives.csv to %s\n", tables_dir))
