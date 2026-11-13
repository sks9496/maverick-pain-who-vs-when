# ===================================================================
# Multiple imputation of the daily diary, wearable and baseline survey data.
#
# Steps: set physiologically implausible wearable values to missing; drop
# participants and then variables with more than 70% missing values; impute
# 50 complete datasets by multilevel multiple imputation by chained equations
# (mice with miceadds methods, 10 iterations): daily variables by two-level
# predictive mean matching (or two-level logistic regression when binary) with
# the participant as the cluster, baseline survey variables at the participant
# level.
#
# Usage: Rscript 01_imputation.R
#   DATA_DIR (default "data") must hold mavdata.csv; the imputed datasets are
#   written to <DATA_DIR>/imputed_datasets/imputed_dataset_<i>.csv.
# ===================================================================
# ===================================================================
# 0. SETUP
# ===================================================================

# Load necessary packages
library(tidyverse)
library(mice)
library(miceadds)

# Data folder
data_dir <- Sys.getenv("DATA_DIR", unset = "data")

# Load data
daily_with_bl <- read.csv(file.path(data_dir, "mavdata.csv"))

# Set seed for reproducibility
set.seed(123)

# ===================================================================
# 1. CONFIGURATION
# ===================================================================

# --- Imputation Settings ---
M                   <- 50    # Number of imputations
maxit               <- 10    # MICE iterations per imputation

# --- Filtering Thresholds ---
# Drop people with missingness *above* this (e.g., 0.7 = 70%)
person_miss_thresh  <- 0.7
# Drop variables with missingness *above* this (after dropping people)
var_miss_thresh     <- 0.7
# Physiologically implausible wearable values are set to missing before imputation (rules in section 3)


# ===================================================================
# 2. DEFINE VARIABLE LISTS
# ===================================================================

all_vars <- colnames(daily_with_bl)

# Identifiers
id_vars <- c("uid", "day_number")

# Baseline variables (end in _0)
baseline_vars <- all_vars[grepl("_0$", all_vars)]

# Daily variables (start with 'd.' or end in '_daily' or '_sum')
# And are NOT the ID or baseline vars
daily_vars <- all_vars[
  (grepl("^d\\.", all_vars) | grepl("_daily$", all_vars) | grepl("_daily_sum$", all_vars)) &
    !(all_vars %in% c(id_vars, baseline_vars))
]

# All variables to be used in imputation
vars_to_impute <- c(daily_vars, baseline_vars)


# ===================================================================
# 3. PRE-IMPUTATION FILTERING
# ===================================================================


# --- Clean Impossible Sensor Values ---

daily_cleaned <- daily_with_bl %>%
  mutate(
    # Rule 1: Sleep duration (Time in Bed) > 24 hours is impossible
    d.sleep.duration_sum_value = if_else(
      d.sleep.duration_sum_value > 24 & !is.na(d.sleep.duration_sum_value),
      NA_real_,
      d.sleep.duration_sum_value
    ),
    
    # Rule 2: Sedentary minutes > 1440 (total mins in a day) is impossible
    d.sedentary_minutes_activity_value_sum_value = if_else(
      d.sedentary_minutes_activity_value_sum_value > 1440 & !is.na(d.sedentary_minutes_activity_value_sum_value),
      NA_real_,
      d.sedentary_minutes_activity_value_sum_value
    ),
    
    # Rule 3: REM sleep > total minutes asleep OR > 1440 (hard cap)
    d.sleep.rem.minutes_sum_value = if_else(
      (d.sleep.rem.minutes_sum_value > d.sleep.minutes.asleep_sum_value) | (d.sleep.rem.minutes_sum_value > 1440),
      NA_real_,                       # Value if TRUE (it's bad)
      d.sleep.rem.minutes_sum_value,  # Value if FALSE (it's good)
      d.sleep.rem.minutes_sum_value   # Value if NA (it's unverifiable, so we keep it)
    ),
    
    # Rule 4: Distance > 100,000 meters (100km/62mi) is impossible
    d.distance_sum_value = if_else(
      d.distance_sum_value > 100000 & !is.na(d.distance_sum_value),
      NA_real_,
      d.distance_sum_value
    ),
    
    # Rule 5: Minutes Asleep > 1440 (total mins in a day) is impossible
    d.sleep.minutes.asleep_sum_value = if_else(
      d.sleep.minutes.asleep_sum_value > 1440 & !is.na(d.sleep.minutes.asleep_sum_value),
      NA_real_,
      d.sleep.minutes.asleep_sum_value
    ),
    
    # Rule 6: Time in Bed > 1440 (total mins in a day) is impossible
    d.sleep.timein.bed_sum_value = if_else(
      d.sleep.timein.bed_sum_value > 1440 & !is.na(d.sleep.timein.bed_sum_value),
      NA_real_,
      d.sleep.timein.bed_sum_value
    )
  )


# --- Filter People (e.g., > 70% missing) ---

initial_n_people <- n_distinct(daily_cleaned$uid)

# Calculate missing proportion for each person
missing_per_person <- daily_cleaned %>%
  select(all_of(vars_to_impute)) %>% # Only check missingness on data columns
  is.na() %>%
  as.data.frame() %>%
  mutate(uid = daily_cleaned$uid) %>% # Add uid back for grouping
  group_by(uid) %>%
  summarise(person_miss_prop = mean(c_across(), na.rm = TRUE), .groups = 'drop')

# Get list of UIDs to *keep*
uids_to_keep <- missing_per_person %>%
  filter(person_miss_prop <= person_miss_thresh) %>%
  pull(uid)

# Filter dataframe
data_filtered <- daily_cleaned %>%
  filter(uid %in% uids_to_keep)

final_n_people <- n_distinct(data_filtered$uid)
cat(sprintf("Dropped %d people. Remaining: %d\n", initial_n_people - final_n_people, final_n_people))


# --- Filter Features (e.g., > 70% missing) ---

initial_n_vars <- length(vars_to_impute)

# Calculate missing proportion for each variable on the *filtered* dataset
missing_per_var_summary <- data_filtered %>%
  select(all_of(vars_to_impute)) %>%
  summarise(across(everything(), ~mean(is.na(.), na.rm = TRUE)))

# Convert the 1-row summary dataframe to the named vector
missing_per_var <- unlist(missing_per_var_summary[1, ])

# Get list of variables to *drop*
vars_to_drop <- names(missing_per_var[missing_per_var > var_miss_thresh])

# Clean potential NAs from the drop list (just in case)
vars_to_drop <- vars_to_drop[!is.na(vars_to_drop)]

if (length(vars_to_drop) > 0) {
  cat("Dropping the following variables:\n")
  print(vars_to_drop)
  
  # Filter the dataset
  data_filtered <- data_filtered %>%
    select(-all_of(vars_to_drop))
  
  # Update the variable lists
  daily_vars <- daily_vars[!(daily_vars %in% vars_to_drop)]
  baseline_vars <- baseline_vars[!(baseline_vars %in% vars_to_drop)]
} else {
  cat("No variables exceeded the missingness threshold.\n")
}

final_n_vars <- length(daily_vars) + length(baseline_vars)
cat(sprintf("Dropped %d variables. Remaining: %d\n", initial_n_vars - final_n_vars, final_n_vars))


# ===================================================================
# 4. ID & NAME MAPPING (Pre-MICE Setup)
# ===================================================================

# --- UID Mapping (Changes UID to simple number - Required for multi-level MICE methods) ---
uid_map <- data.frame(
  uid = unique(data_filtered$uid)
) %>%
  mutate(uid_int = row_number())

# Add the new integer ID to the dataframe
data_to_impute <- data_filtered %>%
  left_join(uid_map, by = "uid")

# --- Variable Name Mapping (Shortens variable names - Recommended for MICE stability) ---
all_impute_vars <- c("uid_int", "day_number", daily_vars, baseline_vars)
variable_map <- data.frame(
  long_name = all_impute_vars
) %>%
  mutate(
    short_name = case_when(
      long_name == "uid_int" ~ "uid_int",
      long_name == "day_number" ~ "day_num",
      TRUE ~ paste0("v", row_number() - 2) # v1, v2, v3...
    )
  )

# Create a "short-name" version of the dataframe for mice
data_mice_ready <- data_to_impute %>%
  select(all_of(c("uid", all_impute_vars))) # Keep original uid for later

# Rename columns to short names
data_mice_ready <- data_mice_ready %>%
  rename_with(
    ~ variable_map$short_name[match(.x, variable_map$long_name)],
    .cols = all_of(variable_map$long_name)
  )


# ===================================================================
# 5. MICE IMPUTATION
# ===================================================================

# --- Helper function to properly set imputation method for each var ---
get_imputation_method <- function(var_name, data, daily_vars_short, baseline_vars_short) {
  
  # Get unique non-NA values
  unique_vals <- unique(na.omit(data[[var_name]]))
  n_unique <- length(unique_vals)
  
  is_binary <- (n_unique == 2) && (all(unique_vals %in% c(0, 1)))
  
  if (var_name %in% daily_vars_short) {
    # It's a daily (Level 1) variable
    if (is_binary) {
      return("2l.bin") # Multilevel method for L1 binary 
    } else {
      return("2l.pmm") # Multilevel method for L1 continuous
    }
    
  } else if (var_name %in% baseline_vars_short) {
    # It's a baseline (Level 2) variable
    return("2lonly.pmm") # Multilevel method for L2 variables
    
  } else {
    # It's an ID variable or something else
    return("") # Don't impute
  }
}


# --- Initialize mice objects ---

# Get short names for method helper
daily_vars_short <- variable_map$short_name[variable_map$long_name %in% daily_vars]
baseline_vars_short <- variable_map$short_name[variable_map$long_name %in% baseline_vars]

# Get a "dry run" mice object to generate default settings
ini <- mice(data_mice_ready %>% select(-uid), m = 1, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# --- Specify Cluster Variable & Methods ---

# Set the cluster variable (uid_int)
pred[, "uid_int"] <- -2 # -2 = Cluster variable
meth["uid_int"] <- ""   # Don't impute the ID itself

# Don't use day_number to predict (or impute it)
pred[, "day_num"] <- 0
meth["day_num"] <- ""

# Set imputation methods for other variables correctly based on helper function 
for (var_short in names(meth)) {
  if (meth[var_short] != "") { # Skip ones we already set
    meth[var_short] <- get_imputation_method(
      var_short,
      data_mice_ready,
      daily_vars_short,
      baseline_vars_short
    )
  }
}

cat("Methods set. Final methods list:\n")
print(table(meth[meth != ""]))

# Convert all columns assigned a binary method to 'factor' (necessary for MICE)
binary_vars_short <- names(meth[meth %in% c("logreg", "2l.bin")])
data_mice_ready <- data_mice_ready %>%
  mutate(across(all_of(binary_vars_short), as.factor))


# --- 5.4: Run Imputation ---

# Run the full imputation
# We select(-uid) because the original character ID is not needed by mice
imputed_data <- mice(
  data_mice_ready %>% select(-uid),
  m = M,
  maxit = maxit,
  predictorMatrix = pred,
  method = meth
)


# ===================================================================
# 6. EXTRACT & SAVE IMPUTATED DATASETS
# ===================================================================

# --- Process and Store Imputations in a List ---

# Create a list to hold the dataframes in R
imputed_datasets_list <- list()

# Create a reverse map (To convert short name back to long name)
reverse_variable_map <- setNames(variable_map$long_name, variable_map$short_name)

# Dataset processing
for (i in 1:M) {
  cat(sprintf("Processing dataset %d/%d...\n", i, M))
  
  # 1. Complete the dataset
  temp_df_short <- mice::complete(imputed_data, i)
  
  # 2. Add back the original character 'uid'
  # We join using the integer ID, then select the original char uid
  temp_df_with_uid <- temp_df_short %>%
    left_join(uid_map, by = "uid_int") %>%
    select(-uid_int) # Drop the integer ID
  
  # 3. Rename columns back to long names for clarity
  temp_df_long <- temp_df_with_uid %>%
    rename_with(
      ~ reverse_variable_map[.x],
      # any_of(): renames the columns present in the map and ignores the rest
      .cols = any_of(names(reverse_variable_map)) 
    )
  
  # 4. Re-order columns to be user-friendly (IDs first)
  final_df <- temp_df_long %>%
    select(uid, day_number, everything())
  
  # 5. Add to the list
  imputed_datasets_list[[paste0("imputation_", i)]] <- final_df
  
}
cat("All datasets processed and stored in 'imputed_datasets_list'.\n")


# --- Save the datasets as CSV files ---

# Create the output directory if it doesn't exist
output_dir <- file.path(data_dir, "imputed_datasets")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through the list just created

for (i in 1:M) {
  # Get the correct list name and file name
  list_name <- paste0("imputation_", i)
  file_name <- file.path(output_dir, paste0("imputed_dataset_", i, ".csv"))
  
  # Write the dataframe from the list
  write.csv(imputed_datasets_list[[list_name]], file_name, row.names = FALSE)
}

cat("All datasets saved to CSVs in the '", output_dir, "' folder.\n")

