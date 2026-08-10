# ==============================================================================
# Example: cox_subgroup_forest() using the survival::lung dataset
# ==============================================================================

# This example demonstrates the function using the built-in lung dataset.
#
# IMPORTANT:
# The lung dataset does not contain a randomised treatment variable.
# For demonstration purposes, sex is used as the two-level comparison
# variable. The resulting HRs are observational Male-vs-Female comparisons,
# NOT causal treatment effects.

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(survival)
library(dplyr)
library(tidyr)

# The function file should be in the same working directory as this script.
source("forestplot_favouring_function.R")


# ------------------------------------------------------------------------------
# 2. Prepare the lung dataset
# ------------------------------------------------------------------------------

data("lung", package = "survival")

df <- lung %>%
  mutate(
    sex        = factor(
      sex,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    ecog       = factor(
      ph.ecog,
      levels = c(0, 1, 2, 3),
      labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3")
    ),
    age_grp    = factor(
      ifelse(age < 63, "Age < 63", "Age >= 63")
    ),
    status_bin = as.numeric(status == 2)
  ) %>%
  drop_na(time, status_bin, sex, ecog, age_grp, age)


# ------------------------------------------------------------------------------
# 3. Quick checks
# ------------------------------------------------------------------------------

table(df$sex)
table(df$ecog)
table(df$age_grp)

summary(df$time)


# ------------------------------------------------------------------------------
# 4. Run the subgroup analysis
# ------------------------------------------------------------------------------

# The function requires exactly two treatment levels.
#
# Here:
#   Arm A = Male
#   Arm B = Female
#
# Consequently, HR = Male / Female and HR < 1 favours Male.

result <- cox_subgroup_forest(
  data = df,
  covariates = c("ecog", "age_grp"),
  labels = c(
    "ECOG performance status",
    "Age group"
  ),
  treatment = "sex",
  time = "time",
  event = "status_bin",
  endpoint_labels = "Overall survival",
  arm_order = c("Male", "Female"),
  treatment_labels = c("Male", "Female"),
  palette = "grey30",
  time_unit = "months"
)


# ------------------------------------------------------------------------------
# 5. Inspect the results
# ------------------------------------------------------------------------------

print(result$results_table)

print(result$median_survival_table)

if (!is.null(result$continuous_effects)) {
  print(result$continuous_effects)
}


# ------------------------------------------------------------------------------
# 6. Display individual forest plots
# ------------------------------------------------------------------------------

print(result$plots[["ECOG performance status"]])

print(result$plots[["Age group"]])


# ------------------------------------------------------------------------------
# 7. Display the combined forest plot
# ------------------------------------------------------------------------------

print(result$combined_plots)


# ------------------------------------------------------------------------------
# 8. Inspect one fitted Cox model
# ------------------------------------------------------------------------------

# View the available model names first:
print(names(result$models[["ECOG performance status"]]))

# Then select a model by one of the names printed above, for example:
#
# fit <- result$models[["ECOG performance status"]][[
#   "ECOG performance status | ECOG 2 | Overall survival"
# ]]
#
# summary(fit)


# ------------------------------------------------------------------------------
# 9. Optional: analyse age as a continuous covariate
# ------------------------------------------------------------------------------

result_age <- cox_subgroup_forest(
  data = df,
  covariates = "age",
  labels = "Age",
  treatment = "sex",
  time = "time",
  event = "status_bin",
  endpoint_labels = "Overall survival",
  continuous_vars = "age",
  arm_order = c("Male", "Female"),
  treatment_labels = c("Male", "Female"),
  palette = "grey30",
  time_unit = "months"
)

print(result_age$plots[["Age"]])

# Per one-unit increase in age:
print(result_age$continuous_effects)


# ------------------------------------------------------------------------------
# 10. Save the tidy results
# ------------------------------------------------------------------------------

# Uncomment if desired:
#
# write.csv(
#   result$results_table,
#   "lung_subgroup_cox_results.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   result$median_survival_table,
#   "lung_median_survival.csv",
#   row.names = FALSE
# )


# ==============================================================================
# End of example
# ==============================================================================
