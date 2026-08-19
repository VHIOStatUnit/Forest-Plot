# Cox Subgroup Forest Plot in R

This repository contains a publication-style subgroup forest-plot function for
univariable Cox proportional-hazards models, together with a reproducible
example based on the built-in `lung` dataset from the **survival** package.

## Contents

- `forestplot_favouring_function.R` — main function and helper functions.
- `example_lung_cox_subgroup_forest.R` — a short, runnable example using `lung`.
- `tutorial_cox_subgroup_forest.Rmd` — a step-by-step tutorial.
- `README.md` — project overview and usage guide.

## What the function does

`cox_subgroup_forest()` fits a treatment-effect Cox model separately within
each subgroup defined by each requested covariate:

```r
Surv(time, event) ~ treatment
```

For every covariate, the function produces:

- an **All patients** row;
- one row for each subgroup;
- the treatment-effect hazard ratio (HR);
- 95% confidence interval;
- p-value;
- sample size and number of events;
- Kaplan–Meier median survival by treatment arm; and
- a forest plot with a log-scaled HR axis and **Favours [Arm A] /
  Favours [Arm B]** annotations.

The implementation fixes the treatment reference level so that the reported
HR is consistently **Arm A vs Arm B**. An HR below 1 therefore favours Arm A,
assuming the endpoint is coded so that higher hazard means a worse outcome.

Continuous covariates can be supplied explicitly through `continuous_vars`.
They are split at their median for the subgroup rows. The function also
returns a separate `continuous_effects` table containing the HR per one-unit
increase from a stratified Cox model.

## Important note about the `lung` example

The `lung` dataset is an observational lung-cancer dataset; it is **not a
randomised clinical trial and it does not contain a treatment-arm variable**.

To make the forest-plot function runnable without inventing a clinical
treatment variable, the example uses `sex` as the two-level `treatment`
variable and analyses `ecog` and `age_grp` as subgroup variables. This is
strictly a **demonstration of the software and plotting workflow**, not a
causal treatment-effect analysis.

The example first performs the requested data preparation:

```r
df <- lung %>%
  mutate(
    sex        = factor(sex, levels = c(1, 2),
                        labels = c("Male", "Female")),
    ecog       = factor(ph.ecog, levels = c(0, 1, 2, 3),
                        labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3")),
    age_grp    = factor(ifelse(age < 63, "Age < 63", "Age >= 63")),
    status_bin = as.numeric(status == 2)
  ) %>%
  drop_na(time, status_bin, sex, ecog, age_grp, age)
```

The code uses ASCII `>=` in the runnable example so that it is robust across
R environments and output devices.

## Installation

Install the required packages once:

```r
install.packages(c(
  "survival",
  "broom",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "forestploter",
  "tibble",
  "stringr"
))
```

Then source the function:

```r
source("forestplot_favouring_function.R")
```

## Basic usage

A minimal call looks like:

```r
result <- cox_subgroup_forest(
  data = df,
  covariates = c("ecog", "age_grp"),
  labels = c("ECOG performance status", "Age group"),
  treatment = "sex",
  time = "time",
  event = "status_bin",
  endpoint_labels = "Overall survival",
  arm_order = c("Male", "Female"),
  treatment_labels = c("Male", "Female"),
  palette = "grey30",
  time_unit = "months"
)
```

Here, `Male` is treated as Arm A and `Female` as Arm B. This is only a
software demonstration using the `lung` dataset.

## Inspecting the results

The returned object is a list. Useful components include:

```r
result$results_table
result$median_survival_table
result$continuous_effects
result$plots
result$combined_plots
result$models
```

For example:

```r
result$results_table
result$plots[["ECOG performance status"]]
result$combined_plots
```

The `results_table` contains one row per covariate × subgroup × endpoint.

## Model interpretation

Suppose a row reports:

```text
HR = 0.70 (95% CI 0.50-0.98)
```

Because the function reports **Arm A vs Arm B**, an HR below 1 indicates a
lower estimated hazard in Arm A than Arm B. The confidence interval gives the
uncertainty around that estimate.

For the `lung` demonstration, however, this should be interpreted as a
**sex-associated hazard comparison**, not as evidence that a treatment causes
a survival benefit.

## Categorical and continuous covariates

By default, numeric covariates with more than five distinct non-missing values
are treated as continuous. You can override this with:

```r
continuous_vars = "age"
```

A continuous covariate is split at its median for the forest-plot subgroup
rows:

- `< median`
- `>= median`

The function separately estimates a per-unit covariate HR in
`continuous_effects`. These are different estimands and should not be
interpreted as interchangeable.

## Multiple endpoints

The function accepts paired `time` and `event` vectors:

```r
time = c("os_time", "pfs_time"),
event = c("os_event", "pfs_event"),
endpoint_labels = c("Overall survival", "Progression-free survival")
```

The endpoint labels are used to distinguish the estimates in the plot.

## Reproducibility and analysis cautions

This function is designed for subgroup presentation and exploration. In a
real clinical-trial analysis, subgroup HRs should be accompanied by an
appropriate prespecified interaction analysis when the scientific question
is whether treatment effects differ between subgroups.

The `lung` example intentionally does not make a causal claim because the
dataset has no randomised treatment assignment.

## Source implementation notes

The supplied function validates that the treatment variable has exactly two
levels and fixes the reference coding so the displayed HR is consistently
Arm A versus Arm B. It also calculates median survival by arm and returns
the fitted Cox models alongside the tidy results.

The current implementation expects an explicit `palette` whose length
matches the number of endpoint labels. The examples therefore provide a
one-colour palette for the single-endpoint demonstration.

