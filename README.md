# Cox Forest Plots in R

Publication-quality forest plots for Cox proportional hazards analyses and treatment-effect subgroup analyses.

This repository provides two R functions for generating publication-style forest plots from Cox proportional hazards models:

* **`cox_forest()`** — standard univariable and multivariable Cox regression forest plots
* **`cox_subgroup_forest()`** — treatment-effect subgroup forest plots for two-arm studies

The plots are designed to resemble the style commonly used in clinical trial publications and can be used across one or multiple survival endpoints.

---

## Features

### `cox_forest()`

`cox_forest()` is intended for standard prognostic-factor or baseline-variable Cox regression analyses.

It supports:

* Univariable Cox models
* Multivariable Cox models
* Both model types simultaneously
* Multiple survival endpoints
* Categorical and continuous covariates
* User-defined reference levels
* Optional exclusion of specific covariate levels
* Optional weighting
* Custom endpoint colours
* Custom axis limits and tick marks
* Publication-style forest plots
* Reference rows for categorical variables
* Optional separate covariate header rows

For univariable models, one Cox model is fitted per covariate. For multivariable models, all specified covariates are included in a single Cox model for each endpoint.

### `cox_subgroup_forest()`

`cox_subgroup_forest()` is specifically designed for treatment-effect subgroup analyses in two-arm studies.

For each covariate, the function estimates the treatment hazard ratio separately within each subgroup and additionally reports the overall treatment effect ("All patients").

It supports:

* Two treatment arms
* Multiple subgroup variables
* Categorical subgroups
* Continuous covariates
* Multiple survival endpoints
* User-defined treatment-arm order
* Custom treatment labels
* Kaplan–Meier median survival estimates per treatment arm
* Optional exclusion of specific subgroup levels
* Endpoint-specific colours
* Combined forest plots
* Supplementary effects for continuous covariates

For continuous covariates, the forest plot uses a median split (`< median` / `>= median`). In addition, a separate stratified Cox model provides the hazard ratio per one-unit increase in the continuous covariate.

---

## Requirements

The functions are written in **R** and use the following packages.

### `cox_forest()`

```r
install.packages(c(
  "survival",
  "broom",
  "dplyr",
  "tidyr",
  "purrr",
  "forestploter",
  "tibble",
  "stringr"
))
```

The function additionally uses `grid` from base R.

### `cox_subgroup_forest()`

```r
install.packages(c(
  "survival",
  "broom",
  "dplyr",
  "tidyr",
  "purrr",
  "forestploter",
  "tibble",
  "stringr"
))
```

The source file also lists `survminer`, `ggplot2`, and `gridExtra` among its requirements.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/your-repository.git
cd your-repository
```

Then source the R files:

```r
source("cox_forest.R")
source("cox_subgroup_forest.R")
```

---

# 1. Standard Cox Forest Plot

## Basic usage

```r
result <- cox_forest(
  data = my_trial_data,
  covariates = c("arm", "sex", "age", "ecog"),
  labels = c(
    "Treatment",
    "Sex",
    "Age",
    "ECOG PS"
  ),
  time = c("os_time", "pfs_time"),
  event = c("os_event", "pfs_event"),
  endpoint_labels = c("OS", "PFS"),
  continuous_vars = "age",
  model_type = "both"
)
```

The function returns a list containing the model results and forest plot objects.

### Display the plots

For both univariable and multivariable analyses:

```r
result$plots$univariable
result$plots$multivariable
```

---

## Univariable Cox regression

To generate only univariable models:

```r
result <- cox_forest(
  data = my_trial_data,
  covariates = c("sex", "age", "ecog"),
  labels = c("Sex", "Age", "ECOG PS"),
  time = "os_time",
  event = "os_event",
  continuous_vars = "age",
  model_type = "univariable"
)
```

Each covariate is analysed in a separate Cox model:

```text
Surv(time, event) ~ covariate
```

---

## Multivariable Cox regression

```r
result <- cox_forest(
  data = my_trial_data,
  covariates = c("sex", "age", "ecog"),
  labels = c("Sex", "Age", "ECOG PS"),
  time = "os_time",
  event = "os_event",
  continuous_vars = "age",
  model_type = "multivariable"
)
```

All covariates are included simultaneously:

```text
Surv(time, event) ~ sex + age + ecog
```

The function automatically displays a reference row for categorical variables.

---

## Reference levels

Reference categories can be explicitly defined:

```r
result <- cox_forest(
  data = my_trial_data,
  covariates = c("sex", "ecog"),
  labels = c("Sex", "ECOG PS"),
  time = "os_time",
  event = "os_event",
  ref_levels = list(
    sex = "Female",
    ecog = "0"
  )
)
```

If no reference level is specified, the first factor level is used.

---

## Continuous variables

Continuous variables can be explicitly specified:

```r
continuous_vars = c("age", "bmi")
```

For continuous covariates, the forest plot reports the hazard ratio per specified increase.

For example:

```r
continuous_increase = list(
  age = 10,
  bmi = 5
)
```

This allows the result to be interpreted as the HR per 10-year increase in age or per 5-unit increase in BMI.

---

## Excluding levels

Specific levels can be excluded:

```r
exclude_levels = list(
  ecog = "Unknown",
  sex = "Unknown"
)
```

For `cox_forest()`, excluded values are treated as missing before model fitting. In univariable models this affects the respective covariate model; in multivariable models the affected observations are removed from the complete-case dataset.

---

# 2. Treatment-Effect Subgroup Forest Plot

## Basic usage

```r
result <- cox_subgroup_forest(
  data = my_trial_data,
  covariates = c("sex", "age", "ecog"),
  labels = c(
    "Sex",
    "Age",
    "ECOG PS"
  ),
  treatment = "arm",
  time = "os_time",
  event = "os_event",
  continuous_vars = "age"
)
```

The treatment variable must contain exactly **two treatment levels**.

### Display an individual subgroup plot

```r
result$plots[["Sex"]]
```

### Display the combined plot

```r
result$combined_plots
```

The returned object also contains the underlying Cox models, results table, median survival table and supplementary continuous-variable effects.

---

## Treatment-arm order

The order of the treatment arms determines the direction of the reported hazard ratio.

```r
arm_order = c("Treatment", "Control")
```

The resulting HR is therefore interpreted as:

```text
Treatment vs Control
```

Custom display labels can be supplied independently:

```r
treatment_labels = c(
  "Experimental treatment",
  "Control"
)
```

---

## Multiple endpoints

Multiple endpoints can be analysed simultaneously:

```r
result <- cox_subgroup_forest(
  data = my_trial_data,
  covariates = c("sex", "age"),
  labels = c("Sex", "Age"),
  treatment = "arm",
  time = c("os_time", "pfs_time"),
  event = c("os_event", "pfs_event"),
  endpoint_labels = c("Overall Survival", "Progression-Free Survival")
)
```

The forest plot uses different colours for the different endpoints and provides a corresponding legend.

---

## Continuous subgroup variables

Continuous variables are handled differently from categorical subgroup variables.

For the forest plot, the variable is split at its median:

```text
< median
>= median
```

The treatment effect is then estimated separately within both groups.

In addition, the function calculates a supplementary hazard ratio per one-unit increase using a stratified Cox model:

```text
Surv(time, event) ~ covariate + strata(treatment)
```

This supplementary estimate is available in:

```r
result$continuous_effects
```

The two estimates should not be interpreted as the same quantity: the forest plot represents treatment effects within median-defined subgroups, whereas `continuous_effects` represents the prognostic effect of the continuous covariate.

---

# Output Objects

## `cox_forest()`

The function returns:

```r
result$results_table
```

A tidy long-format table containing:

* covariate
* covariate level
* endpoint
* model type
* N
* number of events
* hazard ratio
* 95% confidence interval
* p-value

The plots are available as:

```r
result$plots$univariable
result$plots$multivariable
```

depending on the selected `model_type`.

---

## `cox_subgroup_forest()`

The returned object contains:

```r
result$models
result$results_table
result$median_survival_table
result$plots
result$combined_plots
result$continuous_effects
result$arm_order
result$arm_labels
```

The `results_table` contains subgroup-level HRs, confidence intervals, p-values, sample sizes, events and median survival by treatment arm.

---

# Important Considerations

### Missing data

Both functions use complete-case data for the variables required by each Cox model.

### Reference categories

For categorical variables, the reference level is explicitly represented in the forest plot for `cox_forest()`.

### Treatment effect interpretation

For `cox_subgroup_forest()`, the treatment HR is defined according to `arm_order`. Make sure the treatment-arm order is chosen intentionally before interpreting HR values.

### Excluded levels

No subgroup level is automatically excluded. If levels such as `"Unknown"` or `"Not reported"` should be excluded, they must be specified using `exclude_levels`.

### Continuous variables

Continuous variables should be explicitly declared using `continuous_vars` when their intended treatment differs from the automatic detection based on numeric variables with more than five distinct values.

---

# Example Workflow

A typical analysis workflow could look like this:

```r
# Load functions
source("cox_forest.R")
source("cox_subgroup_forest.R")

# Standard Cox analysis
cox_result <- cox_forest(
  data = trial_data,
  covariates = c("arm", "sex", "age", "ecog"),
  labels = c("Treatment", "Sex", "Age", "ECOG PS"),
  time = c("os_time", "pfs_time"),
  event = c("os_event", "pfs_event"),
  endpoint_labels = c("OS", "PFS"),
  continuous_vars = "age",
  model_type = "both"
)

# Subgroup treatment-effect analysis
subgroup_result <- cox_subgroup_forest(
  data = trial_data,
  covariates = c("sex", "age", "ecog"),
  labels = c("Sex", "Age", "ECOG PS"),
  treatment = "arm",
  time = c("os_time", "pfs_time"),
  event = c("os_event", "pfs_event"),
  endpoint_labels = c("OS", "PFS"),
  continuous_vars = "age"
)

# Display plots
cox_result$plots$multivariable
subgroup_result$combined_plots
```

---

# Functions at a Glance

| Function                | Purpose                        | Model                                    | Main Use Case                           |
| ----------------------- | ------------------------------ | ---------------------------------------- | --------------------------------------- |
| `cox_forest()`          | Standard Cox forest plot       | Univariable / Multivariable              | Prognostic factors / baseline variables |
| `cox_subgroup_forest()` | Treatment-effect subgroup plot | Treatment-only Cox model within subgroup | RCT subgroup analysis                   |

---



---

# Author

**Julian Silan**

---

# License

Add the appropriate license for your project, for example:

```text
MIT License
```

If this repository is intended for publication or clinical research, please ensure that the statistical methodology, assumptions of the Cox proportional hazards model, reference-level coding and handling of missing data are reviewed before using the generated plots for final reporting.
