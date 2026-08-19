# coxforest

Publication-quality forest plots for univariable and/or multivariable Cox
proportional hazards models, across one or more endpoints — in a single
function call.

`cox_forest()` fits the models, builds a tidy long-format results table, and
renders a [forestploter](https://cran.r-project.org/package=forestploter)
figure styled like the baseline-covariate / subgroup forest plots you see in
NEJM, Lancet, and JCO: one row per covariate level, an automatic "Reference"
row for every categorical covariate, N / events / HR (95% CI) / p-value
columns, alternating row shading, and — when more than one endpoint is
supplied — colour-coded confidence intervals with a legend underneath the
plot.

![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

## Features

- **Univariable and/or multivariable models** in one call (`model_type =
  "univariable" | "multivariable" | "both"`).
- **Any number of endpoints** (e.g. OS and PFS) fit from paired `time` /
  `event` variable vectors, with automatic colour coding and a legend once
  there is more than one.
- **Automatic covariate typing**: numeric variables with more than 5 distinct
  values are treated as continuous by default; override with
  `continuous_vars`. Continuous HRs can be reported per any custom increment
  via `continuous_increase`.
- **Reference-level control** per covariate (`ref_levels = list(sex =
  "Female")`), with an explicit "Reference" row inserted automatically.
- **Level exclusion** per covariate (`exclude_levels = list(ecog =
  "Unknown")`) that respects complete-case handling correctly for both
  univariable and multivariable models.
- **Optional bold covariate-header rows** (`covariate_label_row = TRUE`) or
  the default "bold first level" style.
- **Adjustable row height, marker style, CI line width/height, plot width,
  axis limits and tick marks** — all common `forestploter` tuning knobs are
  exposed as arguments.
- Returns both the **long-format results table** (for building your own
  supplementary tables) and the **forestploter plot object(s)** (for further
  `edit_plot()` customisation or saving with `ggplot2::ggsave()` /
  `grid::grid.draw()`).

## Installation

This is currently distributed as a single source file rather than a full
installable package. Clone the repository and source the function directly:

```r
source("R/cox_forest.R")
```

Or, if you prefer to depend on it as a lightweight local package, copy
`cox_forest.R` into an `R/` folder alongside the included `DESCRIPTION` and
install with:

```r
# install.packages("devtools")
devtools::install_local(".")
library(coxforest)
```

### Dependencies

```r
install.packages(c(
  "survival", "broom", "dplyr", "tidyr", "purrr",
  "forestploter", "tibble", "stringr", "gtable"
))
```

## Quick start

```r
library(survival)
library(dplyr)

lung2 <- survival::lung %>%
  mutate(
    event_os = status - 1,
    time_os  = time,
    sex      = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    ph.ecog  = factor(ph.ecog)
  ) %>%
  filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))

result <- cox_forest(
  data       = lung2,
  covariates = c("sex", "age", "ph.ecog"),
  labels     = c("Sex", "Age", "ECOG PS"),
  time       = "time_os",
  event      = "event_os",
  endpoint_labels = "OS",
  continuous_vars = "age",
  model_type = "univariable"
)

result$plots$univariable
result$results_table
```

See [`tutorial.md`](tutorial.md) for a complete, narrated walkthrough
covering multi-endpoint plots, multivariable models, custom reference
levels, level exclusion, and styling options — all built on
`survival::lung`.

## Function reference (summary)

```r
cox_forest(
  data, covariates, labels, time, event,
  endpoint_labels = NULL, continuous_vars = NULL, continuous_increase = NULL,
  ref_levels = NULL, exclude_levels = NULL,
  model_type = c("univariable", "multivariable", "both"),
  digits = 2, palette = NULL, xlim = NULL, ticks_at = NULL,
  length_plot = 40, weights = NULL, ci_pch = 15, ci_lty = 1,
  ci_lwd = 3, ci_Theight = 0.4, sizes = 0.5, risk_labels = NULL,
  title = NULL, row_height = 1, covariate_label_row = FALSE
)
```

| Argument | Description |
|---|---|
| `data` | A data.frame / tibble containing every variable referenced below. |
| `covariates`, `labels` | Variable names to model and their display labels (same length/order). Include a treatment/arm variable here like any other covariate. |
| `time`, `event` | Character vectors of survival time and event-indicator variable names, paired by position — one pair per endpoint. |
| `endpoint_labels` | Display labels per endpoint; defaults to `time`. |
| `continuous_vars` | Force specific covariates to be treated as continuous rather than auto-detected. |
| `continuous_increase` | Named list, e.g. `list(age = 10)`, to report a continuous HR per 10-unit increase instead of per 1 unit. |
| `ref_levels` | Named list of reference levels for categorical covariates. |
| `exclude_levels` | Named list of level(s) to drop per covariate (handled as missing data, so complete-case logic stays correct for multivariable models). |
| `model_type` | `"univariable"`, `"multivariable"`, or `"both"`. |
| `digits` | Decimal places for HR / CI. |
| `palette` | Endpoint colours (only used when there is more than one endpoint). |
| `xlim`, `ticks_at` | Axis limits / tick positions; auto-computed if omitted. |
| `row_height` | Row-height multiplier (`1.3` = 30% taller rows). |
| `covariate_label_row` | `TRUE` inserts a separate bold header row per covariate instead of bolding the first level row. |

Returns a list with `results_table` (long-format tibble) and `plots`
(named list with up to `$univariable` and `$multivariable` entries).

## Testing

Unit tests live in [`test-cox_forest.R`](test-cox_forest.R) and use
`testthat`. Run them with:

```r
testthat::test_file("test-cox_forest.R")
```

The tests cover input validation, categorical and continuous covariate
handling, univariable/multivariable/both model types, custom reference
levels, level exclusion (including its interaction with multivariable
complete-case handling), single- and multi-endpoint plots, and the plot
styling arguments (`row_height`, `covariate_label_row`).

## Author

Julian Silan

## License

MIT — see [`LICENSE`](LICENSE).
