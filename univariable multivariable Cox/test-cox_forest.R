# ==============================================================================
# test-cox_forest.R
#
# testthat unit tests for cox_forest().
#
# Run with:
#   testthat::test_file("test-cox_forest.R")
#
# or, if this file lives in tests/testthat/ inside an installed/loaded
# "coxforest" package:
#   devtools::test()
#
# Requires: testthat (>= 3.0.0), survival, dplyr, and every package listed
# as an Import in DESCRIPTION.
# ==============================================================================
rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(testthat)
  library(survival)
  library(dplyr)
})

# ------------------------------------------------------------------------
# Locate and source cox_forest.R regardless of whether this file is run
# standalone, via testthat::test_file(), or from within a package's
# tests/testthat directory.
# ------------------------------------------------------------------------
source("/Users/julian/Desktop/VHIO/R code/Forestplot/Univariable_Multivariable/cox.forest_function.R")

# ==============================================================================
# Shared test data: survival::lung, prepared as in the tutorial
# ==============================================================================

set.seed(42)

lung2 <- survival::lung %>%
  dplyr::mutate(
    event_os  = status - 1,
    time_os   = time,
    sex       = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    ph.ecog   = factor(ph.ecog),
    treatment = sample(c("Arm A", "Arm B"), dplyr::n(), replace = TRUE),
    time_pfs  = pmin(time_os, time_os * stats::runif(dplyr::n(), 0.4, 0.9)),
    event_pfs = ifelse(time_pfs < time_os, 1, event_os)
  ) %>%
  dplyr::filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))

# ==============================================================================
# Input validation
# ==============================================================================

test_that("validates that `data` is a data.frame", {
  expect_error(
    cox_forest(
      data = list(a = 1), covariates = "age", labels = "Age",
      time = "time_os", event = "event_os"
    ),
    "must be a data.frame"
  )
})

test_that("validates that `covariates` is non-empty", {
  expect_error(
    cox_forest(
      data = lung2, covariates = character(0), labels = character(0),
      time = "time_os", event = "event_os"
    ),
    "at least one variable"
  )
})

test_that("validates matching lengths of `covariates` and `labels`", {
  expect_error(
    cox_forest(
      data = lung2, covariates = c("age", "sex"), labels = "Age",
      time = "time_os", event = "event_os"
    ),
    "same length"
  )
})

test_that("validates matching lengths of `time` and `event`", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "age", labels = "Age",
      time = c("time_os", "time_pfs"), event = "event_os"
    ),
    "same length"
  )
})

test_that("validates that referenced variables exist in `data`", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "not_a_column", labels = "Nope",
      time = "time_os", event = "event_os"
    ),
    "not found in `data`"
  )
})

test_that("validates that `endpoint_labels` matches `time` in length", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "age", labels = "Age",
      time = c("time_os", "time_pfs"), event = c("event_os", "event_pfs"),
      endpoint_labels = "OS only"
    ),
    "same length as .time."
  )
})

test_that("validates that `palette` matches `endpoint_labels` in length", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "age", labels = "Age",
      time = c("time_os", "time_pfs"), event = c("event_os", "event_pfs"),
      endpoint_labels = c("OS", "PFS"), palette = "#000000"
    ),
    "same length as .endpoint_labels."
  )
})

test_that("`row_height` must be a single positive number", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "age", labels = "Age",
      time = "time_os", event = "event_os", row_height = -1
    ),
    "positive numeric value"
  )
  expect_error(
    cox_forest(
      data = lung2, covariates = "age", labels = "Age",
      time = "time_os", event = "event_os", row_height = c(1, 2)
    ),
    "positive numeric value"
  )
})

# ==============================================================================
# Univariable models: continuous covariate
# ==============================================================================

test_that("a continuous covariate produces a single results row with the expected columns", {
  res <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os", endpoint_labels = "OS",
    continuous_vars = "age", model_type = "univariable"
  )

  expect_s3_class(res$results_table, "data.frame")
  expect_equal(nrow(res$results_table), 1)
  expect_equal(res$results_table$covariate, "age")
  expect_equal(res$results_table$covariate_label, "Age")
  expect_equal(res$results_table$level, "per 1-unit increase")
  expect_false(res$results_table$is_reference)
  expect_true(is.finite(res$results_table$hr))
  expect_true(res$results_table$hr > 0)
  expect_true(res$results_table$conf_low <= res$results_table$hr)
  expect_true(res$results_table$conf_high >= res$results_table$hr)
  expect_true(is.finite(res$results_table$p_value))

  # Cross-check against a manually fit model
  fit <- survival::coxph(survival::Surv(time_os, event_os) ~ age, data = lung2)
  expect_equal(unname(exp(stats::coef(fit))), res$results_table$hr, tolerance = 1e-8)
})

test_that("`continuous_increase` rescales the HR to the requested increment", {
  res1 <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os",
    continuous_vars = "age", model_type = "univariable"
  )
  res10 <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os",
    continuous_vars = "age", continuous_increase = list(age = 10),
    model_type = "univariable"
  )

  expect_equal(res10$results_table$level, "per 10-unit increase")
  expect_equal(res10$results_table$hr, res1$results_table$hr^10, tolerance = 1e-8)
})

test_that("numeric covariates with > 5 distinct values are auto-detected as continuous", {
  res <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os", model_type = "univariable"
  )
  expect_equal(res$results_table$level, "per 1-unit increase")
})

# ==============================================================================
# Univariable models: categorical covariate
# ==============================================================================

test_that("a categorical covariate produces one row per level, including a Reference row", {
  res <- cox_forest(
    data = lung2, covariates = "sex", labels = "Sex",
    time = "time_os", event = "event_os", model_type = "univariable"
  )

  expect_equal(nrow(res$results_table), nlevels(lung2$sex))
  expect_equal(sort(res$results_table$level), sort(levels(lung2$sex)))

  ref_row <- res$results_table[res$results_table$is_reference, ]
  expect_equal(nrow(ref_row), 1)
  expect_equal(ref_row$level, "Male") # first factor level by default
  expect_true(is.na(ref_row$hr))
  expect_true(is.na(ref_row$conf_low))
  expect_true(is.na(ref_row$conf_high))

  non_ref_row <- res$results_table[!res$results_table$is_reference, ]
  expect_equal(non_ref_row$level, "Female")
  expect_true(is.finite(non_ref_row$hr))
})

test_that("`ref_levels` overrides the default (first-level) reference", {
  res <- cox_forest(
    data = lung2, covariates = "sex", labels = "Sex",
    time = "time_os", event = "event_os", model_type = "univariable",
    ref_levels = list(sex = "Female")
  )
  ref_row <- res$results_table[res$results_table$is_reference, ]
  expect_equal(ref_row$level, "Female")

  non_ref_row <- res$results_table[!res$results_table$is_reference, ]
  expect_equal(non_ref_row$level, "Male")

  # HR should be (approximately) the reciprocal of the default-reference fit
  res_default <- cox_forest(
    data = lung2, covariates = "sex", labels = "Sex",
    time = "time_os", event = "event_os", model_type = "univariable"
  )
  hr_default <- res_default$results_table$hr[!res_default$results_table$is_reference]
  hr_flipped <- res$results_table$hr[!res$results_table$is_reference]
  expect_equal(hr_flipped, 1 / hr_default, tolerance = 1e-8)
})

test_that("an invalid `ref_levels` entry errors with an informative message", {
  expect_error(
    cox_forest(
      data = lung2, covariates = "sex", labels = "Sex",
      time = "time_os", event = "event_os",
      ref_levels = list(sex = "Not a level")
    ),
    "not found among observed levels"
  )
})

test_that("`exclude_levels` drops the specified level from the univariable model", {
  res <- cox_forest(
    data = lung2, covariates = "ph.ecog", labels = "ECOG PS",
    time = "time_os", event = "event_os", model_type = "univariable",
    exclude_levels = list(ph.ecog = "3")
  )
  expect_false("3" %in% res$results_table$level)
})

# ==============================================================================
# Multiple covariates in one call (univariable = one model per covariate)
# ==============================================================================

test_that("multiple covariates each get their own univariable model", {
  res <- cox_forest(
    data = lung2, covariates = c("sex", "age", "ph.ecog"),
    labels = c("Sex", "Age", "ECOG PS"),
    time = "time_os", event = "event_os",
    continuous_vars = "age", model_type = "univariable"
  )

  expect_setequal(unique(res$results_table$covariate), c("sex", "age", "ph.ecog"))
  # age (continuous) contributes exactly 1 row
  expect_equal(sum(res$results_table$covariate == "age"), 1)
  # sex and ph.ecog (categorical) contribute one row per observed level
  expect_equal(sum(res$results_table$covariate == "sex"), nlevels(lung2$sex))
  expect_equal(sum(res$results_table$covariate == "ph.ecog"), nlevels(lung2$ph.ecog))
})

# ==============================================================================
# Multivariable models
# ==============================================================================

test_that("multivariable model fits a single model containing all covariates", {
  res <- cox_forest(
    data = lung2, covariates = c("sex", "age", "ph.ecog"),
    labels = c("Sex", "Age", "ECOG PS"),
    time = "time_os", event = "event_os",
    continuous_vars = "age", model_type = "multivariable"
  )

  expect_true(all(res$results_table$model_type == "multivariable"))

  fit <- survival::coxph(
    survival::Surv(time_os, event_os) ~ sex + age + ph.ecog, data = lung2
  )
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  age_row <- res$results_table[res$results_table$covariate == "age", ]
  age_hr_manual <- tt$estimate[tt$term == "age"]
  expect_equal(age_row$hr, age_hr_manual, tolerance = 1e-8)
})

test_that("`model_type = \"both\"` returns both model types and both plots", {
  res <- cox_forest(
    data = lung2, covariates = c("sex", "age"), labels = c("Sex", "Age"),
    time = "time_os", event = "event_os",
    continuous_vars = "age", model_type = "both"
  )

  expect_setequal(unique(res$results_table$model_type), c("univariable", "multivariable"))
  expect_named(res$plots, c("univariable", "multivariable"))
  expect_true(all(vapply(res$plots, function(p) "gtable" %in% class(p), logical(1))))
})

test_that("`exclude_levels` drops affected patients from the whole multivariable model", {
  res_full <- cox_forest(
    data = lung2, covariates = c("sex", "ph.ecog"), labels = c("Sex", "ECOG PS"),
    time = "time_os", event = "event_os", model_type = "multivariable"
  )
  res_excl <- cox_forest(
    data = lung2, covariates = c("sex", "ph.ecog"), labels = c("Sex", "ECOG PS"),
    time = "time_os", event = "event_os", model_type = "multivariable",
    exclude_levels = list(ph.ecog = "3")
  )

  n_excl_ecog <- res_excl$results_table$n[res_excl$results_table$covariate == "ph.ecog" &
                                             res_excl$results_table$is_reference]
  n_full_ecog <- res_full$results_table$n[res_full$results_table$covariate == "ph.ecog" &
                                             res_full$results_table$is_reference]
  # Excluding an ECOG level should never increase the analysable N.
  expect_true(n_excl_ecog <= n_full_ecog)
  expect_false("3" %in% res_excl$results_table$level)
})

# ==============================================================================
# Multiple endpoints
# ==============================================================================

test_that("multiple endpoints produce results for each endpoint separately", {
  res <- cox_forest(
    data = lung2, covariates = c("sex", "age"), labels = c("Sex", "Age"),
    time = c("time_os", "time_pfs"), event = c("event_os", "event_pfs"),
    endpoint_labels = c("OS", "PFS"),
    continuous_vars = "age", model_type = "univariable"
  )

  expect_setequal(unique(res$results_table$endpoint_label), c("OS", "PFS"))
  n_single_endpoint <- cox_forest(
    data = lung2, covariates = c("sex", "age"), labels = c("Sex", "Age"),
    time = "time_os", event = "event_os",
    continuous_vars = "age", model_type = "univariable"
  )$results_table |> nrow()
  expect_equal(nrow(res$results_table), 2 * n_single_endpoint)
})

test_that("multi-endpoint plots render without error and use a custom palette", {
  res <- cox_forest(
    data = lung2, covariates = c("sex", "age"), labels = c("Sex", "Age"),
    time = c("time_os", "time_pfs"), event = c("event_os", "event_pfs"),
    endpoint_labels = c("OS", "PFS"), continuous_vars = "age",
    model_type = "univariable", palette = c("#1B9E77", "#D95F02")
  )
  expect_true("gtable" %in% class(res$plots$univariable))
})

# ==============================================================================
# Plot styling arguments
# ==============================================================================

test_that("`row_height` scales the body row heights", {
  res_default <- cox_forest(
    data = lung2, covariates = "sex", labels = "Sex",
    time = "time_os", event = "event_os", model_type = "univariable"
  )
  res_tall <- cox_forest(
    data = lung2, covariates = "sex", labels = "Sex",
    time = "time_os", event = "event_os", model_type = "univariable",
    row_height = 1.5
  )

  p1 <- res_default$plots$univariable
  p2 <- res_tall$plots$univariable
  core_fg1 <- which(p1$layout$name == "core-fg")
  core_fg2 <- which(p2$layout$name == "core-fg")
  body_rows1 <- unique(p1$layout$t[core_fg1])
  body_rows2 <- unique(p2$layout$t[core_fg2])

  h1 <- sum(grid::convertHeight(p1$heights[body_rows1], "in", valueOnly = TRUE))
  h2 <- sum(grid::convertHeight(p2$heights[body_rows2], "in", valueOnly = TRUE))
  expect_gt(h2, h1)
  expect_equal(h2 / h1, 1.5, tolerance = 0.05)
})

test_that("`covariate_label_row = TRUE` adds one extra header row per covariate", {
  covariates <- c("sex", "ph.ecog")
  res <- cox_forest(
    data = lung2, covariates = covariates, labels = c("Sex", "ECOG PS"),
    time = "time_os", event = "event_os", model_type = "univariable",
    covariate_label_row = TRUE
  )
  res_no_header <- cox_forest(
    data = lung2, covariates = covariates, labels = c("Sex", "ECOG PS"),
    time = "time_os", event = "event_os", model_type = "univariable",
    covariate_label_row = FALSE
  )

  p_header <- res$plots$univariable
  p_plain  <- res_no_header$plots$univariable

  n_body_header <- length(unique(p_header$layout$t[p_header$layout$name == "core-fg"]))
  n_body_plain  <- length(unique(p_plain$layout$t[p_plain$layout$name == "core-fg"]))

  expect_equal(n_body_header, n_body_plain + length(covariates))
})

test_that("a custom title is honoured on the plot object", {
  res <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os", continuous_vars = "age",
    model_type = "univariable", title = "My custom title"
  )
  expect_true("gtable" %in% class(res$plots$univariable))
})

test_that("custom `xlim`/`ticks_at` are accepted and produce a valid plot", {
  res <- cox_forest(
    data = lung2, covariates = "age", labels = "Age",
    time = "time_os", event = "event_os", continuous_vars = "age",
    model_type = "univariable", xlim = c(0.5, 2), ticks_at = c(0.5, 1, 1.5, 2)
  )
  expect_true("gtable" %in% class(res$plots$univariable))
})

# ==============================================================================
# Internal helper functions
# ==============================================================================

test_that(".format_hr_ci formats numeric estimates and flags NA as Reference", {
  expect_equal(.format_hr_ci(1.234, 1.01, 1.50, digits = 2), "1.23 (1.01-1.50)")
  expect_equal(.format_hr_ci(NA, NA, NA), "Reference")
})

test_that(".format_pvalue applies the < threshold and fixed-digit formatting", {
  expect_equal(.format_pvalue(0.0001, digits = 3), "<0.001")
  expect_equal(.format_pvalue(0.0456, digits = 3), "0.046")
  expect_equal(.format_pvalue(NA), "")
})

test_that(".detect_covariate_type distinguishes categorical vs continuous", {
  expect_equal(.detect_covariate_type(lung2, "sex"), "categorical")
  expect_equal(.detect_covariate_type(lung2, "age"), "continuous")
  # explicit override forces categorical treatment even for numeric data
  expect_equal(.detect_covariate_type(lung2, "age", continuous_vars = character(0)), "continuous")
})

test_that(".factor_with_ref relevels correctly and errors on an unknown level", {
  f <- .factor_with_ref(lung2$sex, ref_level = "Female")
  expect_equal(levels(f)[1], "Female")
  expect_error(.factor_with_ref(lung2$sex, ref_level = "Unknown"), "not found among observed levels")
})

test_that(".apply_exclude_levels sets matching values to NA", {
  d <- .apply_exclude_levels(lung2, list(ph.ecog = "3"))
  expect_true(all(is.na(d$ph.ecog[as.character(lung2$ph.ecog) == "3"])))
  expect_equal(sum(is.na(d$ph.ecog)), sum(as.character(lung2$ph.ecog) == "3"))
})

test_that(".apply_exclude_levels is a no-op when `exclude_levels` is NULL", {
  d <- .apply_exclude_levels(lung2, NULL)
  expect_identical(d, lung2)
})
