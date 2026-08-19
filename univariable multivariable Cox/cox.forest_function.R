# ==============================================================================
# cox_forest.R
#
# Publication-quality forest plots for standard univariable and/or
# multivariable Cox proportional hazards models, across one or more
# endpoints.
#
# For every covariate (which
# may now include a treatment/arm variable, if included in `covariates`), the
# HR for each
# non-reference level is estimated either
#   - univariable: one Cox model per covariate, i.e. Surv(time, event) ~ covariate
#   - multivariable: one Cox model per endpoint containing ALL covariates
#     together, i.e. Surv(time, event) ~ cov1 + cov2 + ... + covk
# and the results are rendered as a forest plot with one row per covariate
# level (a "Reference" row with HR = 1 is shown for the reference level of
# every categorical covariate).
#
# Author: Julian Silan
# Requires: survival, broom, dplyr, tidyr, purrr, forestploter, tibble,
#           stringr, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(broom)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(forestploter)
  library(tibble)
  library(stringr)
  library(grid)
})


# ==============================================================================
# Helper functions
# ==============================================================================

#' @keywords internal
.validate_cox_inputs <- function(data, covariates, labels, time, event) {
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame (or tibble).", call. = FALSE)
  }
  if (length(covariates) == 0) {
    stop("`covariates` must contain at least one variable name.", call. = FALSE)
  }
  if (length(labels) != length(covariates)) {
    stop(sprintf(
      "`labels` (length %d) and `covariates` (length %d) must have the same length.",
      length(labels), length(covariates)
    ), call. = FALSE)
  }
  if (length(time) != length(event)) {
    stop(sprintf(
      "`time` (length %d) and `event` (length %d) must have the same length.",
      length(time), length(event)
    ), call. = FALSE)
  }
  
  all_vars <- unique(c(covariates, time, event))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variable(s) were not found in `data`: %s",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Decide whether a covariate should be treated as categorical or continuous
#' @keywords internal
.detect_covariate_type <- function(data, covariate, continuous_vars = NULL) {
  if (!is.null(continuous_vars) && covariate %in% continuous_vars) return("continuous")
  x <- data[[covariate]]
  if (is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) > 5) return("continuous")
  "categorical"
}

#' Coerce a variable to a factor, moving `ref_level` to the front if supplied
#' @keywords internal
.factor_with_ref <- function(x, ref_level = NULL) {
  f <- if (is.factor(x)) droplevels(x) else factor(as.character(x))
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(f)) {
      stop(sprintf("Reference level '%s' not found among observed levels (%s).",
                   ref_level, paste(levels(f), collapse = ", ")), call. = FALSE)
    }
    f <- stats::relevel(f, ref = ref_level)
  }
  f
}

#' Format an HR (95% CI) text string, e.g. "1.34 (1.05-1.72)"
#' @keywords internal
.format_hr_ci <- function(est, lower, upper, digits = 2) {
  # ASCII hyphen (not an en dash) is used deliberately: the base R "pdf"
  # device's default fonts use WinAnsi/Latin-1 encoding and silently corrupt
  # true en dashes. Use a cairo-based device if a typographic dash is wanted.
  fmt <- paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)")
  ifelse(is.na(est) | is.na(lower) | is.na(upper), "Reference",
         sprintf(fmt, est, lower, upper))
}

#' Format a p-value to `digits` decimals, showing "<0.001" style thresholds
#' @keywords internal
.format_pvalue <- function(p, digits = 3) {
  threshold <- 10^(-digits)
  dplyr::case_when(
    is.na(p)      ~ "",
    p < threshold ~ paste0("<", format(threshold, scientific = FALSE)),
    TRUE          ~ formatC(p, digits = digits, format = "f")
  )
}

#' Exclude specific levels of specific covariates from analysis, by setting
#' matching values to NA. This lets each covariate's excluded level(s) drop
#' out via the existing complete-case filtering: for univariable models only
#' that covariate's own model is affected, while for the multivariable model
#' the affected patients are dropped from the whole model (since it requires
#' complete cases across all covariates simultaneously).
#' @keywords internal
.apply_exclude_levels <- function(data, exclude_levels) {
  if (is.null(exclude_levels) || length(exclude_levels) == 0) return(data)
  for (cov in names(exclude_levels)) {
    if (!cov %in% names(data)) next
    lvls_to_drop <- exclude_levels[[cov]]
    data[[cov]][as.character(data[[cov]]) %in% lvls_to_drop] <- NA
  }
  data
}

#' N / events for one level of one covariate within a given (complete-case)
#' dataset
#' @keywords internal
.level_counts <- function(data, covariate, event_var, level = NULL) {
  d <- if (is.null(level)) data else data[as.character(data[[covariate]]) == level, , drop = FALSE]
  list(n = nrow(d), n_event = suppressWarnings(sum(d[[event_var]], na.rm = TRUE)))
}

#' Build the rows (incl. a "Reference" row) for one categorical covariate
#' from a fitted coxph model, given the term-name prefix used by R for that
#' covariate's factor levels (e.g. "sex" -> "sexMale").
#' @keywords internal
.rows_from_categorical_term <- function(tidy_fit, model_data, covariate, covariate_label,
                                        event_var, fac_levels, ref_level) {
  
  non_ref_levels <- setdiff(fac_levels, ref_level)
  term_names <- paste0(covariate, non_ref_levels)
  
  level_rows <- purrr::map_dfr(seq_along(non_ref_levels), function(i) {
    lev <- non_ref_levels[i]
    trm <- term_names[i]
    hit <- tidy_fit[tidy_fit$term == trm, , drop = FALSE]
    cnt <- .level_counts(model_data, covariate, event_var, lev)
    tibble::tibble(
      covariate = covariate, covariate_label = covariate_label,
      level = lev, is_reference = FALSE,
      n = cnt$n, n_event = cnt$n_event,
      hr = if (nrow(hit) == 1) hit$estimate[1] else NA_real_,
      conf_low = if (nrow(hit) == 1) hit$conf.low[1] else NA_real_,
      conf_high = if (nrow(hit) == 1) hit$conf.high[1] else NA_real_,
      p_value = if (nrow(hit) == 1) hit$p.value[1] else NA_real_
    )
  })
  
  ref_cnt <- .level_counts(model_data, covariate, event_var, ref_level)
  ref_row <- tibble::tibble(
    covariate = covariate, covariate_label = covariate_label,
    level = ref_level, is_reference = TRUE,
    n = ref_cnt$n, n_event = ref_cnt$n_event,
    hr = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_
  )
  
  dplyr::bind_rows(ref_row, level_rows)
}

#' Build the (single) row for a continuous covariate from a fitted coxph
#' model
#' @keywords internal
.row_from_continuous_term <- function(tidy_fit, model_data, covariate, covariate_label, event_var,
                                      increase = 1) {
  hit <- tidy_fit[tidy_fit$term == covariate, , drop = FALSE]
  cnt <- .level_counts(model_data, covariate, event_var, level = NULL)
  
  tibble::tibble(
    covariate = covariate, covariate_label = covariate_label,
    level = paste0("per ", increase, "-unit increase"), 
    is_reference = FALSE,
    n = cnt$n, n_event = cnt$n_event,
    hr = if (nrow(hit) == 1) hit$estimate[1]^increase else NA_real_,
    conf_low = if (nrow(hit) == 1) hit$conf.low[1]^increase else NA_real_,
    conf_high = if (nrow(hit) == 1) hit$conf.high[1]^increase else NA_real_,
    p_value = if (nrow(hit) == 1) hit$p.value[1] else NA_real_
  )
}

#' Fit ONE univariable Cox model (Surv(time, event) ~ covariate) and return
#' its result rows
#' @keywords internal
.fit_univariable <- function(data, time_var, event_var, covariate, cov_type,
                             ref_level, weights, continuous_increase = 1) {
  
  needed <- c(time_var, event_var, covariate)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  
  if (cov_type == "categorical") {
    d[[covariate]] <- .factor_with_ref(d[[covariate]], ref_level)
    fac_levels <- levels(d[[covariate]])
    ref_level  <- fac_levels[1]
  }
  
  if (nrow(d) == 0) return(tibble::tibble())
  
  fmla <- stats::as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", covariate))
  fit <- tryCatch(
    if (is.null(weights)) survival::coxph(fmla, data = d)
    else survival::coxph(fmla, data = d, weights = get(weights)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble::tibble())
  
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  if (cov_type == "categorical") {
    .rows_from_categorical_term(tt, d, covariate, covariate, event_var, fac_levels, ref_level)
  } else {
    .row_from_continuous_term(tt, d, covariate, covariate, event_var, increase = continuous_increase)
  }
}

#' Fit ONE multivariable Cox model (Surv(time, event) ~ cov1 + ... + covk)
#' and return the result rows for ALL covariates at once
#' @keywords internal
.fit_multivariable <- function(data, time_var, event_var, covariates, cov_types,
                               ref_levels, weights, continuous_increases = NULL) {
  
  needed <- c(time_var, event_var, covariates)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  
  fac_level_map <- list()
  for (cov in covariates) {
    if (cov_types[[cov]] == "categorical") {
      d[[cov]] <- .factor_with_ref(d[[cov]], ref_levels[[cov]])
      fac_level_map[[cov]] <- levels(d[[cov]])
    }
  }
  
  if (nrow(d) == 0) return(tibble::tibble())
  
  fmla <- stats::as.formula(paste0(
    "survival::Surv(", time_var, ", ", event_var, ") ~ ", paste(covariates, collapse = " + ")
  ))
  fit <- tryCatch(
    if (is.null(weights)) survival::coxph(fmla, data = d)
    else survival::coxph(fmla, data = d, weights = get(weights)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble::tibble())
  
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  purrr::map_dfr(covariates, function(cov) {
    if (cov_types[[cov]] == "categorical") {
      .rows_from_categorical_term(tt, d, cov, cov, event_var,
                                  fac_level_map[[cov]], fac_level_map[[cov]][1])
    } else {
      increase <- if (!is.null(continuous_increases[[cov]])) {
        continuous_increases[[cov]]
      } else {
        1
      }
      
      .row_from_continuous_term(
        tt, d, cov, cov, event_var,
        increase = increase
      )
    }
  })
}

#' Apply a row-height multiplier to the body rows of a forestploter object
#'
#' @param p A forestploter/gtable object.
#' @param row_height Numeric multiplier. 1 = default height, 1.3 = 30% taller.
#' @keywords internal
.apply_row_height <- function(p, row_height = 1) {
  
  if (length(row_height) != 1 || !is.numeric(row_height) ||
      is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }
  
  # Locate the body rows. forestploter stores the table as a gtable.
  core_fg <- which(p$layout$name == "core-fg")
  
  if (length(core_fg) == 0) {
    return(p)
  }
  
  body_rows <- unique(p$layout$t[core_fg])
  
  # Only enlarge body rows, not the header / title / legend.
  p$heights[body_rows] <- p$heights[body_rows] * row_height
  
  p
}


#' Build the full long-format results table for one model type
#' ("univariable" or "multivariable"), across all endpoints
#' @keywords internal
.build_cox_results <- function(data, covariates, labels, time, event, endpoint_labels,
                               cov_types, ref_levels, weights, model_type, 
                               continuous_increases = NULL) {
  
  label_lookup <- stats::setNames(labels, covariates)
  
  purrr::map_dfr(seq_along(time), function(ep) {
    time_var  <- time[ep]
    event_var <- event[ep]
    ep_label  <- endpoint_labels[ep]
    
    if (model_type == "univariable") {
      rows <- purrr::map_dfr(covariates, function(cov) {
        increase <- if (!is.null(continuous_increases[[cov]])) {continuous_increases[[cov]]} else {1}

        .fit_univariable(data, time_var, event_var, cov, cov_types[[cov]],
                         ref_levels[[cov]], weights, continuous_increase = increase)
      })
    } else {
      rows <- .fit_multivariable(data, time_var, event_var, covariates, cov_types,
                                 ref_levels, weights, continuous_increases = continuous_increases)
    }
    
    if (nrow(rows) == 0) return(rows)
    
    rows |>
      dplyr::mutate(
        covariate_label = unname(label_lookup[covariate]),
        endpoint = time_var, endpoint_label = ep_label,
        model_type = model_type
      )
  })
}

#' Build a forestploter forest plot object for one model's results table
#' (rows grouped covariate-major; a bold "Reference" row per categorical
#' covariate; colour-coded by endpoint when there is more than one)
#' @keywords internal
.build_cox_forest_plot <- function(results, digits, endpoint_palette, xlim, ticks_at,
                                   length_plot, ci_pch, ci_lty, ci_lwd, ci_Theight, sizes,
                                   risk_labels, title, row_height, covariate_label_row = FALSE) {
  
  if (length(row_height) != 1 || !is.numeric(row_height) || is.na(row_height) || row_height <= 0) {
    stop("`row_height` must be a single positive numeric value.", call. = FALSE)
  }
  
  if (is.null(endpoint_palette)) {
    endpoint_palette <- c("#E67E22", "#16A085", "#8E44AD", "#C0392B", "#2980B9")
  }
  
  # Preserve covariate order as first encountered, then order rows
  # covariate-major so each covariate's levels stay adjacent across endpoints.
  cov_order <- unique(results$covariate_label)
  df <- results |>
    dplyr::mutate(covariate_label = factor(covariate_label, levels = cov_order)) |>
    dplyr::arrange(covariate_label, desc(is_reference), level, endpoint_label) |>
    # dplyr::arrange(covariate_label, level, endpoint_label) |>
    dplyr::group_by(covariate_label, level, is_reference) |>
    dplyr::mutate(endpoint_order = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # is_first_of_block = covariate_label != dplyr::lag(covariate_label, default = factor(NA, levels = cov_order)),
      # covariate_display  = ifelse(is_first_of_block, as.character(covariate_label), ""),
      is_first_of_block = dplyr::row_number() == 1 | covariate_label != dplyr::lag(covariate_label),
      covariate_display = ifelse(is_first_of_block, as.character(covariate_label), ""),
      level_display = ifelse(endpoint_order == 1, level, ""),
      n_fmt = ifelse(endpoint_order == 1 & !is.na(n), as.character(n), ""),
     # level_display      = level, #ifelse(is_reference, paste0(level, " (ref.)"), level),
      block_id           = cumsum(is_first_of_block),
      hr_ci        = .format_hr_ci(hr, conf_low, conf_high, digits),
      p_fmt        = .format_pvalue(p_value),
      # n_fmt        = ifelse(is.na(n), "", as.character(n)),
      event_fmt    = ifelse(is.na(n_event), "", as.character(n_event)),
      blank_col    = strrep(" ", length_plot)
    )

  
  # Add optional covariate label rows
  if (covariate_label_row) {
    df$is_covariate_header <- FALSE
    # Split by covariate and insert one header row immediately before
    # the corresponding covariate block.
    df <- purrr::map_dfr(
      split(df, df$covariate_label),
      function(block) {
        
        header <- block[1, , drop = FALSE]
        
        header$covariate_display <- as.character(block$covariate_label[1])
        header$level_display <- ""
        header$n_fmt <- ""
        header$event_fmt <- ""
        header$blank_col <- strrep(" ", length_plot)
        header$hr_ci <- ""
        header$p_fmt <- ""
        header$hr <- NA_real_
        header$conf_low <- NA_real_
        header$conf_high <- NA_real_
        header$is_covariate_header <- TRUE
        
        block$covariate_display <- ""
        block$is_covariate_header <- FALSE
        
        dplyr::bind_rows(header, block)
      }
    )} else {
      df$is_covariate_header <- FALSE
    }
  
  endpoints      <- unique(df$endpoint_label)
  multi_endpoint <- length(endpoints) > 1
  endpoint_col   <- stats::setNames(rep(endpoint_palette, length.out = length(endpoints)), endpoints)
  
  plot_cols <- c("covariate_display", "level_display", "n_fmt", "event_fmt",
                 "blank_col", "hr_ci", "p_fmt")
  tbl <- as.data.frame(df[, plot_cols])
  names(tbl) <- c("Variable", "Level", "N", "Events", " ", "HR (95% CI)", "p-value")
  ci_col <- which(names(tbl) == " ")
  
  # Reference rows carry no estimate -> forestploter draws no point/CI for
  # them automatically, which is exactly what we want.
  if (is.null(xlim)) {
    finite_bounds <- c(df$conf_low, df$conf_high, 1)
    finite_bounds <- finite_bounds[is.finite(finite_bounds) & finite_bounds > 0]
    x_range <- range(finite_bounds, na.rm = TRUE)
    xlim <- c(max(0.05, x_range[1] * 0.8), x_range[2] * 1.25)
  }
  if (is.null(ticks_at)) {
    ticks_at <- pretty(xlim, n = 5)
    ticks_at <- ticks_at[ticks_at > 0 & ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    if (length(ticks_at) < 2) ticks_at <- signif(exp(seq(log(xlim[1]), log(xlim[2]), length.out = 4)), 2)
  }
  
  
  # Alternating row shading
  # IMPORTANT:
  # With covariate_label_row = TRUE, the header row belongs to the
  # corresponding covariate block. Therefore the complete block
  # (header + all levels) gets the same background.

  if (covariate_label_row) {
    covariate_block_id <- cumsum(df$is_covariate_header)
    shaded_rows <- which(covariate_block_id %% 2 == 1)
  } else {
    shaded_rows <- which(df$block_id %% 2 == 1)
  }
  
  ft_theme <- forestploter::forest_theme(
    base_size  = 11,
    ci_pch     = ci_pch,
    ci_col     = if (multi_endpoint) endpoint_col[[1]] else endpoint_palette[[1]],
    ci_lty     = ci_lty,
    ci_lwd     = ci_lwd,
    ci_Theight = ci_Theight,
    axis_extend = c(0.2, 0.05),
    vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20",
    refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "grey30"),
    footnote_gp = grid::gpar(col = "grey20", fontface = "italic", cex = 0.8),
    core = list(bg_params = list(fill = "white"))
  )
  
  arrow_lab <- if (is.null(risk_labels)) NULL else risk_labels
  
  p <- forestploter::forest(
    tbl,
    est       = df$hr,
    lower     = df$conf_low,
    upper     = df$conf_high,
    sizes     = sizes,
    ci_column = ci_col,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = xlim,
    ticks_at  = ticks_at,
    arrow_lab = arrow_lab,
    title     = title,
    theme     = ft_theme
  )

  p <- .apply_row_height(p, row_height = row_height)
  
  # ==========================================================================
  # Bold covariate header rows
  # ==========================================================================
  
  if (covariate_label_row) {
    covariate_header_rows <- which(df$is_covariate_header)
    if (length(covariate_header_rows) > 0) {
      p <- forestploter::edit_plot(p, row = covariate_header_rows, col = 1, 
                                   which = "text", gp = grid::gpar(fontface = "bold"))
    }
  } else {
    # Original behaviour:
    # bold the first row of each covariate block.
    covlabel_rows <- which(df$is_first_of_block)
    if (length(covlabel_rows) > 0) {
      p <- forestploter::edit_plot(p, row = covlabel_rows, col = 1,
                                   gp = grid::gpar(fontface = "bold"))
    }
  }
  
  
  
  if (length(shaded_rows) > 0) {
    p <- forestploter::edit_plot(p, row = shaded_rows, which = "background",
                                 gp = grid::gpar(fill = "grey95", col = "grey95"))
  }
  
  # ref_rows <- which(df$is_reference)
  # if (length(ref_rows) > 0) {
  #   p <- forestploter::edit_plot(p, row = ref_rows, gp = grid::gpar(fontface = "italic"))
  # }

  if (multi_endpoint) {
    for (ep_label in endpoints) {
      rows_i <- which(df$endpoint_label == ep_label & !df$is_reference)
      col_i  <- endpoint_col[[ep_label]]
      if (length(rows_i) > 0) {
        p <- forestploter::edit_plot(p, row = rows_i, col = ci_col, which = "ci",
                                     gp = grid::gpar(col = col_i, fill = col_i))
      }
    }
  }
  
  # Vertical separator lines between table columns
  for (col_i in c(0, ci_col-1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p,
      row = 1:nrow(tbl),
      col = col_i,
      part = "body",
      where = "right",
      gp = grid::gpar(
        col = "grey30",
        lwd = 1
      )
    )
  }
  
  # Vertical separator lines in header
  for (col_i in c(0, ci_col-1, ci_col, length(names(tbl)))) {
    p <- forestploter::add_border(
      p,
      row = 1,
      col = col_i,
      part = "header",
      where = "right",
      gp = grid::gpar(
        col = "grey30",
        lwd = 1
      )
    )
  }
  
  # Horizontal lines 
  for (row_i in c(0, nrow(tbl))) {
    p <- forestploter::add_border(
      p,
      row = row_i,
      col = 1:ncol(tbl),
      part = "body",
      where = "bottom",
      gp = grid::gpar(
        col = "grey30",
        lwd = 1
      )
    )
  }
  
  for (row_i in c(0)) {
    p <- forestploter::add_border(
      p,
      row = row_i,
      col = 1:ncol(tbl),
      part = "header",
      where = "bottom",
      gp = grid::gpar(
        col = "grey30",
        lwd = 1
      )
    )
  }
  
  # Manual legend for multiple endpoints
  if (multi_endpoint) {

    legend_cols <- unname(endpoint_col[endpoints])
    n_legend <- length(endpoints)

    # Breite pro Legendeneintrag
    entry_width <- 1.4

    # Gesamtbreite der Legende
    total_width <- n_legend * entry_width

    # Erzeuge Punkt + Label für jeden Endpoint
    legend_grobs <- list()

    for (i in seq_along(endpoints)) {

      # x-Position des Eintrags
      x_pos <- (i - 0.5) / n_legend

      # Punkt
      legend_grobs[[length(legend_grobs) + 1]] <- grid::pointsGrob(
        x = grid::unit(x_pos - 0.035, "npc"),
        y = grid::unit(0.5, "npc"),
        pch = ci_pch,
        size = grid::unit(2.2, "mm"),
        gp = grid::gpar(
          col = legend_cols[i],
          fill = legend_cols[i]
        )
      )

      # Label
      legend_grobs[[length(legend_grobs) + 1]] <- grid::textGrob(
        label = endpoints[i],
        x = grid::unit(x_pos - 0.005, "npc"),
        y = grid::unit(0.5, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(
          fontsize = 9,
          col = "grey20"
        )
      )
    }

    legend_grob <- grid::gTree(
      children = do.call(grid::gList, legend_grobs)
    )

    # Zusätzliche Zeile unterhalb des Forest Plots
    p <- gtable::gtable_add_rows(
      p,
      heights = grid::unit(0.3, "in"),
      pos = nrow(p)
    )

    # Legende über die gesamte Breite
    p <- gtable::gtable_add_grob(
      p,
      legend_grob,
      t = nrow(p),
      l = 1,
      r = ncol(p),
      name = "manual_endpoint_legend"
    )
  }
  
  p
}


# ==============================================================================
# Main function
# ==============================================================================

#' Publication-quality forest plots for univariable and/or multivariable Cox
#' models
#'
#' For every covariate in \code{covariates} (optionally including
#' a treatment/arm variable, if included), estimates the hazard ratio for each
#' non-reference level, for every requested endpoint, and renders the result
#' as a forest plot in the style of a standard baseline/prognostic-factor
#' Cox regression table.
#'
#' Two model types are supported (see \code{model_type}):
#' \describe{
#'   \item{univariable}{One Cox model per covariate:
#'     \code{Surv(time, event) ~ covariate}.}
#'   \item{multivariable}{One Cox model per endpoint containing all
#'     covariates simultaneously:
#'     \code{Surv(time, event) ~ cov1 + cov2 + ... + covk}.}
#' }
#'
#' @param data A data.frame or tibble.
#' @param covariates Character vector of covariate names in \code{data}.
#' @param labels Character vector of display labels, same length/order as
#'   \code{covariates}.
#' @param time Character vector of one or more survival time variable names.
#' @param event Character vector of event indicator variable names, paired
#'   by position with \code{time}.
#' @param endpoint_labels Optional display labels for each (time, event)
#'   pair. Defaults to \code{time}.
#' @param continuous_vars Optional character vector identifying which
#'   \code{covariates} should be treated as continuous. If \code{NULL},
#'   numeric variables with more than 5 distinct values are auto-detected as
#'   continuous. To include a treatment/arm variable, just add it to
#'   \code{covariates} (and \code{labels}) like any other covariate.
#' @param ref_levels Optional named list, e.g. \code{list(sex = "Female")},
#'   giving the reference level for one or more categorical covariates.
#'   Defaults to each variable's first factor level (or first level in
#'   sorted order).
#' @param exclude_levels Optional named list, e.g.
#'   \code{list(ecog = "Unknown", smoking = c("Unknown", "Not reported"))},
#'   giving one or more levels to exclude per covariate. Matching values are
#'   treated as missing: for that covariate's univariable model only those
#'   patients are dropped, while for the multivariable model they are
#'   dropped from the whole model (since it requires complete cases across
#'   all covariates simultaneously).
#' @param model_type One of \code{"univariable"}, \code{"multivariable"}, or
#'   \code{"both"}. Default \code{"univariable"}.
#' @param digits Number of decimal places for HR/CI display. Default 2.
#' @param palette Colours used per endpoint (only relevant with >1
#'   endpoint). Must be the same length as \code{endpoint_labels}.
#' @param risk_labels Optional length-2 character vector for the axis arrow
#'   annotations, e.g. \code{c("Lower risk", "Higher risk")}. \code{NULL}
#'   (default) omits the arrows, since with multiple, unrelated covariates
#'   there is no single "favours A / favours B" framing.
#' @param row_height Numeric multiplier controlling the height of body rows.
#'   `1` uses the default forestploter row height; values above 1 increase
#'   the row height proportionally, e.g. `1.3` makes rows approximately 30%
#'   taller.
#' @param covariate_label_row Logical. If `TRUE`, a separate bold row
#'   containing the covariate label is inserted before each covariate block.
#'   Default is `FALSE`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{results_table}{Long-format tibble, one row per covariate level x
#'     endpoint x model type, with HR, 95% CI, p-value, N and events.}
#'   \item{plots}{Named list with up to two forestploter plot objects:
#'     \code{$univariable} and/or \code{$multivariable}, depending on
#'     \code{model_type}.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cox_forest(
#'   data = my_trial_data,
#'   covariates = c("arm", "sex", "age", "ecog"),
#'   labels = c("Treatment", "Sex", "Age", "ECOG PS"),
#'   time = c("os_time", "pfs_time"),
#'   event = c("os_event", "pfs_event"),
#'   endpoint_labels = c("OS", "PFS"),
#'   continuous_vars = "age",
#'   exclude_levels = list(ecog = "Unknown"),
#'   model_type = "both",
#'   palette = c("#E67E22", "#16A085")
#' )
#' result$plots$univariable
#' result$plots$multivariable
#' }
#' @export
cox_forest <- function(data, covariates, labels, time, event, endpoint_labels = NULL,
                       continuous_vars = NULL, continuous_increase = NULL,
                       ref_levels = NULL, exclude_levels = NULL,
                       model_type = c("univariable", "multivariable", "both"),
                       digits = 2, palette = NULL, xlim = NULL, ticks_at = NULL,
                       length_plot = 40, weights = NULL, ci_pch = 15, ci_lty = 1,
                       ci_lwd = 3, ci_Theight = 0.4, sizes = 0.5, risk_labels = NULL,
                       title = NULL, row_height = 1, covariate_label_row = FALSE) {
  
  model_type <- match.arg(model_type, several.ok = TRUE)
  if ("both" %in% model_type) model_type <- c("univariable", "multivariable")
  model_type <- unique(model_type)
  
  .validate_cox_inputs(data, covariates, labels, time, event)
  
  if (is.null(endpoint_labels)) endpoint_labels <- time
  if (length(endpoint_labels) != length(time)) {
    stop("`endpoint_labels` must be the same length as `time`.", call. = FALSE)
  }
  if (!is.null(palette) && length(endpoint_labels) != length(palette)) {
    stop("`palette` must be the same length as `endpoint_labels`.", call. = FALSE)
  }
  
  data <- .apply_exclude_levels(data, exclude_levels)
  
  cov_types <- stats::setNames(
    lapply(covariates, function(cov) .detect_covariate_type(data, cov, continuous_vars)),
    covariates
  )
  if (is.null(ref_levels)) ref_levels <- list()
  
  results_table <- purrr::map_dfr(model_type, function(mt) {
    .build_cox_results(data, covariates, labels, time, event, endpoint_labels,
                       cov_types, ref_levels, weights, mt, 
                       continuous_increases = continuous_increase)
  })
  
  plots <- list()
  for (mt in model_type) {
    res_mt <- dplyr::filter(results_table, model_type == mt)
    if (nrow(res_mt) == 0) next
    plot_title <- if (is.null(title)) {
      if (mt == "univariable") "Univariable Cox regression" else "Multivariable Cox regression"
    } else title
    plots[[mt]] <- .build_cox_forest_plot(
      results = res_mt, digits = digits, endpoint_palette = palette, xlim = xlim,
      ticks_at = ticks_at, length_plot = length_plot, ci_pch = ci_pch, ci_lty = ci_lty,
      ci_lwd = ci_lwd, ci_Theight = ci_Theight, sizes = sizes, risk_labels = risk_labels,
      title = plot_title, row_height = row_height, covariate_label_row = covariate_label_row
    )
  }
  
  list(
    results_table = results_table,
    plots         = plots
  )
}