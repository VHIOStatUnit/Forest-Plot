rm(list = ls(all = TRUE))

library(tidyverse)

source("~/Desktop/VHIO/R code/Forestplot/forestplot_favouring_function.R")

lung <- survival::lung

lung2 <- lung %>%
  mutate(
    event_os = status - 1,
    time_os  = time,
    sex      = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    ph.ecog  = factor(ph.ecog),
    treatment = sample(c("Placebo", "Experimental"), n(), replace = TRUE),
    time_pfs  = pmin(time_os, time_os * runif(n(), 0.4, 0.9)),
    event_pfs = ifelse(time_pfs < time_os, 1, event_os)
  ) %>%
  filter(!is.na(ph.ecog), !is.na(sex), !is.na(age))

ps.formula <- as.formula(paste("treatment", "~", paste(c("sex", "age", "ph.ecog"), collapse = "+")))


W.out <- weightit(ps.formula, data = lung2, method = "ps", estimand = "ATE")

lung2$weights_ATE <- W.out$weights

res2 <- cox_subgroup_forest(
  data = lung2,
  covariates = c("sex", "ph.ecog"),
  labels = c("Sex", "ECOG PS"),
  treatment = "treatment",
  time = c("time_os", "time_pfs"),
  event = c("event_os", "event_pfs"),
  endpoint_labels = c("Overall Survival", "Progression-Free Survival"),
  time_unit = "days",
  palette = c("darkorange", "darkgreen"),
  xlim = c(0.25, 4),
  ticks_at = c(0.5, 1, 2),
  length_plot = 60,
  weights = "weights_ATE", # optional if PS used
  ci_pch = 15, 
  ci_lty = 1,
  ci_lwd = 3, 
  ci_Theight = 0.4, 
  sizes = 1
)

res2$combined_plots


