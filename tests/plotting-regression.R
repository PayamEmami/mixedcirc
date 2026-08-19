library(mixedcirc)

plot_has_data <- function(object) {
  inherits(object, "ggplot") &&
    length(suppressMessages(ggplot2::ggplot_build(object))$data) > 0L
}

# The S4 plot method must continue to dispatch to mixedcirc_fit_plot(), retain
# the fitted period, and draw the same population curve as interpolation.
single_simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 8,
  n_features = 2,
  groups = "A",
  period = 12,
  seed = 2201
)
single_fit <- suppressWarnings(mixedcirc_detect(
  single_simulation$data_input,
  single_simulation$meta_input,
  time = "time",
  id = NULL,
  period = 12,
  lm_method = "lm",
  random_effect_type = "none"
))
direct_plot <- mixedcirc_fit_plot(single_fit[[1L]], n_points = 31)
dispatch_plot <- plot(single_fit[[1L]], n_points = 31)
interpolation <- mixedcirc_interpolate(single_fit[[1L]], npoints = 31)
direct_build <- suppressMessages(ggplot2::ggplot_build(direct_plot))
dispatch_build <- suppressMessages(ggplot2::ggplot_build(dispatch_plot))
direct_line <- direct_build$data[[length(direct_build$data)]]
dispatch_line <- dispatch_build$data[[length(dispatch_build$data)]]

stopifnot(
  plot_has_data(direct_plot),
  plot_has_data(dispatch_plot),
  identical(single_fit[[1L]]@params$period, 12),
  nrow(direct_line) == 31L,
  nrow(dispatch_line) == 31L,
  isTRUE(all.equal(direct_line$x, interpolation$time, tolerance = 1e-12)),
  isTRUE(all.equal(direct_line$y, interpolation$Y.hat, tolerance = 1e-10)),
  isTRUE(all.equal(dispatch_line$x, direct_line$x, tolerance = 1e-12)),
  isTRUE(all.equal(dispatch_line$y, direct_line$y, tolerance = 1e-10))
)

# Grouped mixed models with nuisance covariates exercise the factor averaging,
# fixed-value, linetype, facet, and label branches of the original plot helper.
grouped_simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 2,
  groups = c("control", "treated"),
  include_covariates = TRUE,
  seed = 2202
)
grouped_fit <- suppressWarnings(mixedcirc_detect(
  grouped_simulation$data_input,
  grouped_simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  formula_extra = ~ age + batch,
  lm_method = "lme",
  random_effect_type = "mesor"
))
expanded_plots <- lapply(
  c("linetype", "facet", "label"),
  function(display) mixedcirc_fit_plot(
    grouped_fit[[1L]],
    n_points = 19,
    fixed_values = list(age = 50),
    expand_vars = "batch",
    show_expanded = display
  )
)
stopifnot(all(vapply(expanded_plots, plot_has_data, logical(1))))

# The RRBS plotting branch must pair methylated/unmethylated observations even
# after fitted-value, residual, and weight audit columns were added.
rrbs_simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 2,
  groups = c("A", "B"),
  rrbs = TRUE,
  seed = 2203
)
rrbs_fit <- suppressWarnings(mixedcirc_detect(
  rrbs_simulation$data_input,
  rrbs_simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = rrbs_simulation$obs_weights,
  lm_method = "lm",
  random_effect_type = "none"
))
rrbs_plot <- plot(rrbs_fit[[1L]], n_points = 23)
rrbs_interpolation <- mixedcirc_interpolate(rrbs_fit[[1L]], npoints = 23)
stopifnot(
  plot_has_data(rrbs_plot),
  nrow(rrbs_interpolation) == 2L * 23L,
  !"scaler" %in% names(rrbs_interpolation),
  all(is.finite(rrbs_interpolation$Y.hat))
)

# Compact results deliberately cannot recreate model-based plots, but the
# error must be explicit rather than a broken dispatch or obscure prediction
# failure.
compact_fit <- suppressWarnings(mixedcirc_detect(
  single_simulation$data_input,
  single_simulation$meta_input,
  time = "time",
  id = NULL,
  period = 12,
  lm_method = "lm",
  random_effect_type = "none",
  store_models = FALSE
))
compact_plot_error <- tryCatch(
  plot(compact_fit[[1L]]),
  error = conditionMessage
)
stopifnot(grepl("store_models = TRUE", compact_plot_error, fixed = TRUE))
