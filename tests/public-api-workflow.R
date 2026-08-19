library(mixedcirc)

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 4,
  n_timepoints = 6,
  n_features = 4,
  groups = c("A", "B"),
  include_covariates = TRUE,
  seed = 801
)

fit <- suppressWarnings(mixedcirc_detect(
  data_input = sim$data_input,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  formula_extra = ~ age + batch,
  lm_method = "lme",
  random_effect_type = "mesor",
  ncores = 1
))

stopifnot(
  methods::is(fit, "mixedcirc_fit_list"),
  length(fit) == 4L,
  identical(names(fit@results), colnames(sim$data_input)),
  all(vapply(fit@results, function(z) inherits(z@fit, "merMod"), logical(1))),
  all(vapply(fit@results, function(z) identical(z@params$period, 24), logical(1)))
)

fitted_values <- mixedcirc_fitted(fit)
detrended_values <- suppressWarnings(mixedcirc_detrend(fit, per_group = TRUE))
interpolated <- mixedcirc_interpolate(fit, npoints = 17)
plot_object <- mixedcirc_fit_plot(
  fit[[1L]],
  n_points = 17,
  expand_vars = "batch",
  show_expanded = "linetype"
)

stopifnot(
  identical(dim(fitted_values), c(nrow(sim$data_input), 4L)),
  identical(dim(detrended_values), c(nrow(sim$data_input), 4L)),
  all(is.finite(fitted_values)),
  all(is.finite(detrended_values)),
  nrow(interpolated) == 4L * 2L * 17L,
  all(c("Y.hat", "feature", ".expanded_label") %in% names(interpolated)),
  inherits(plot_object, "ggplot"),
  length(ggplot2::ggplot_build(plot_object)$data) > 0L
)

adjusted <- mixedcirc_adjust(fit, method = "BH", pattern = "^p.value_global$")
stopifnot(all(vapply(
  adjusted@results,
  function(z) "p.value_global_adjusted_BH" %in% names(z@results),
  logical(1)
)))

lagged <- mixedcirc_add_lag(
  fit[[1L]],
  lags = c(-1L, 0L, 1L),
  type = "fitted"
)
comparison <- suppressWarnings(mixedcirc_compare(
  fit[[1L]],
  fit[[2L]],
  lags = c(-1L, 0L, 1L),
  analysis = "association",
  type = "fitted",
  return_models = TRUE,
  return_data = TRUE
))

stopifnot(
  nrow(lagged) == nrow(sim$data_input),
  all(c("measure_lead1", "measure_lag0", "measure_lag1") %in% names(lagged)),
  nrow(comparison$results) == 6L,
  all(comparison$results$method %in% c("lmer", "lm")),
  length(comparison$models) == 6L,
  is.data.frame(comparison$data)
)

pearson_comparison <- mixedcirc_compare(
  fit[[1L]],
  fit[[2L]],
  lags = c(-1L, 0L, 1L),
  analysis = "pearson",
  type = "fitted"
)
stopifnot(
  nrow(pearson_comparison$results) == 6L,
  all(pearson_comparison$results$analysis == "pearson")
)

if (requireNamespace("rmcorr", quietly = TRUE)) {
  repeated_comparison <- mixedcirc_compare(
    fit[[1L]],
    fit[[2L]],
    lags = c(-1L, 0L, 1L),
    analysis = "rmcorr",
    type = "fitted"
  )
  stopifnot(
    nrow(repeated_comparison$results) == 6L,
    all(repeated_comparison$results$analysis == "rmcorr")
  )
}

mixed_bootstrap <- suppressWarnings(mixedcirc_boottest(
  fit[[1L]],
  nsim = 4,
  seed = 802,
  parallel = "no"
))
stopifnot(
  methods::is(mixed_bootstrap, "mixedcirc_boot"),
  identical(mixed_bootstrap@metadata$model_class, "lme4"),
  nrow(mixed_bootstrap@results) > 0L,
  any(is.finite(mixed_bootstrap@results$boot.estimate))
)

# Missing responses and nuisance covariates must remain aligned with the
# original observation rows rather than shifting fitted values.
missing_sim <- sim
missing_sim$data_input[2L, 1L] <- NA_real_
missing_sim$meta_input$age[3L] <- NA_real_
missing_fit <- suppressWarnings(mixedcirc_detect(
  missing_sim$data_input[, 1:2, drop = FALSE],
  missing_sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  formula_extra = ~ age,
  lm_method = "lm",
  random_effect_type = "none"
))
missing_fitted <- mixedcirc_fitted(missing_fit)
stopifnot(
  identical(dim(missing_fitted), c(nrow(missing_sim$data_input), 2L)),
  is.na(missing_fitted[2L, 1L]),
  all(is.na(missing_fitted[3L, ])),
  all(is.finite(missing_fitted[-c(2L, 3L), 1L]))
)
