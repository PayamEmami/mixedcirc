library(mixedcirc)

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 3,
  groups = c("A", "B"),
  rrbs = TRUE,
  cross_lag_spec = list(
    list(from = 1, to = 2, lag = 1, beta = 1)
  ),
  seed = 811
)

fit <- suppressWarnings(mixedcirc_detect(
  sim$data_input,
  sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = sim$obs_weights,
  lm_method = "lm",
  random_effect_type = "none"
))

biological_n <- nrow(sim$meta_input)
interpolated <- mixedcirc_interpolate(fit, npoints = 19)
plot_object <- mixedcirc_fit_plot(fit[[1L]], n_points = 19)
lagged <- mixedcirc_add_lag(fit[[1L]], lags = c(-1L, 1L))
comparison <- mixedcirc_compare(
  fit[[1L]],
  fit[[2L]],
  lags = c(-1L, 0L, 1L),
  analysis = "pearson"
)

stopifnot(
  length(fit) == 3L,
  all(vapply(fit@results, function(z) identical(z@type, "RRBS"), logical(1))),
  identical(dim(mixedcirc_fitted(fit)), c(2L * biological_n, 3L)),
  nrow(interpolated) == 3L * 2L * 19L,
  inherits(plot_object, "ggplot"),
  length(ggplot2::ggplot_build(plot_object)$data) > 0L,
  nrow(lagged) == 2L * biological_n,
  sum(!is.na(lagged$measure_lag1)) == 2L * (biological_n - 6L),
  nrow(comparison$results) == 6L,
  all(comparison$results$n.used <= biological_n / 2L),
  any(is.finite(comparison$results$estimate))
)

removed <- mixedcirc_remove_trend(
  sim$data_input,
  sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = sim$obs_weights,
  lm_method = "lm",
  random_effect_type = "none"
)
stopifnot(identical(dim(removed), dim(sim$data_input)))

