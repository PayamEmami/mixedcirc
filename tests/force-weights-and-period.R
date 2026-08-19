library(mixedcirc)

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 8,
  n_features = 3,
  groups = c("A", "B"),
  period = 12,
  seed = 821
)

set.seed(822)
counts <- matrix(
  stats::rpois(nrow(sim$meta_input) * 3L, lambda = 30),
  nrow = nrow(sim$meta_input),
  dimnames = list(NULL, paste0("count_", 1:3))
)

weighted_fit <- suppressWarnings(mixedcirc_detect(
  counts,
  sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  period = 12,
  lm_method = "lm",
  random_effect_type = "none",
  force_weight_estimation = TRUE
))
weighted_removed <- mixedcirc_remove_trend(
  counts,
  sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lm",
  random_effect_type = "none",
  force_weight_estimation = TRUE
)

implicit_period <- mixedcirc_interpolate(weighted_fit[[1L]], npoints = 31)
explicit_period <- mixedcirc_interpolate(weighted_fit[[1L]], period = 12, npoints = 31)
plot_object <- mixedcirc_fit_plot(weighted_fit[[1L]], n_points = 31)

stopifnot(
  identical(names(weighted_fit@results), colnames(counts)),
  all(vapply(weighted_fit@results, function(z) identical(z@params$period, 12), logical(1))),
  all(vapply(weighted_fit@results, function(z) all(is.finite(z@exp_design$.weight)), logical(1))),
  identical(dim(weighted_removed), dim(counts)),
  all(is.finite(weighted_removed)),
  isTRUE(all.equal(implicit_period, explicit_period, tolerance = 1e-12)),
  inherits(plot_object, "ggplot"),
  length(ggplot2::ggplot_build(plot_object)$data) > 0L
)

phase_columns <- grep("^phase_", names(weighted_fit[[1L]]@results), value = TRUE)
peak_times <- mixedcirc_detectpeaks(
  unlist(weighted_fit[[1L]]@results[phase_columns], use.names = FALSE),
  period = weighted_fit[[1L]]@params$period
)
stopifnot(all(peak_times >= 0 & peak_times < 12))

# The simulator, fitter, stored period, phase convention, and peak conversion
# must agree for a non-24-hour rhythm.
exact_sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 2,
  n_timepoints = 12,
  n_features = 2,
  groups = c("A", "B"),
  period = 12,
  subject_random_intercept_sd = 0,
  residual_sd = 0,
  seed = 823
)
exact_fit <- suppressWarnings(mixedcirc_detect(
  exact_sim$data_input,
  exact_sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  period = 12,
  lm_method = "lm",
  random_effect_type = "none"
))
circular_difference <- function(x, y, cycle) {
  ((x - y + cycle / 2) %% cycle) - cycle / 2
}
for (i in seq_along(exact_fit)) {
  estimated_phase <- unlist(
    exact_fit[[i]]@results[c("phase_A", "phase_B")],
    use.names = FALSE
  )
  true_phase <- unlist(
    exact_sim$truth[i, c("phase_A", "phase_B")],
    use.names = FALSE
  )
  stopifnot(
    max(abs(circular_difference(estimated_phase, true_phase, 2 * pi))) < 1e-10,
    all(mixedcirc_detectpeaks(estimated_phase, period = 12) >= 0),
    all(mixedcirc_detectpeaks(estimated_phase, period = 12) < 12)
  )
}

copy_sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 2,
  n_timepoints = 5,
  n_features = 2,
  groups = "A",
  cross_lag_spec = list(list(from = 1, to = 2, lag = 0, beta = 1)),
  seed = 824
)
stopifnot(identical(copy_sim$data_input[[1L]], copy_sim$data_input[[2L]]))

# Nonzero simulated lags must rotate both the reported phase and harmonic
# coefficients by one sampling step, not by one raw time unit.
lag_truth <- mixedcirc_simulate_data(
  n_subjects_per_group = 2,
  n_timepoints = 8,
  n_features = 2,
  groups = c("A", "B"),
  cross_lag_spec = list(
    list(from = 1, to = 2, lag = c(A = 1, B = 2), beta = -0.7)
  ),
  period = 12,
  seed = 825
)

signed_phase <- function(value) {
  wrapped <- value %% (2 * pi)
  ifelse(abs(wrapped) < 1e-12, -2 * pi, wrapped - 2 * pi)
}

for (group_name in c("A", "B")) {
  source_phase <- lag_truth$truth[[paste0("phase_", group_name)]][1L]
  target_phase <- lag_truth$truth[[paste0("phase_", group_name)]][2L]
  lag_steps <- lag_truth$truth[[paste0("copied_lag_", group_name)]][2L]
  expected_phase <- signed_phase(
    source_phase + pi - 2 * pi * lag_steps / 8
  )
  target_amplitude <- lag_truth$truth[[paste0("amplitude_", group_name)]][2L]
  target_cosine <- lag_truth$truth[[paste0("beta_cos_", group_name)]][2L]
  target_sine <- lag_truth$truth[[paste0("beta_sin_", group_name)]][2L]

  stopifnot(
    abs(circular_difference(target_phase, expected_phase, 2 * pi)) < 1e-12,
    abs(target_cosine - target_amplitude * cos(target_phase)) < 1e-12,
    abs(target_sine + target_amplitude * sin(target_phase)) < 1e-12,
    target_phase >= -2 * pi,
    target_phase < 0
  )
}
