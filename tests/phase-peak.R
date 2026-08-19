library(mixedcirc)

circular_difference <- function(x, y, cycle) {
  ((x - y + cycle / 2) %% cycle) - cycle / 2
}

# Phase is now always the signed model shift. The old absolute-phase switch is
# intentionally absent from both public entry points.
stopifnot(
  !"abs_phase" %in% names(formals(mixedcirc_detect)),
  !"abs_phase" %in% names(formals(mixedcirc_boottest))
)

# Direct conversion covers every quadrant, midnight, values outside one cycle,
# vector inputs, missing phases, unwrapped output, and the minimum-time option.
signed_phases <- c(-2 * pi, -pi / 4, -pi / 2, -pi, -3 * pi / 2, -7 * pi / 4)
expected_peaks <- c(0, 3, 6, 12, 18, 21)

stopifnot(
  isTRUE(all.equal(
    mixedcirc_detectpeaks(signed_phases, period = 24),
    expected_peaks,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    mixedcirc_detectpeaks(c(-5 * pi / 2, NA_real_), period = 24),
    c(6, NA_real_),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    mixedcirc_detectpeaks(-5 * pi / 2, period = 24, correct_period = FALSE),
    30,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    mixedcirc_detectpeaks(
      c(-2 * pi, -pi / 2),
      period = 24,
      min_time = 2,
      correct_min = TRUE
    ),
    c(24, 6),
    tolerance = 1e-12
  ))
)

invalid_period_error <- tryCatch(
  mixedcirc_detectpeaks(-pi / 2, period = 0),
  error = conditionMessage
)
stopifnot(identical(
  invalid_period_error,
  "`period` must be a single positive finite number."
))

# Exact rhythms with known peak times verify the full path from the fitted
# cosine/sine coefficients through calc_phase() and mixedcirc_detectpeaks().
period <- 24
times <- 0:23
true_peaks <- c(0, 3, 6, 12, 18, 21)
exact_data <- vapply(
  true_peaks,
  function(peak) 10 + 2 * cos(2 * pi * (times - peak) / period),
  numeric(length(times))
)
colnames(exact_data) <- paste0("peak_", true_peaks)

exact_fit <- suppressWarnings(mixedcirc_detect(
  data_input = exact_data,
  meta_input = data.frame(time = times),
  time = "time",
  period = period,
  lm_method = "lm",
  random_effect_type = "none",
  ncores = 1
))

detected_phases <- vapply(
  exact_fit@results,
  function(x) x@results$phase[1],
  numeric(1)
)
expected_phases <- -2 * pi * true_peaks / period
detected_peaks <- mixedcirc_detectpeaks(detected_phases, period = period)

stopifnot(
  max(abs(circular_difference(detected_phases, expected_phases, 2 * pi))) < 1e-10,
  max(abs(circular_difference(detected_peaks, true_peaks, period))) < 1e-10,
  all(detected_phases <= 1e-10)
)

# Repeat the end-to-end check with noisy replicated observations. This checks
# estimation, not only an algebraically exact conversion.
set.seed(710)
noisy_times <- rep(0:23, times = 12)
noisy_true_peaks <- c(2, 7.5, 13, 20)
noisy_data <- vapply(
  noisy_true_peaks,
  function(peak) {
    7 + 2.5 * cos(2 * pi * (noisy_times - peak) / period) +
      rnorm(length(noisy_times), sd = 0.3)
  },
  numeric(length(noisy_times))
)
colnames(noisy_data) <- paste0("noisy_peak_", noisy_true_peaks)

noisy_fit <- suppressWarnings(mixedcirc_detect(
  data_input = noisy_data,
  meta_input = data.frame(time = noisy_times),
  time = "time",
  period = period,
  lm_method = "lm",
  random_effect_type = "none",
  ncores = 1
))

noisy_phases <- vapply(
  noisy_fit@results,
  function(x) x@results$phase[1],
  numeric(1)
)
noisy_peaks <- mixedcirc_detectpeaks(noisy_phases, period = period)

stopifnot(
  max(abs(circular_difference(
    noisy_phases,
    -2 * pi * noisy_true_peaks / period,
    2 * pi
  ))) < 0.08,
  max(abs(circular_difference(noisy_peaks, noisy_true_peaks, period))) < 0.3
)

# The package simulator's truth and grouped detection must use the same signed
# phase convention. Zero noise makes this a precise implementation check.
simulated <- mixedcirc_simulate_data(
  n_subjects_per_group = 2,
  n_timepoints = 12,
  n_features = 3,
  groups = c("A", "B"),
  subject_random_intercept_sd = 0,
  residual_sd = 0,
  seed = 711
)

simulated_fit <- mixedcirc_detect(
  data_input = simulated$data_input,
  meta_input = simulated$meta_input,
  time = "time",
  group = "group",
  id = "id",
  period = period,
  lm_method = "lm",
  random_effect_type = "none",
  ncores = 1
)

for (feature_i in seq_along(simulated_fit@results)) {
  detected <- simulated_fit@results[[feature_i]]@results
  truth <- simulated$truth[feature_i, ]

  fitted_phase <- c(detected$phase_A, detected$phase_B)
  true_phase <- c(truth$phase_A, truth$phase_B)
  fitted_peak <- mixedcirc_detectpeaks(fitted_phase, period = period)
  true_peak <- (
    atan2(
      c(truth$beta_sin_A, truth$beta_sin_B),
      c(truth$beta_cos_A, truth$beta_cos_B)
    ) %% (2 * pi)
  ) * period / (2 * pi)

  stopifnot(
    max(abs(circular_difference(fitted_phase, true_phase, 2 * pi))) < 1e-10,
    max(abs(circular_difference(fitted_peak, true_peak, period))) < 1e-10
  )
}

grouped_boot_stats <- mixedcirc:::.mc_extract_stats_from_fit(
  fit = simulated_fit[[1]]@fit,
  object = simulated_fit[[1]],
  wrap_phase = TRUE
)
grouped_detected <- simulated_fit[[1]]@results
expected_phase_difference <- circular_difference(
  grouped_detected$phase_A,
  grouped_detected$phase_B,
  2 * pi
)
expected_peak_difference <- circular_difference(
  mixedcirc_detectpeaks(grouped_detected$phase_A, period),
  mixedcirc_detectpeaks(grouped_detected$phase_B, period),
  period
)
phase_difference_as_peak_hours <- circular_difference(
  -expected_phase_difference * period / (2 * pi),
  0,
  period
)

stopifnot(
  isTRUE(all.equal(
    unname(grouped_boot_stats[c("phase_A", "phase_B")]),
    unname(c(grouped_detected$phase_A, grouped_detected$phase_B)),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    unname(grouped_boot_stats["phase_AvsB"]),
    unname(expected_phase_difference),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    unname(phase_difference_as_peak_hours),
    unname(expected_peak_difference),
    tolerance = 1e-10
  ))
)

# Bootstrap point estimates are recomputed from the fitted coefficients and
# must retain the signed phase rather than silently taking its absolute value.
bootstrap <- mixedcirc_boottest(
  noisy_fit[[2]],
  nsim = 8,
  seed = 712,
  parallel = "no"
)

stopifnot(
  isTRUE(all.equal(
    unname(bootstrap@results["phase", "boot.estimate"]),
    unname(noisy_phases[2]),
    tolerance = 1e-10
  )),
  bootstrap@results["phase", "boot.estimate"] < 0,
  !"abs_phase" %in% names(bootstrap@metadata)
)
