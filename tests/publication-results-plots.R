library(mixedcirc)

# Three-group inference: the omnibus differential test must be distinct from
# the omnibus test that any group is rhythmic.
set.seed(1901)
time <- rep(seq(0, 21, by = 3), 3)
group <- factor(rep(c("A", "B", "C"), each = 8), levels = c("A", "B", "C"))
meta <- data.frame(time = time, group = group)
signal <- c(
  5 + 2 * cos(2 * pi * time[group == "A"] / 24),
  5 + 2 * cos(2 * pi * time[group == "B"] / 24),
  5 + 2 * cos(2 * pi * time[group == "C"] / 24 + pi / 2)
) + stats::rnorm(length(time), sd = 0.15)
null <- 4 + stats::rnorm(length(time), sd = 0.8)
fit_three <- suppressWarnings(mixedcirc_detect(
  cbind(phase_shift = signal, null = null),
  meta,
  time = "time",
  group = "group",
  id = NULL,
  lm_method = "lm"
))

wide <- mixedcirc_results(fit_three, "wide")
estimates <- mixedcirc_results(fit_three, "estimates")
contrasts <- mixedcirc_contrasts(fit_three)
classes <- mixedcirc_classify(fit_three)
diagnostics <- mixedcirc_results(fit_three, "diagnostics")
analysis_summary <- mixedcirc_summary(fit_three)

stopifnot(
  nrow(wide) == 2L,
  nrow(estimates) == 6L,
  nrow(contrasts) == 6L,
  nrow(classes) == 6L,
  nrow(diagnostics) == 2L,
  all(c(
    "beta_cos", "beta_sin", "peak_time", "q_rhythm",
    "q_differential_global"
  ) %in% names(estimates)),
  all(c(
    "harmonic_distance", "delta_peak_time", "q_differential"
  ) %in% names(contrasts)),
  is.finite(wide$p.value_differential_global[1L]),
  wide$p.value_differential_global[1L] < 0.01,
  all(abs(contrasts$delta_phase[is.finite(contrasts$delta_phase)]) <= pi),
  all(abs(contrasts$delta_peak_time[is.finite(contrasts$delta_peak_time)]) <= 12),
  identical(as.data.frame(fit_three), wide),
  all(diagnostics$fit_status == "ok"),
  all(diagnostics$converged),
  inherits(analysis_summary, "mixedcirc_summary"),
  analysis_summary$overview$features == 2L,
  analysis_summary$overview$groups == 3L,
  identical(summary(fit_three), analysis_summary)
)

for (group_name in unique(estimates$group)) {
  current <- estimates[estimates$group == group_name, , drop = FALSE]
  rhythmic <- is.finite(current$q_rhythm) & current$q_rhythm < 0.05
  reported <- analysis_summary$groups$median_amplitude[
    analysis_summary$groups$group == group_name
  ]
  expected <- if (any(rhythmic)) {
    stats::median(current$amplitude[rhythmic])
  } else {
    NA_real_
  }
  stopifnot(isTRUE(all.equal(reported, expected)))
}

# With one feature failing, on_error = "keep" must preserve successful fits and
# make the failure inspectable. Default fail-fast behavior is tested elsewhere.
bad_data <- cbind(good = signal, all_missing = rep(NA_real_, length(signal)))
kept <- suppressWarnings(mixedcirc_detect(
  bad_data,
  meta,
  time = "time",
  group = "group",
  id = NULL,
  lm_method = "lm",
  on_error = "keep"
))
kept_diagnostics <- mixedcirc_results(kept, "diagnostics")
stopifnot(
  identical(kept_diagnostics$fit_status, c("ok", "failed")),
  is.character(kept_diagnostics$fit_message[2L]),
  nzchar(kept_diagnostics$fit_message[2L])
)

# Compact results keep statistics and diagnostics while dropping heavy models.
compact <- suppressWarnings(mixedcirc_detect(
  cbind(phase_shift = signal, null = null),
  meta,
  time = "time",
  group = "group",
  id = NULL,
  lm_method = "lm",
  store_models = FALSE
))
compact_diagnostics <- mixedcirc_results(compact, "diagnostics")
stopifnot(
  all(vapply(compact@results, function(object) is.null(object@fit), logical(1))),
  all(compact_diagnostics$model_class == "lm"),
  all(compact_diagnostics$n_obs == nrow(meta)),
  all(compact_diagnostics$converged),
  all(is.finite(compact_diagnostics$AIC)),
  as.numeric(utils::object.size(compact)) < as.numeric(utils::object.size(fit_three))
)
compact_interpolation_error <- tryCatch(
  mixedcirc_interpolate(compact[[1L]]),
  error = conditionMessage
)
compact_bootstrap_error <- tryCatch(
  mixedcirc_boottest(compact[[1L]], nsim = 2, seed = 9),
  error = conditionMessage
)
stopifnot(
  grepl("store_models = TRUE", compact_interpolation_error, fixed = TRUE),
  grepl("store_models = TRUE", compact_bootstrap_error, fixed = TRUE)
)

# Dataset-level plots and subject estimates on a repeated-measures fit.
simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 4,
  n_timepoints = 8,
  n_features = 5,
  groups = c("control", "treated"),
  residual_sd = 0.45,
  seed = 1902
)
mixed_fit <- suppressWarnings(mixedcirc_detect(
  simulation$data_input,
  simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lme",
  random_effect_type = "mesor_phase"
))
subject_estimates <- mixedcirc_subject_rhythms(
  mixed_fit,
  features = colnames(simulation$data_input)[1:2]
)
stopifnot(
  nrow(subject_estimates) == 16L,
  all(subject_estimates$subject_peak_time >= 0),
  all(subject_estimates$subject_peak_time < 24)
)

plots <- list(
  mixedcirc_plot_heatmap(mixed_fit, npoints = 24),
  mixedcirc_plot_differential(mixed_fit, label_n = 2),
  mixedcirc_plot_polar(mixed_fit, label_n = 2),
  mixedcirc_plot_design(mixed_fit),
  mixedcirc_plot_diagnostics(
    mixed_fit,
    features = colnames(simulation$data_input)[1:2]
  ),
  mixedcirc_plot_diagnostics(
    mixed_fit,
    features = colnames(simulation$data_input)[1:2],
    type = "qq"
  ),
  mixedcirc_plot_subjects(
    mixed_fit,
    features = colnames(simulation$data_input)[1:2],
    npoints = 24
  )
)
stopifnot(all(vapply(plots, inherits, logical(1), what = "ggplot")))
stopifnot(all(vapply(
  plots,
  function(plot) length(suppressMessages(ggplot2::ggplot_build(plot))$data) > 0L,
  logical(1)
)))
polar_alpha_error <- tryCatch(
  mixedcirc_plot_polar(mixed_fit, alpha = 2),
  error = conditionMessage
)
stopifnot(grepl("between 0 and 1", polar_alpha_error, fixed = TRUE))
polar_plot <- mixedcirc_plot_polar(mixed_fit, label_n = 2)
stopifnot(
  identical(polar_plot$labels$alpha, "Rhythmic"),
  identical(levels(polar_plot$data$rhythmic), c("Yes", "No"))
)

bootstrap <- suppressWarnings(mixedcirc_boottest(
  mixed_fit[[1L]], nsim = 3, seed = 1903, parallel = "no"
))
bootstrap_plot <- mixedcirc_plot_bootstrap(bootstrap)
stopifnot(
  inherits(bootstrap_plot, "ggplot"),
  length(ggplot2::ggplot_build(bootstrap_plot)$data) > 0L
)

circular_boot <- methods::new(
  "mixedcirc_boot",
  results = data.frame(
    boot.lower = 3.05,
    boot.upper = -3.05,
    boot.estimate = 3.12,
    boot.SE = 0.08,
    pvalue = 0.1,
    row.names = "phase_A"
  ),
  metadata = list(conf_level = 0.95),
  original_model = mixed_fit[[1L]]
)
circular_plot_data <- ggplot2::ggplot_build(
  mixedcirc_plot_bootstrap(circular_boot)
)$data[[2L]]
stopifnot(abs(circular_plot_data$xend - circular_plot_data$x) < 0.3)

duplicate_match_error <- tryCatch(
  mixedcirc_compare(
    mixed_fit[[1L]],
    mixed_fit[[2L]],
    lags = 0,
    analysis = "pearson",
    match_by = "time"
  ),
  error = conditionMessage
)
stopifnot(grepl("uniquely identify", duplicate_match_error, fixed = TRUE))
