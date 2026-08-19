library(mixedcirc)

expect_error <- function(expr, pattern) {
  message <- tryCatch({
    force(expr)
    NULL
  }, error = conditionMessage)
  stopifnot(is.character(message), length(message) == 1L, grepl(pattern, message))
}

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 2,
  n_timepoints = 6,
  n_features = 2,
  groups = c("A", "B"),
  seed = 841
)

expect_error(
  mixedcirc_detect(sim$data_input, sim$meta_input, time = "time", period = 0),
  "period"
)
expect_error(
  mixedcirc_detect(sim$data_input, sim$meta_input, time = "time", ncores = 0),
  "ncores"
)
expect_error(
  mixedcirc_detect(sim$data_input, sim$meta_input, time = "time", verbose = NA),
  "verbose"
)
expect_error(
  mixedcirc_detect(
    transform(sim$data_input, bad = letters[seq_len(nrow(sim$data_input))]),
    sim$meta_input,
    time = "time"
  ),
  "numeric"
)
expect_error(
  mixedcirc_detect(
    sim$data_input,
    transform(sim$meta_input, inphase = 1),
    time = "time"
  ),
  "reserved"
)

switched <- suppressWarnings(mixedcirc_detect(
  sim$data_input[, 1L, drop = FALSE],
  sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lme",
  random_effect_type = "none"
))
stopifnot(inherits(switched[[1L]]@fit, "lm"))

expect_error(
  mixedcirc_boottest(switched[[1L]], nsim = 1.5, seed = 1),
  "positive integer"
)
expect_error(
  mixedcirc_boottest(switched[[1L]], nsim = 2, seed = Inf),
  "finite integer"
)
expect_error(
  mixedcirc_add_lag(switched[[1L]], lags = 1.5),
  "integers"
)
expect_error(
  mixedcirc_interpolate(methods::new("mixedcirc_fit_list", results = list())),
  "at least one"
)
expect_error(
  mixedcirc_integrate(
    list(one = as.matrix(sim$data_input)),
    sim$meta_input,
    time = "time"
  ),
  "at least two"
)

