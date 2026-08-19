library(mixedcirc)

simulated <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 2,
  groups = c("A", "B"),
  seed = 201
)

lm_result <- mixedcirc_remove_trend(
  data_input = simulated$data_input,
  meta_input = simulated$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lm",
  random_effect_type = "none",
  ncores = 1
)

stopifnot(
  identical(dim(lm_result), dim(simulated$data_input)),
  all(is.finite(lm_result)),
  inherits(future::plan("list")[[1L]], "sequential")
)

mixed_result <- mixedcirc_remove_trend(
  data_input = simulated$data_input,
  meta_input = simulated$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lme",
  random_effect_type = "mesor",
  ncores = 2
)

stopifnot(
  identical(dim(mixed_result), dim(simulated$data_input)),
  all(is.finite(mixed_result)),
  inherits(future::plan("list")[[1L]], "sequential")
)

input_error <- tryCatch(
  mixedcirc_remove_trend(
    data_input = NULL,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    ncores = 2
  ),
  error = conditionMessage
)

stopifnot(
  identical(input_error, "data_input must be a data frame, matrix, or DGEList"),
  inherits(future::plan("list")[[1L]], "sequential")
)
