library(mixedcirc)

expression_data <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 2,
  groups = c("A", "B"),
  seed = 501
)

expression_fit <- mixedcirc_detect(
  expression_data$data_input,
  expression_data$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lm",
  random_effect_type = "none"
)

stopifnot(
  all(vapply(
    expression_fit@results,
    function(x) identical(x@type, "expression"),
    logical(1)
  )),
  all(vapply(
    expression_fit@results,
    function(x) identical(x@params$RRBS, FALSE),
    logical(1)
  ))
)

rrbs_data <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 2,
  groups = c("A", "B"),
  rrbs = TRUE,
  seed = 502
)

rrbs_fit <- mixedcirc_detect(
  rrbs_data$data_input,
  rrbs_data$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lm",
  random_effect_type = "none",
  obs_weights = rrbs_data$obs_weights,
  RRBS = TRUE,
  replicate_id = "replicate_id"
)

stopifnot(
  all(vapply(
    rrbs_fit@results,
    function(x) identical(x@type, "RRBS"),
    logical(1)
  )),
  all(vapply(
    rrbs_fit@results,
    function(x) identical(x@params$RRBS, TRUE),
    logical(1)
  ))
)
