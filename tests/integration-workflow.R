library(mixedcirc)

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 7,
  groups = c("A", "B"),
  seed = 831
)
blocks <- list(
  transcript = as.matrix(sim$data_input[, 1:3]),
  protein = as.matrix(sim$data_input[, 4:7])
)

fit <- mixedcirc_integrate(
  data_input = blocks,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  ncomp = 2,
  variables = list(transcript = c(2, 2), protein = c(3, 3)),
  max.iter = 2000
)

stopifnot(
  methods::is(fit, "mixedcirc_integration"),
  identical(dim(mixedcirc_getscore(fit, "average")), c(nrow(sim$meta_input), 2L)),
  identical(dim(mixedcirc_getscore(fit, "partial", merge = TRUE)), c(nrow(sim$meta_input), 4L)),
  identical(dim(mixedcirc_getscore(fit, "loadings", merge = TRUE)), c(7L, 2L)),
  identical(dim(mixedcirc_getscore(fit, "partial", "transcript")), c(nrow(sim$meta_input), 2L)),
  is.list(mixedcirc_getscore(fit, "model"))
)

merged_fit <- mixedcirc_integrate(
  data_input = blocks,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  ncomp = 1,
  merge = TRUE,
  no_interaction = TRUE,
  center = TRUE,
  scale = TRUE,
  max.iter = 2000
)
decomposed_fit <- mixedcirc_integrate(
  data_input = blocks,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  ncomp = 1,
  decompose = TRUE,
  max.iter = 2000
)
uncorrelated_fit <- mixedcirc_integrate(
  data_input = blocks,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  ncomp = 1,
  no_correlation = TRUE,
  max.iter = 2000
)

stopifnot(
  methods::is(merged_fit, "mixedcirc_integration"),
  methods::is(decomposed_fit, "mixedcirc_integration"),
  methods::is(uncorrelated_fit, "mixedcirc_integration"),
  ncol(mixedcirc_getscore(merged_fit, "average")) == 1L,
  ncol(mixedcirc_getscore(decomposed_fit, "average")) == 1L,
  ncol(mixedcirc_getscore(uncorrelated_fit, "average")) == 1L
)
