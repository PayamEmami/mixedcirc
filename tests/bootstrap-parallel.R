library(mixedcirc)

simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 3,
  n_timepoints = 6,
  n_features = 1,
  groups = c("A", "B"),
  seed = 941
)

ordinary_fit <- suppressWarnings(mixedcirc_detect(
  simulation$data_input,
  simulation$meta_input,
  time = "time",
  group = "group",
  id = NULL,
  lm_method = "lm",
  random_effect_type = "none"
))

serial <- suppressMessages(suppressWarnings(mixedcirc_boottest(
  ordinary_fit[[1L]],
  nsim = 4,
  seed = 942,
  parallel = "no"
)))
snow <- suppressMessages(suppressWarnings(mixedcirc_boottest(
  ordinary_fit[[1L]],
  nsim = 4,
  seed = 942,
  parallel = "snow",
  ncpus = 2
)))

stopifnot(
  methods::is(snow, "mixedcirc_boot"),
  isTRUE(all.equal(serial@results, snow@results, tolerance = 1e-10))
)

if (.Platform$OS.type != "windows") {
  multicore <- suppressMessages(suppressWarnings(mixedcirc_boottest(
    ordinary_fit[[1L]],
    nsim = 4,
    seed = 942,
    parallel = "multicore",
    ncpus = 2
  )))
  stopifnot(isTRUE(all.equal(serial@results, multicore@results, tolerance = 1e-10)))
}

mixed_fit <- suppressWarnings(mixedcirc_detect(
  simulation$data_input,
  simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  lm_method = "lme",
  random_effect_type = "mesor"
))

cluster <- parallel::makePSOCKcluster(2)
mixed_snow <- suppressMessages(suppressWarnings(mixedcirc_boottest(
  mixed_fit[[1L]],
  nsim = 3,
  seed = 943,
  parallel = "snow",
  ncpus = 2,
  cl = cluster
)))

worker_check <- parallel::clusterCall(cluster, function() TRUE)
parallel::stopCluster(cluster)
stopifnot(
  methods::is(mixed_snow, "mixedcirc_boot"),
  any(is.finite(mixed_snow@results$boot.estimate)),
  all(vapply(worker_check, isTRUE, logical(1)))
)
