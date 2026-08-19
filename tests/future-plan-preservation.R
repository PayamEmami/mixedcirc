library(mixedcirc)

test_future_plan_preservation <- function() {
  original_plan <- future::plan("list")
  on.exit(future::plan(original_plan), add = TRUE)
  original_foreach_backend <- c(
    name = foreach::getDoParName(),
    workers = as.character(foreach::getDoParWorkers())
  )

  plan_signature <- function(plan_list) {
    lapply(plan_list, function(strategy) {
      list(
        class = class(strategy),
        call = paste(deparse(attr(strategy, "call")), collapse = ""),
        tweaks = attr(strategy, "tweaks")
      )
    })
  }

  assert_plan_unchanged <- function(before, expected_workers) {
    after <- future::plan("list")
    stopifnot(
      identical(plan_signature(after), plan_signature(before)),
      future::nbrOfWorkers() == expected_workers,
      identical(
        c(
          name = foreach::getDoParName(),
          workers = as.character(foreach::getDoParWorkers())
        ),
        original_foreach_backend
      )
    )

    probe <- future::future(42L)
    stopifnot(identical(future::value(probe), 42L))
  }

  expect_input_error <- function(code) {
    error_message <- tryCatch(code, error = conditionMessage)
    stopifnot(identical(
      error_message,
      "data_input must be a data frame, matrix, or DGEList"
    ))
  }

  expect_runtime_error <- function(code) {
    error_message <- tryCatch(code, error = conditionMessage)
    stopifnot(is.character(error_message), length(error_message) == 1L)
  }

  simulated <- mixedcirc_simulate_data(
    n_subjects_per_group = 3,
    n_timepoints = 6,
    n_features = 2,
    groups = c("A", "B"),
    seed = 601
  )

  # A sequential caller must remain sequential after both serial and parallel
  # package operations.
  future::plan(future::sequential)

  before <- future::plan("list")
  detect_serial <- mixedcirc_detect(
    simulated$data_input,
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 1
  )
  assert_plan_unchanged(before, expected_workers = 1L)

  before <- future::plan("list")
  detect_parallel <- mixedcirc_detect(
    simulated$data_input,
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 2
  )
  assert_plan_unchanged(before, expected_workers = 1L)

  serial_results <- lapply(detect_serial@results, function(x) x@results)
  parallel_results <- lapply(detect_parallel@results, function(x) x@results)
  stopifnot(isTRUE(all.equal(serial_results, parallel_results, tolerance = 1e-10)))

  before <- future::plan("list")
  detrend_serial <- mixedcirc_remove_trend(
    simulated$data_input,
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 1
  )
  assert_plan_unchanged(before, expected_workers = 1L)

  before <- future::plan("list")
  detrend_parallel <- mixedcirc_remove_trend(
    simulated$data_input,
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 2
  )
  assert_plan_unchanged(before, expected_workers = 1L)

  stopifnot(isTRUE(all.equal(
    detrend_serial,
    detrend_parallel,
    tolerance = 1e-10
  )))

  before <- future::plan("list")
  expect_input_error(mixedcirc_detect(
    data_input = NULL,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    ncores = 2
  ))
  assert_plan_unchanged(before, expected_workers = 1L)

  before <- future::plan("list")
  expect_input_error(mixedcirc_remove_trend(
    data_input = NULL,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    ncores = 2
  ))
  assert_plan_unchanged(before, expected_workers = 1L)

  # Errors raised inside a parallel feature fit must also restore the plan.
  invalid_data <- simulated$data_input[, 1L, drop = FALSE]
  invalid_data[, 1L] <- NA_real_

  before <- future::plan("list")
  expect_runtime_error(mixedcirc_detect(
    data_input = invalid_data,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 2
  ))
  assert_plan_unchanged(before, expected_workers = 1L)

  before <- future::plan("list")
  expect_runtime_error(mixedcirc_remove_trend(
    data_input = invalid_data,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lm",
    random_effect_type = "none",
    ncores = 2
  ))
  assert_plan_unchanged(before, expected_workers = 1L)

  # A caller-configured multisession plan and its worker count must survive
  # serial and parallel package operations.
  future::plan(future::multisession, workers = 2)

  before <- future::plan("list")
  mixed_fit <- mixedcirc_detect(
    simulated$data_input[, 1L, drop = FALSE],
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lme",
    random_effect_type = "mesor",
    ncores = 1
  )
  stopifnot(inherits(mixed_fit[[1L]]@fit, "merMod"))
  assert_plan_unchanged(before, expected_workers = 2L)

  before <- future::plan("list")
  mixed_detrended <- mixedcirc_remove_trend(
    simulated$data_input,
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lme",
    random_effect_type = "mesor",
    ncores = 2
  )
  stopifnot(
    identical(dim(mixed_detrended), dim(simulated$data_input)),
    all(is.finite(mixed_detrended))
  )
  assert_plan_unchanged(before, expected_workers = 2L)

  before <- future::plan("list")
  expect_input_error(mixedcirc_detect(
    data_input = NULL,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    ncores = 1
  ))
  assert_plan_unchanged(before, expected_workers = 2L)

  # Preserve multi-level plans as a complete strategy list, not merely the
  # class of the first strategy.
  future::plan(list(
    future::tweak(future::multisession, workers = 2),
    future::sequential
  ))
  before <- future::plan("list")

  expect_input_error(mixedcirc_remove_trend(
    data_input = NULL,
    meta_input = simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    ncores = 1
  ))
  assert_plan_unchanged(before, expected_workers = 2L)
}

test_future_plan_preservation()
