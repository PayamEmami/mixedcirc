library(mixedcirc)

test_dgelist_dream_weights_once <- function() {
  counter_name <- ".mixedcirc_test_voom_call_counter"
  assign(counter_name, 0L, envir = .GlobalEnv)
  on.exit(rm(list = counter_name, envir = .GlobalEnv), add = TRUE)

  suppressMessages(
    trace(
      "voomWithDreamWeights",
      tracer = quote(
        assign(
          ".mixedcirc_test_voom_call_counter",
          get(".mixedcirc_test_voom_call_counter", envir = .GlobalEnv) + 1L,
          envir = .GlobalEnv
        )
      ),
      where = asNamespace("mixedcirc"),
      print = FALSE
    )
  )
  on.exit(
    suppressMessages(
      untrace("voomWithDreamWeights", where = asNamespace("mixedcirc"))
    ),
    add = TRUE
  )

  simulated <- mixedcirc_simulate_data(
    n_subjects_per_group = 3,
    n_timepoints = 6,
    n_features = 2,
    groups = c("A", "B"),
    seed = 401
  )

  set.seed(402)
  counts <- matrix(
    stats::rnbinom(
      2L * nrow(simulated$meta_input),
      mu = 100,
      size = 10
    ),
    nrow = 2L,
    dimnames = list(
      paste0("count_feature_", 1:2),
      rownames(simulated$meta_input)
    )
  )

  fit <- mixedcirc_detect(
    edgeR::DGEList(counts = counts),
    simulated$meta_input,
    time = "time",
    group = "group",
    id = "id",
    lm_method = "lme",
    random_effect_type = "mesor",
    ncores = 1
  )

  stopifnot(
    get(counter_name, envir = .GlobalEnv) == 1L,
    methods::is(fit, "mixedcirc_fit_list"),
    length(fit) == 2L,
    all(vapply(fit@results, function(x) identical(x@type, "expression"), logical(1))),
    all(vapply(fit@results, function(x) inherits(x@fit, "merMod"), logical(1)))
  )
}

test_dgelist_dream_weights_once()
