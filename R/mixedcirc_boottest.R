# =========================================================
# Helper functions for mixedcirc bootstrap
# =========================================================

.mc_get_call_arg_name <- function(object, arg_name) {
  if (!"call" %in% methods::slotNames(object)) {
    return(NULL)
  }

  cl <- object@call
  if (is.null(cl) || is.null(cl[[arg_name]])) {
    return(NULL)
  }

  out <- tryCatch(as.character(cl[[arg_name]]), error = function(e) NULL)
  if (length(out) == 0) {
    return(NULL)
  }
  out[1]
}

.mc_get_boot_spec <- function(object) {
  params <- if ("params" %in% methods::slotNames(object)) object@params else list()

  list(
    group_var = .mc_get_call_arg_name(object, "group"),
    lmer_df = .mc_get_call_arg_name(object, "lmer.df"),
    RRBS = isTRUE(params$RRBS),
    multiple_groups = isTRUE(params$multiple_groups)
  )
}

.mc_extract_estimate_col <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)

  trend_cols <- grep("\\.trend$", names(x), value = TRUE)
  if (length(trend_cols) > 0) {
    return(as.numeric(x[[trend_cols[1]]]))
  }

  if ("emmean" %in% names(x)) {
    return(as.numeric(x[["emmean"]]))
  }

  num_cols <- which(vapply(x, is.numeric, logical(1)))
  if (length(num_cols) == 0) {
    stop("Could not identify numeric estimate column.")
  }

  as.numeric(x[[num_cols[1]]])
}

.mc_pair_names <- function(x) {
  if (length(x) < 2) return(list())
  combn(x, 2, simplify = FALSE)
}

.mc_safe_calc_phase <- function(beta_in, beta_out) {
  out <- tryCatch(calc_phase(beta_in, beta_out), error = function(e) NA_real_)
  as.numeric(out)[1]
}

.mc_wrap_phase_diff <- function(phi1, phi2) {
  d <- phi1 - phi2
  ((d + pi) %% (2 * pi)) - pi
}

.mc_extract_stats_from_fit <- function(fit, object, wrap_phase = TRUE) {
  spec <- .mc_get_boot_spec(object)
  rrbs <- spec$RRBS
  grouped <- spec$multiple_groups
  group_var <- spec$group_var
  lmer_df <- spec$lmer_df

  if (is.null(lmer_df) || is.na(lmer_df) || !nzchar(lmer_df)) {
    lmer_df <- "Satterthwaite"
  }

  out <- list()

  if (!grouped) {
    if (!rrbs) {
      em_mesor <- as.data.frame(
        emmeans::emmeans(
          fit,
          ~ 1,
          at = list(
            inphase = 0,
            outphase = 0
          ),
          lmer.df = lmer_df
        )
      )

      em_in <- as.data.frame(
        emmeans::emtrends(
          fit,
          ~ 1,
          var = "inphase",
          lmer.df = lmer_df
        )
      )

      em_out <- as.data.frame(
        emmeans::emtrends(
          fit,
          ~ 1,
          var = "outphase",
          lmer.df = lmer_df
        )
      )
    } else {
      em_mesor <- as.data.frame(
        emmeans::emtrends(
          fit,
          ~ 1,
          var = "scaler",
          at = list(
            inphase = 0,
            outphase = 0,
            scaler = 1
          ),
          lmer.df = lmer_df
        )
      )

      em_in <- as.data.frame(
        emmeans::emtrends(
          fit,
          ~ scaler,
          var = "inphase",
          at = list(scaler = 1),
          lmer.df = lmer_df
        )
      )

      em_out <- as.data.frame(
        emmeans::emtrends(
          fit,
          ~ scaler,
          var = "outphase",
          at = list(scaler = 1),
          lmer.df = lmer_df
        )
      )
    }

    mesor <- .mc_extract_estimate_col(em_mesor)[1]
    beta_in <- .mc_extract_estimate_col(em_in)[1]
    beta_out <- .mc_extract_estimate_col(em_out)[1]

    amp <- sqrt(beta_in^2 + beta_out^2)
    phase <- .mc_safe_calc_phase(beta_in, beta_out)

    out[["mesor"]] <- mesor
    out[["amps"]] <- amp
    out[["phase"]] <- phase

    return(unlist(out, use.names = TRUE))
  }

  # grouped case
  formula_for_group <- stats::as.formula(paste0("~ 1 + ", group_var))

  if (!rrbs) {
    em_mesor <- as.data.frame(
      emmeans::emmeans(
        fit,
        formula_for_group,
        at = list(
          inphase = 0,
          outphase = 0
        ),
        lmer.df = lmer_df
      )
    )

    em_in <- as.data.frame(
      emmeans::emtrends(
        fit,
        formula_for_group,
        var = "inphase",
        lmer.df = lmer_df
      )
    )

    em_out <- as.data.frame(
      emmeans::emtrends(
        fit,
        formula_for_group,
        var = "outphase",
        lmer.df = lmer_df
      )
    )
  } else {
    em_mesor <- as.data.frame(
      emmeans::emtrends(
        fit,
        formula_for_group,
        var = "scaler",
        at = list(
          inphase = 0,
          outphase = 0,
          scaler = 1
        ),
        lmer.df = lmer_df
      )
    )

    em_in <- as.data.frame(
      emmeans::emtrends(
        fit,
        formula_for_group,
        var = "inphase",
        at = list(scaler = 1),
        lmer.df = lmer_df
      )
    )

    em_out <- as.data.frame(
      emmeans::emtrends(
        fit,
        formula_for_group,
        var = "outphase",
        at = list(scaler = 1),
        lmer.df = lmer_df
      )
    )
  }

  groups <- as.character(em_mesor[[group_var]])
  mesors <- .mc_extract_estimate_col(em_mesor)
  beta_in <- .mc_extract_estimate_col(em_in)
  beta_out <- .mc_extract_estimate_col(em_out)

  amps <- sqrt(beta_in^2 + beta_out^2)
  phases <- mapply(.mc_safe_calc_phase, beta_in, beta_out)

  for (i in seq_along(groups)) {
    g <- groups[i]
    out[[paste0("mesor_", g)]] <- mesors[i]
    out[[paste0("amps_", g)]] <- amps[i]
    out[[paste0("phase_", g)]] <- phases[i]
  }

  # pairwise differences
  for (pair in .mc_pair_names(groups)) {
    g1 <- pair[1]
    g2 <- pair[2]
    i1 <- match(g1, groups)
    i2 <- match(g2, groups)

    out[[paste0("mesor_", g1, "vs", g2)]] <- mesors[i1] - mesors[i2]
    out[[paste0("amps_", g1, "vs", g2)]] <- amps[i1] - amps[i2]

    if (wrap_phase) {
      out[[paste0("phase_", g1, "vs", g2)]] <- .mc_wrap_phase_diff(phases[i1], phases[i2])
    } else {
      out[[paste0("phase_", g1, "vs", g2)]] <- phases[i1] - phases[i2]
    }
  }

  unlist(out, use.names = TRUE)
}

## Wrap a value (or vector) onto the circular interval [-pi, pi).
.mc_circular_wrap <- function(x) ((x + pi) %% (2 * pi)) - pi

## Names beginning with "phase" represent circular quantities. Their intervals
## are calculated from wrapped deviations around the point estimate so that an
## interval remains local when its endpoints cross the -pi/pi boundary.
.mc_is_circular_param <- function(param_names) grepl("^phase", param_names)

.mc_boot_ci <- function(t0, tmat, type = c("perc", "basic", "norm"), conf_level = 0.95,
                         circular = NULL) {
  type <- match.arg(type)
  alpha <- 1 - conf_level

  res <- matrix(NA_real_, nrow = length(t0), ncol = 2)
  colnames(res) <- c("boot.lower", "boot.upper")
  rownames(res) <- names(t0)

  if (is.null(circular)) circular <- .mc_is_circular_param(names(t0))

  for (j in seq_along(t0)) {
    x <- tmat[, j]
    x <- x[is.finite(x)]

    if (length(x) < 5) next

    if (isTRUE(circular[j])) {
      ## Work on wrapped deviations, then map the endpoints back to the circle.
      dev <- .mc_circular_wrap(x - t0[j])

      if (type == "perc") {
        qs <- stats::quantile(dev, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
        res[j, ] <- .mc_circular_wrap(t0[j] + qs)
      } else if (type == "basic") {
        qs <- stats::quantile(dev, probs = c(1 - alpha / 2, alpha / 2), na.rm = TRUE, names = FALSE)
        res[j, ] <- .mc_circular_wrap(t0[j] - qs)
      } else if (type == "norm") {
        se <- stats::sd(dev, na.rm = TRUE)
        z <- stats::qnorm(1 - alpha / 2)
        res[j, ] <- .mc_circular_wrap(t0[j] + c(-z * se, z * se))
      }
      next
    }

    if (type == "perc") {
      qs <- stats::quantile(x, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
      res[j, ] <- qs
    } else if (type == "basic") {
      qs <- stats::quantile(x, probs = c(1 - alpha / 2, alpha / 2), na.rm = TRUE, names = FALSE)
      res[j, ] <- c(2 * t0[j] - qs[1], 2 * t0[j] - qs[2])
    } else if (type == "norm") {
      se <- stats::sd(x, na.rm = TRUE)
      z <- stats::qnorm(1 - alpha / 2)
      res[j, ] <- c(t0[j] - z * se, t0[j] + z * se)
    }
  }

  as.data.frame(res, check.names = FALSE)
}

.mc_boot_pvalues <- function(t0, tmat, circular = NULL) {
  out <- numeric(length(t0))
  names(out) <- names(t0)

  if (is.null(circular)) circular <- .mc_is_circular_param(names(t0))

  for (j in seq_along(t0)) {
    x <- tmat[, j]
    x <- x[is.finite(x)]

    if (length(x) == 0) {
      out[j] <- NA_real_
      next
    }

    if (isTRUE(circular[j])) {
      ## Compare the circular distance to zero with deviations around the
      ## estimated phase; an arithmetic mean is unsuitable at the branch cut.
      dev <- .mc_circular_wrap(x - t0[j])
      null_dist <- abs(.mc_circular_wrap(-t0[j]))
      out[j] <- (sum(abs(dev) >= null_dist, na.rm = TRUE) + 1) / (length(x) + 1)
      next
    }

    z_underH0 <- x - mean(x, na.rm = TRUE)
    out[j] <- (sum(abs(z_underH0) > abs(t0[j]), na.rm = TRUE) + 1) / (length(x) + 1)
  }

  out
}

# =========================================================
# Main bootstrap function
# =========================================================

#' Parametric Bootstrap Uncertainty for Circadian Rhythm Summaries
#'
#' Uses a parametric bootstrap to estimate uncertainty for mesor, amplitude,
#' and phase summaries from a fitted `mixedcirc_detect()` model, and
#' for the supported between-group contrasts of these quantities when a
#' grouped (multi-condition) model is used. Supports `lm` and `lmer` fits,
#' RRBS and non-RRBS analyses, and grouped and ungrouped models.
#'
#' The bootstrap is **parametric**, not case resampling: for each replicate,
#' new response values are simulated from the fitted model (`stats::simulate()`
#' on the retained `lm`/`merMod` object, using its estimated fixed effects and,
#' for mixed models, its estimated random-effect and residual variance
#' components), the model is refit to the simulated data, and the requested
#' summaries are recomputed from the refit. Subjects and observation times are
#' not resampled; only the response is resimulated. This differs from a
#' case-resampling (nonparametric) bootstrap, in which subjects/rows are
#' resampled with replacement from the observed data.
#'
#' Phase is circular. Confidence intervals, standard errors, and bootstrap
#' p-values for parameters whose names begin with `"phase"` are therefore
#' calculated from wrapped deviations around the point estimate. This prevents
#' intervals near the `-pi`/`pi` boundary from spanning almost the full cycle.
#' Mesor and amplitude are handled on their ordinary linear scales.
#'
#' Pairwise group differences follow the convention `<group_1>vs<group_2> =
#' <group_1> - <group_2>` for mesor and amplitude; the corresponding phase
#' difference is the shortest signed circular difference
#' `wrap(phase_group_1 - phase_group_2)` when `wrap_phase = TRUE`. Individual
#' phase estimates retain the signed model convention used by
#' `mixedcirc_detect()`.
#'
#' @param object A fitted model object of class `mixedcirc_fit`.
#' @param nsim Number of parametric bootstrap replicates.
#' @param seed Seed for reproducibility.
#' @param parallel The type of parallel processing to use.
#' @param ncpus Number of CPU cores to use.
#' @param boot.type Type of bootstrap confidence interval (`"perc"`, `"basic"`,
#'   or `"norm"`); applied on the circular deviation scale for phase
#'   parameters and on the raw scale for mesor/amplitude.
#' @param conf_level Confidence level for the bootstrap intervals.
#' @param wrap_phase Logical; if `TRUE` (default), pairwise phase differences
#'   are wrapped to the interval `[-pi, pi)` before being returned. This gives
#'   the shortest signed circular phase difference between groups. Ignored for
#'   ungrouped analyses.
#' @param ... Additional arguments passed to `boot_lm()`.
#'
#' @return A `mixedcirc_boot` object. Its `results` data frame reports
#'   `boot.lower`/`boot.upper` (bootstrap confidence interval), `boot.estimate`
#'   (point estimate), `boot.SE` (bootstrap standard deviation - the circular
#'   deviation-based SD for phase parameters), and `pvalue` (bootstrap p-value
#'   for the null that the parameter is zero) for every returned mesor,
#'   amplitude, and phase summary, including the between-group
#'   contrasts when a grouped model is used. Phase-scale parameters are in
#'   radians; convert a phase difference to time units using the fitted
#'   model's period if a peak-time difference is preferred, e.g.
#'   `delta_t <- ((-delta_phi * period / (2 * pi) + period / 2) %% period) -
#'   period / 2` for a result wrapped to `[-period/2, period/2)`.
#' @export
mixedcirc_boottest <- function(object,
                               nsim,
                               seed,
                               parallel = c("no", "multicore", "snow"),
                               ncpus = getOption("boot.ncpus", 1L),
                               boot.type = c("perc", "basic", "norm")[1],
                               conf_level = 0.95,
                               wrap_phase = TRUE,
                               ...) {

  if (!inherits(object, "mixedcirc_fit")) {
    stop("object must be a mixedcirc_fit", call. = FALSE)
  }

  if (missing(nsim) || length(nsim) != 1L || !is.numeric(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim %% 1 != 0) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (missing(seed) || length(seed) != 1L || !is.numeric(seed) ||
      !is.finite(seed) || seed %% 1 != 0) {
    stop("`seed` must be a single finite integer.")
  }
  seed <- as.integer(seed)
  if (!is.character(parallel) || length(parallel) < 1L || anyNA(parallel)) {
    stop("`parallel` must be one of 'no', 'multicore', or 'snow'.")
  }
  parallel <- match.arg(parallel, c("no", "multicore", "snow"))
  if (!is.numeric(ncpus) || length(ncpus) != 1L || !is.finite(ncpus) ||
      ncpus < 1 || ncpus %% 1 != 0) {
    stop("`ncpus` must be a positive integer.")
  }
  ncpus <- as.integer(ncpus)
  boot.type <- match.arg(boot.type, c("perc", "basic", "norm"))
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be between 0 and 1.")
  }
  if (!is.logical(wrap_phase) || length(wrap_phase) != 1L || is.na(wrap_phase)) {
    stop("`wrap_phase` must be TRUE or FALSE.")
  }

  fit <- object@fit
  if (is.null(fit)) {
    stop(
      "The fitted model was not retained. Refit with store_models = TRUE ",
      "before running the bootstrap.",
      call. = FALSE
    )
  }

  model_type <- if (inherits(fit, "merMod")) "lme4" else "lm"

  coeff_fun <- function(x) {
    .mc_extract_stats_from_fit(
      fit = x,
      object = object,
      wrap_phase = wrap_phase
    )
  }

  boot_results <- boot_lm(
    x = fit,
    FUN = coeff_fun,
    nsim = nsim,
    seed = seed,
    parallel = parallel,
    ncpus = ncpus,
    ...
  )

  t0 <- boot_results$t0
  tmat <- boot_results$t
  colnames(tmat) <- names(t0)

  circular_param <- .mc_is_circular_param(names(t0))

  ci_df <- .mc_boot_ci(
    t0 = t0,
    tmat = tmat,
    type = boot.type,
    conf_level = conf_level,
    circular = circular_param
  )

  se_vec <- vapply(seq_along(t0), function(j) {
    x <- tmat[, j]
    if (isTRUE(circular_param[j])) {
      stats::sd(.mc_circular_wrap(x - t0[j]), na.rm = TRUE)
    } else {
      stats::sd(x, na.rm = TRUE)
    }
  }, numeric(1))
  names(se_vec) <- names(t0)
  p_vals <- .mc_boot_pvalues(t0, tmat, circular = circular_param)

  results <- cbind(
    ci_df,
    boot.estimate = t0,
    boot.SE = se_vec,
    pvalue = p_vals
  )

  mixedcirc_boot(
    results = data.frame(results, check.names = FALSE),
    metadata = list(
      nsim = nsim,
      seed = seed,
      parallel = parallel,
      ncpus = ncpus,
      conf_level = conf_level,
      boot.type = boot.type,
      wrap_phase = wrap_phase,
      model_class = model_type
    ),
    original_model = object
  )
}

#' @export
refit.lm <- function(object, newresp, ...) {
  if (!inherits(object, "lm")) {
    stop("`object` must inherit from class \"lm\".", call. = FALSE)
  }
  if (inherits(object, "glm") || inherits(object, "mlm")) {
    stop(
      "`refit.lm()` supports univariate Gaussian `lm` objects only.",
      call. = FALSE
    )
  }

  if (is.data.frame(newresp) || is.matrix(newresp)) {
    if (NCOL(newresp) != 1L) {
      stop("`newresp` must contain a single response.", call. = FALSE)
    }
    newresp <- newresp[, 1L]
  }

  model_matrix <- stats::model.matrix(object)
  n_used <- nrow(model_matrix)

  # simulate.lm() restores excluded observations for an na.exclude fit. Remove
  # those positions again so that the response matches the fitted model matrix.
  if (length(newresp) != n_used && !is.null(object$na.action)) {
    omitted <- as.integer(object$na.action)
    if (length(newresp) - length(omitted) == n_used) {
      newresp <- newresp[-omitted]
    }
  }

  if (length(newresp) != n_used) {
    stop(
      sprintf(
        "`newresp` has length %d, but the fitted model uses %d observations.",
        length(newresp),
        n_used
      ),
      call. = FALSE
    )
  }

  tol <- if (!is.null(object$qr$tol)) object$qr$tol else 1e-7
  refitted_components <- if (is.null(object$weights)) {
    stats::lm.fit(
      x = model_matrix,
      y = newresp,
      offset = object$offset,
      tol = tol
    )
  } else {
    stats::lm.wfit(
      x = model_matrix,
      y = newresp,
      w = object$weights,
      offset = object$offset,
      tol = tol
    )
  }

  # Keep the original terms, call, contrasts, factor levels, and model frame.
  # Downstream methods such as emmeans need this metadata from a full lm object.
  refitted <- object
  for (component in names(refitted_components)) {
    refitted[[component]] <- refitted_components[[component]]
  }

  if (!is.null(refitted$model)) {
    refitted$model[[1L]] <- newresp
  }
  if (!is.null(refitted$y)) {
    refitted$y <- newresp
  }

  refitted
}

initialize.parallel <- expression({
  have_mc <- have_snow <- FALSE
  if (length(parallel)>1) parallel <- match.arg(parallel)
  do_parallel <- (parallel != "no" && ncpus > 1L)
  if (do_parallel) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!(have_mc || have_snow))
      do_parallel <- FALSE
  }
})

# Capitalize the first character without relying on the unexported
# utils:::.simpleCap() helper.
.mc_simple_cap <- function(x) {
  paste0(toupper(substr(x, 1L, 1L)), substring(x, 2L))
}

# Evaluate one bootstrap replicate while retaining messages, warnings, and
# errors as attributes. This is the small piece of behaviour mixedcirc needs
# from lme4's private factory() utility, implemented locally so that lme4
# internals are not part of our API contract.
.mc_boot_factory <- function(fun, errval) {
  force(fun)
  force(errval)

  function(...) {
    messages <- warnings <- error <- NULL
    error_occurred <- FALSE

    result <- withCallingHandlers(
      tryCatch(
        fun(...),
        error = function(e) {
          error <<- conditionMessage(e)
          error_occurred <<- TRUE
          NULL
        }
      ),
      warning = function(w) {
        warnings <<- append(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- append(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )

    if (error_occurred) {
      result <- errval
    }
    attr(result, "factory-message") <- messages
    attr(result, "factory-warning") <- warnings
    attr(result, "factory-error") <- error
    result
  }
}

.mc_boot_snow_apply <- function(cl, indices, fun, ncpus) {
  owns_cluster <- is.null(cl)
  if (owns_cluster) {
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  parallel::clusterCall(cl, function(packages) {
    invisible(lapply(packages, loadNamespace))
    TRUE
  }, c("mixedcirc", "lme4", "emmeans"))

  if (owns_cluster && RNGkind()[1L] == "L'Ecuyer-CMRG") {
    parallel::clusterSetRNGStream(cl)
  }
  parallel::parLapply(cl, indices, fun)
}

boot_lm <-
  function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
            type = c("parametric", "semiparametric"), verbose = FALSE,
            .progress = "none", PBargs = list(), parallel = c("no", "multicore",
                                                              "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL)
  {
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    if (.progress != "none") {
      pbfun <- get(paste0(.progress, "ProgressBar"))
      setpbfun <- get(paste0("set", .mc_simple_cap(.progress),
                             "ProgressBar"))
      pb <- do.call(pbfun, PBargs)
    }
    do_parallel <- have_mc <- have_snow <- NULL
    eval(initialize.parallel)
    if (do_parallel && .progress != "none")
      message("progress bar disabled for parallel operations")
    FUN <- match.fun(FUN)
    type <- match.arg(type)
    if (!is.null(seed))
      set.seed(seed)
    else if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
      stop("This function currently only handles functions that return numeric vectors")

    if (is(x, "lm")) {
      if (type == "parametric") {
        argList <- list(x, nsim = nsim, na.action = na.exclude)
        if (!missing(re.form)) {
          argList <- c(argList)
        } else {
          argList <- c(argList)
        }
        ss <- do.call(simulate, argList)
      } else {
        stop("semiparametric bootstrapping with lm not yet implemented")
      }
      control <- if (!is(x, "lm"))
        NULL
      else eval.parent(x$call$control)
      ffun <- local({
        FUN
        refit_fun <- refit.lm
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- .mc_boot_factory(
          function(i) FUN(refit_fun(x, ss[[i]])),
          errval = rep(NA, length.t0)
        )
        function(i) {
          ret <- f1(i)
          if (verbose) {
            cat(sprintf("%5d :", i))
            str(ret)
          }
          if (!do_parallel && .progress != "none") {
            setpbfun(pb, i/nsim)
          }
          ret
        }
      })
      simvec <- seq_len(nsim)
      res <- if (do_parallel) {
        if (have_mc) {
          parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        } else if (have_snow) {
          .mc_boot_snow_apply(cl, simvec, ffun, ncpus)
        }
      } else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning", "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      numFail <- sum(msgs[["factory-error"]])
      fail.msgs <- if (numFail == 0) NULL else msgs[["factory-error"]]
      if (numFail >= nsim) {
        error_details <- names(fail.msgs)
        detail <- if (length(error_details) > 0) {
          paste0(" Error: ", paste(error_details, collapse = "; "))
        } else {
          ""
        }
        stop(
          "All bootstrap runs failed (", numFail, "/", nsim, ").",
          detail,
          call. = FALSE
        )
      }
      if (numFail > 0) {
        warning("some bootstrap runs failed (", numFail, "/", nsim, ")")
      }
      s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = model.frame(x),
                          seed = .Random.seed, statistic = FUN, sim = "parametric",
                          call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = FALSE),
                     class = c("bootMer", "boot"))
      attr(s, "bootFail") <- numFail
      attr(s, "boot.fail.msgs") <- fail.msgs
      attr(s, "boot.all.msgs") <- msgs
      attr(s, "boot_type") <- "boot"
      s
    } else {
      mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
      if (isLMM(x))
        mle <- c(mle, list(sigma = sigma(x)))
      if (type == "parametric") {
        argList <- list(x, nsim = nsim, na.action = na.exclude)
        if (!missing(re.form)) {
          argList <- c(argList, list(re.form = re.form))
        } else {
          argList <- c(argList, list(use.u = use.u))
        }
        ss <- do.call(simulate, argList)
      } else {
        if (!missing(re.form))
          stop(paste(sQuote("re.form")), "cannot be used with semiparametric bootstrapping")
        if (use.u) {
          if (isGLMM(x))
            warning("semiparametric bootstrapping is questionable for GLMMs")
          ss <- replicate(nsim, fitted(x) + sample(residuals(x, "response"), replace = TRUE), simplify = FALSE)
        } else {
          stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
        }
      }
      control <- if (!is(x, "merMod"))
        NULL
      else eval.parent(x@call$control)
      ffun <- local({
        FUN
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- .mc_boot_factory(
          function(i) FUN(lme4::refit(x, ss[[i]], control = control)),
          errval = rep(NA, length.t0)
        )
        function(i) {
          ret <- f1(i)
          if (verbose) {
            cat(sprintf("%5d :", i))
            str(ret)
          }
          if (!do_parallel && .progress != "none") {
            setpbfun(pb, i/nsim)
          }
          ret
        }
      })
      simvec <- seq_len(nsim)
      res <- if (do_parallel) {
        if (have_mc) {
          parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        } else if (have_snow) {
          .mc_boot_snow_apply(cl, simvec, ffun, ncpus)
        }
      } else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning", "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      numFail <- sum(msgs[["factory-error"]])
      fail.msgs <- if (numFail == 0) NULL else msgs[["factory-error"]]
      if (numFail >= nsim) {
        error_details <- names(fail.msgs)
        detail <- if (length(error_details) > 0) {
          paste0(" Error: ", paste(error_details, collapse = "; "))
        } else {
          ""
        }
        stop(
          "All bootstrap runs failed (", numFail, "/", nsim, ").",
          detail,
          call. = FALSE
        )
      }
      if (numFail > 0) {
        warning("some bootstrap runs failed (", numFail, "/", nsim, ")")
      }
      s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = model.frame(x),
                          seed = .Random.seed, statistic = FUN, sim = "parametric",
                          call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                     class = c("bootMer", "boot"))
      attr(s, "bootFail") <- numFail
      attr(s, "boot.fail.msgs") <- fail.msgs
      attr(s, "boot.all.msgs") <- msgs
      attr(s, "boot_type") <- "boot"
      s
    }
  }
