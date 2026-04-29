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

.mc_extract_stats_from_fit <- function(fit, object, abs_phase = TRUE, wrap_phase = TRUE) {
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

    if (abs_phase) {
      phase <- abs(phase)
    }

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

  if (abs_phase) {
    phases <- abs(phases)
  }

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

.mc_boot_ci <- function(t0, tmat, type = c("perc", "basic", "norm"), conf_level = 0.95) {
  type <- match.arg(type)
  alpha <- 1 - conf_level

  res <- matrix(NA_real_, nrow = length(t0), ncol = 2)
  colnames(res) <- c("boot.lower", "boot.upper")
  rownames(res) <- names(t0)

  for (j in seq_along(t0)) {
    x <- tmat[, j]
    x <- x[is.finite(x)]

    if (length(x) < 5) next

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

.mc_boot_pvalues <- function(t0, tmat) {
  out <- numeric(length(t0))
  names(out) <- names(t0)

  for (j in seq_along(t0)) {
    x <- tmat[, j]
    x <- x[is.finite(x)]

    if (length(x) == 0) {
      out[j] <- NA_real_
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

#' Bootstrap Testing for Circadian Data
#'
#' This function performs bootstrap testing on fitted circadian models
#' from `mixedcirc_detect()`, supporting:
#' - lm and lmer fits
#' - RRBS and non-RRBS analyses
#' - grouped and ungrouped models
#'
#' The function uses `boot_lm()` to bootstrap mesor, amplitude, and phase
#' summaries derived from the fitted model. When grouped analyses are used,
#' pairwise group differences are also returned.
#'
#' By default, phase differences between groups are wrapped onto the circular
#' interval [-pi, pi], which represents the shortest signed phase difference.
#' This is usually the most meaningful representation in circadian analyses.
#'
#' @param object A fitted model object of class `mixedcirc_fit`.
#' @param nsim Number of bootstrap simulations.
#' @param seed Seed for reproducibility.
#' @param parallel The type of parallel processing to use.
#' @param ncpus Number of CPU cores to use.
#' @param boot.type Type of bootstrap confidence interval.
#' @param conf_level Confidence level for the bootstrap intervals.
#' @param abs_phase Logical; if `TRUE`, absolute phases are returned.
#' @param wrap_phase Logical; if `TRUE` (default), pairwise phase differences
#'   are wrapped to the interval [-pi, pi] before being returned. This gives
#'   the shortest signed circular phase difference between groups. Ignored for
#'   ungrouped analyses.
#' @param ... Additional arguments passed to `boot_lm()`.
#'
#' @return A `mixedcirc_boot` object.
#' @import stats
#' @import methods
#' @import emmeans
#' @import parallel
#' @import lme4
#' @export
mixedcirc_boottest <- function(object,
                               nsim,
                               seed,
                               parallel = c("no", "multicore", "snow"),
                               ncpus = getOption("boot.ncpus", 1L),
                               boot.type = c("perc", "basic", "norm")[1],
                               conf_level = 0.95,
                               abs_phase = TRUE,
                               wrap_phase = TRUE,
                               ...) {

  if (!inherits(object, "mixedcirc_fit")) {
    stop("Object must be mixedcirc_fit!")
  }

  if (missing(nsim) || length(nsim) != 1 || !is.numeric(nsim) || nsim < 1) {
    stop("`nsim` must be a positive integer.")
  }

  if (missing(seed) || length(seed) != 1 || !is.numeric(seed)) {
    stop("`seed` must be a single numeric value.")
  }

  fit <- object@fit

  model_type <- if (inherits(fit, "merMod")) "lme4" else "lm"

  coeff_fun <- function(x) {
    .mc_extract_stats_from_fit(
      fit = x,
      object = object,
      abs_phase = abs_phase,
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

  ci_df <- .mc_boot_ci(
    t0 = t0,
    tmat = tmat,
    type = boot.type,
    conf_level = conf_level
  )

  se_vec <- apply(tmat, 2, stats::sd, na.rm = TRUE)
  p_vals <- .mc_boot_pvalues(t0, tmat)

  results <- cbind(
    ci_df,
    boot.estimate = t0,
    boot.SE = se_vec,
    pvalue = p_vals
  )

  mixedcirc:::mixedcirc_boot(
    results = data.frame(results, check.names = FALSE),
    metadata = list(
      nsim = nsim,
      seed = seed,
      parallel = parallel,
      ncpus = ncpus,
      conf_level = conf_level,
      boot.type = boot.type,
      abs_phase = abs_phase,
      wrap_phase = wrap_phase,
      model_class = model_type
    ),
    original_model = object
  )
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

boot_lm <-
  function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
            type = c("parametric", "semiparametric"), verbose = FALSE,
            .progress = "none", PBargs = list(), parallel = c("no", "multicore",
                                                              "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL)
  {
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    if (.progress != "none") {
      pbfun <- get(paste0(.progress, "ProgressBar"))
      setpbfun <- get(paste0("set", .simpleCap(.progress),
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
        refit
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- lme4:::factory(function(i) FUN(refit(x, ss[[i]])),
                             errval = rep(NA, length.t0))
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
          if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
            parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
              parallel::clusterSetRNGStream(cl)
            res <- parallel::parLapply(cl, simvec, ffun)
            parallel::stopCluster(cl)
            res
          } else parallel::parLapply(cl, simvec, ffun)
        }
      } else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning", "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      if ((numFail <- sum(msgs[["factory-error"]])) > 0) {
        warning("some bootstrap runs failed (", numFail, "/", nsim, ")")
      }
      fail.msgs <- if (numFail == 0) NULL else msgs[["factory-error"]]
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
        refit
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- lme4:::factory(function(i) FUN(refit(x, ss[[i]], control = control)),
                             errval = rep(NA, length.t0))
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
          if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
            parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
              parallel::clusterSetRNGStream(cl)
            res <- parallel::parLapply(cl, simvec, ffun)
            parallel::stopCluster(cl)
            res
          } else parallel::parLapply(cl, simvec, ffun)
        }
      } else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning", "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      if ((numFail <- sum(msgs[["factory-error"]])) > 0) {
        warning("some bootstrap runs failed (", numFail, "/", nsim, ")")
      }
      fail.msgs <- if (numFail == 0) NULL else msgs[["factory-error"]]
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
