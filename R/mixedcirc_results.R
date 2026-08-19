# Publication-oriented result extraction ------------------------------------

.mc_fit_items <- function(x) {
  if (methods::is(x, "mixedcirc_fit")) {
    return(list(x))
  }
  if (methods::is(x, "mixedcirc_fit_list")) {
    return(x@results)
  }
  stop("x must be a mixedcirc_fit or mixedcirc_fit_list object", call. = FALSE)
}

.mc_feature_name <- function(object, fallback = "feature") {
  nm <- tryCatch(rownames(object@results), error = function(e) NULL)
  if (is.null(nm) || length(nm) == 0L || is.na(nm[1L]) || !nzchar(nm[1L])) {
    return(fallback)
  }
  as.character(nm[1L])
}

.mc_result_number <- function(object, name) {
  if (!name %in% names(object@results) || nrow(object@results) == 0L) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(object@results[[name]][1L]))
  if (length(value) == 0L) NA_real_ else value[1L]
}

.mc_group_levels <- function(object) {
  params <- object@params
  if (!isTRUE(params$multiple_groups)) {
    return("all")
  }

  group_var <- params$group
  if (is.null(group_var) || !group_var %in% names(object@exp_design)) {
    return(character())
  }

  values <- object@exp_design[[group_var]]
  if (is.factor(values)) {
    observed <- levels(droplevels(values[!is.na(values)]))
  } else {
    observed <- unique(as.character(values[!is.na(values)]))
  }
  as.character(observed)
}

.mc_wrap_radians <- function(x) {
  ((x + pi) %% (2 * pi)) - pi
}

.mc_wrap_period <- function(x, period) {
  ((x + period / 2) %% period) - period / 2
}

.mc_bind_rows_union <- function(x) {
  if (length(x) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing_names <- setdiff(all_names, names(z))
    for (nm in missing_names) z[[nm]] <- NA
    z[all_names]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.mc_adjust_column <- function(x, method) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (any(keep)) out[keep] <- stats::p.adjust(x[keep], method = method)
  out
}

.mc_model_diagnostics <- function(object, fallback = "feature") {
  fit <- object@fit
  params <- object@params
  status <- if (!is.null(params$fit_status)) {
    as.character(params$fit_status)[1L]
  } else if (is.null(fit)) {
    "failed"
  } else {
    "ok"
  }

  message <- if (!is.null(params$fit_message)) {
    paste(as.character(params$fit_message), collapse = "; ")
  } else {
    NA_character_
  }
  warnings <- if (!is.null(params$fit_warnings) && length(params$fit_warnings)) {
    paste(unique(as.character(params$fit_warnings)), collapse = "; ")
  } else {
    NA_character_
  }

  model_class <- if (is.null(fit)) {
    if (is.null(params$model_class)) NA_character_ else as.character(params$model_class)[1L]
  } else {
    class(fit)[1L]
  }
  converged <- NA
  singular <- NA
  rank_deficient <- NA

  if (inherits(fit, "merMod")) {
    conv_messages <- fit@optinfo$conv$lme4$messages
    opt_code <- fit@optinfo$conv$opt
    converged <- is.null(conv_messages) &&
      (is.null(opt_code) || all(as.numeric(opt_code) == 0))
    singular <- tryCatch(lme4::isSingular(fit), error = function(e) NA)
    rank_deficient <- anyNA(lme4::fixef(fit, add.dropped = TRUE))
  } else if (inherits(fit, "glmmTMB")) {
    converged <- isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess)
    rank_deficient <- anyNA(stats::coef(fit)$cond)
  } else if (inherits(fit, "glm")) {
    converged <- isTRUE(fit$converged)
    rank_deficient <- fit$rank < length(stats::coef(fit))
  } else if (inherits(fit, "lm")) {
    converged <- TRUE
    rank_deficient <- fit$rank < length(stats::coef(fit))
  } else if (is.null(fit)) {
    converged <- if (is.null(params$converged)) NA else params$converged
    singular <- if (is.null(params$singular)) NA else params$singular
    rank_deficient <- if (is.null(params$rank_deficient)) NA else params$rank_deficient
  }

  safe_number <- function(expr) {
    value <- tryCatch(suppressWarnings(as.numeric(expr)[1L]), error = function(e) NA_real_)
    if (length(value) == 0L || !is.finite(value)) NA_real_ else value
  }

  data.frame(
    feature = .mc_feature_name(object, fallback),
    fit_status = status,
    fit_message = message,
    fit_warnings = warnings,
    model_class = model_class,
    converged = as.logical(converged),
    singular = as.logical(singular),
    rank_deficient = as.logical(rank_deficient),
    n_obs = if (is.null(fit)) {
      if (is.null(params$n_obs)) NA_integer_ else as.integer(params$n_obs)[1L]
    } else {
      value <- tryCatch(stats::nobs(fit), error = function(e) NA_integer_)
      as.integer(value)[1L]
    },
    AIC = if (is.null(fit)) safe_number(params$AIC) else safe_number(stats::AIC(fit)),
    BIC = if (is.null(fit)) safe_number(params$BIC) else safe_number(stats::BIC(fit)),
    logLik = if (is.null(fit)) safe_number(params$logLik) else safe_number(stats::logLik(fit)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.mc_estimate_rows <- function(object, fallback = "feature") {
  groups <- .mc_group_levels(object)
  if (length(groups) == 0L) groups <- "all"
  grouped <- isTRUE(object@params$multiple_groups)
  period <- object@params$period
  if (is.null(period) || !is.finite(period)) period <- 24
  diagnostic <- .mc_model_diagnostics(object, fallback)

  rows <- lapply(groups, function(group) {
    suffix <- if (grouped) paste0("_", group) else ""
    mesor <- .mc_result_number(object, paste0("mesor", suffix))
    amplitude <- .mc_result_number(object, paste0("amps", suffix))
    phase <- .mc_result_number(object, paste0("phase", suffix))
    beta_cos <- .mc_result_number(object, paste0("beta_cos", suffix))
    beta_sin <- .mc_result_number(object, paste0("beta_sin", suffix))
    if (!is.finite(beta_cos) && is.finite(amplitude) && is.finite(phase)) {
      beta_cos <- amplitude * cos(phase)
    }
    if (!is.finite(beta_sin) && is.finite(amplitude) && is.finite(phase)) {
      beta_sin <- -amplitude * sin(phase)
    }

    data.frame(
      feature = .mc_feature_name(object, fallback),
      group = group,
      mesor = mesor,
      amplitude = amplitude,
      beta_cos = beta_cos,
      beta_sin = beta_sin,
      phase = phase,
      peak_time = if (is.finite(phase)) mixedcirc_detectpeaks(phase, period) else NA_real_,
      p_rhythm = .mc_result_number(object, paste0("p.value", suffix)),
      p_rhythm_global = .mc_result_number(object, "p.value_global"),
      p_differential_global = .mc_result_number(object, "p.value_differential_global"),
      period = as.numeric(period)[1L],
      analysis_type = object@type[1L],
      fit_status = diagnostic$fit_status,
      converged = diagnostic$converged,
      singular = diagnostic$singular,
      n_obs = diagnostic$n_obs,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.mc_contrast_rows <- function(object, fallback = "feature") {
  groups <- .mc_group_levels(object)
  if (length(groups) < 2L) return(data.frame())

  estimates <- .mc_estimate_rows(object, fallback)
  period <- estimates$period[1L]
  pairs <- utils::combn(groups, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    first <- estimates[match(pair[1L], estimates$group), , drop = FALSE]
    second <- estimates[match(pair[2L], estimates$group), , drop = FALSE]
    p_name <- paste0("p.value_", pair[1L], " - ", pair[2L])
    p_value <- .mc_result_number(object, p_name)
    if (!is.finite(p_value)) {
      reverse_name <- paste0("p.value_", pair[2L], " - ", pair[1L])
      p_value <- .mc_result_number(object, reverse_name)
    }
    phase_difference <- if (is.finite(first$phase) && is.finite(second$phase)) {
      .mc_wrap_radians(first$phase - second$phase)
    } else {
      NA_real_
    }
    peak_difference <- if (is.finite(first$peak_time) && is.finite(second$peak_time)) {
      .mc_wrap_period(first$peak_time - second$peak_time, period)
    } else {
      NA_real_
    }

    data.frame(
      feature = first$feature,
      contrast = paste0(pair[1L], " - ", pair[2L]),
      group_1 = pair[1L],
      group_2 = pair[2L],
      delta_mesor = first$mesor - second$mesor,
      delta_amplitude = first$amplitude - second$amplitude,
      amplitude_ratio = if (is.finite(second$amplitude) && second$amplitude > 0) {
        first$amplitude / second$amplitude
      } else {
        NA_real_
      },
      delta_beta_cos = first$beta_cos - second$beta_cos,
      delta_beta_sin = first$beta_sin - second$beta_sin,
      harmonic_distance = sqrt(
        (first$beta_cos - second$beta_cos)^2 +
          (first$beta_sin - second$beta_sin)^2
      ),
      delta_phase = phase_difference,
      delta_peak_time = peak_difference,
      p_differential = p_value,
      p_differential_global = first$p_differential_global,
      period = period,
      analysis_type = first$analysis_type,
      fit_status = first$fit_status,
      converged = first$converged,
      singular = first$singular,
      n_obs = first$n_obs,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Extract analysis-ready mixedcirc result tables
#'
#' Converts fitted objects into publication-oriented tables. Unlike the compact
#' one-row result stored in each fitted feature, these tables explicitly identify
#' features, groups, contrasts, effect sizes, multiplicity-adjusted p-values and
#' model diagnostics.
#'
#' @param x A [mixedcirc_fit-class] or [mixedcirc_fit_list-class] object.
#' @param type Table to return: `"wide"`, `"estimates"`, `"contrasts"`, or
#'   `"diagnostics"`.
#' @param adjust_method Method passed to [stats::p.adjust()]. Use `NULL` to omit
#'   adjusted p-values. Adjustment is performed across features separately for
#'   each group or contrast in long tables.
#' @param row.names,optional,... Standard [base::as.data.frame()] arguments;
#'   accepted for method compatibility and otherwise ignored.
#'
#' @return A data frame.
#' @export
mixedcirc_results <- function(x,
                              type = c("wide", "estimates", "contrasts", "diagnostics"),
                              adjust_method = "BH") {
  type <- match.arg(type)
  objects <- .mc_fit_items(x)
  if (length(objects) == 0L) return(data.frame())
  if (!is.null(adjust_method)) {
    if (!is.character(adjust_method) || length(adjust_method) != 1L ||
        is.na(adjust_method) || !adjust_method %in% stats::p.adjust.methods) {
      stop("adjust_method must be NULL or a method accepted by p.adjust", call. = FALSE)
    }
  }

  if (type == "diagnostics") {
    return(do.call(rbind, lapply(seq_along(objects), function(i) {
      .mc_model_diagnostics(objects[[i]], paste0("feature_", i))
    })))
  }

  if (type == "wide") {
    rows <- lapply(seq_along(objects), function(i) {
      object <- objects[[i]]
      diagnostic <- .mc_model_diagnostics(object, paste0("feature_", i))
      result <- as.data.frame(object@results, check.names = FALSE)
      if (nrow(result) == 0L) result <- data.frame(.placeholder = NA_real_)
      result$feature <- diagnostic$feature
      result$fit_status <- diagnostic$fit_status
      result$converged <- diagnostic$converged
      result$singular <- diagnostic$singular
      result$n_obs <- diagnostic$n_obs
      result <- result[, c("feature", setdiff(names(result), "feature")), drop = FALSE]
      result$.placeholder <- NULL
      result
    })
    out <- .mc_bind_rows_union(rows)
    if (!is.null(adjust_method)) {
      p_columns <- grep("^p\\.value", names(out), value = TRUE)
      for (column in p_columns) {
        q_column <- sub("^p\\.value", "q.value", column)
        out[[q_column]] <- .mc_adjust_column(as.numeric(out[[column]]), adjust_method)
      }
    }
    return(out)
  }

  if (type == "estimates") {
    out <- do.call(rbind, lapply(seq_along(objects), function(i) {
      .mc_estimate_rows(objects[[i]], paste0("feature_", i))
    }))
    rownames(out) <- NULL
    if (!is.null(adjust_method)) {
      out$q_rhythm <- NA_real_
      for (group in unique(out$group)) {
        take <- out$group == group
        out$q_rhythm[take] <- .mc_adjust_column(out$p_rhythm[take], adjust_method)
      }
      global_by_feature <- !duplicated(out$feature)
      global_q <- .mc_adjust_column(out$p_rhythm_global[global_by_feature], adjust_method)
      out$q_rhythm_global <- global_q[match(out$feature, out$feature[global_by_feature])]
      diff_global_q <- .mc_adjust_column(
        out$p_differential_global[global_by_feature], adjust_method
      )
      out$q_differential_global <- diff_global_q[
        match(out$feature, out$feature[global_by_feature])
      ]
    }
    return(out)
  }

  rows <- lapply(seq_along(objects), function(i) {
    .mc_contrast_rows(objects[[i]], paste0("feature_", i))
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) {
    return(data.frame(
      feature = character(), contrast = character(), group_1 = character(),
      group_2 = character(), delta_mesor = numeric(),
      delta_amplitude = numeric(), amplitude_ratio = numeric(),
      delta_beta_cos = numeric(), delta_beta_sin = numeric(),
      harmonic_distance = numeric(), delta_phase = numeric(),
      delta_peak_time = numeric(), p_differential = numeric(),
      q_differential = numeric(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(adjust_method)) {
    out$q_differential <- NA_real_
    for (contrast in unique(out$contrast)) {
      take <- out$contrast == contrast
      out$q_differential[take] <- .mc_adjust_column(
        out$p_differential[take], adjust_method
      )
    }
    first_by_feature <- !duplicated(out$feature)
    global_q <- .mc_adjust_column(
      out$p_differential_global[first_by_feature], adjust_method
    )
    out$q_differential_global <- global_q[
      match(out$feature, out$feature[first_by_feature])
    ]
  }
  out
}

#' Extract pairwise circadian contrasts
#'
#' Convenience wrapper for `mixedcirc_results(x, type = "contrasts")`.
#'
#' @inheritParams mixedcirc_results
#' @return A data frame with one row per feature and group contrast.
#' @export
mixedcirc_contrasts <- function(x, adjust_method = "BH") {
  mixedcirc_results(x, type = "contrasts", adjust_method = adjust_method)
}

#' Classify changes in rhythmicity between groups
#'
#' Uses direct between-group harmonic tests together with within-group
#' rhythmicity tests. A feature is never labelled as changed merely because it
#' is significant in one group and not another.
#'
#' @param x A fitted mixedcirc object.
#' @param alpha FDR threshold.
#' @param min_amplitude Minimum amplitude required to call a group rhythmic.
#' @param min_amplitude_change Minimum absolute amplitude difference considered
#'   biologically relevant.
#' @param min_phase_shift Minimum absolute peak-time shift considered
#'   biologically relevant, in the original time units.
#' @param min_harmonic_distance Minimum Euclidean distance between the two
#'   group harmonic coefficient vectors.
#' @param adjust_method Multiple-testing method.
#'
#' @return The pairwise contrast table with group rhythmicity indicators and a
#'   conservative `rhythm_class` column.
#' @export
mixedcirc_classify <- function(x,
                               alpha = 0.05,
                               min_amplitude = 0,
                               min_amplitude_change = 0,
                               min_phase_shift = 0,
                               min_harmonic_distance = 0,
                               adjust_method = "BH") {
  numeric_args <- list(
    alpha = alpha,
    min_amplitude = min_amplitude,
    min_amplitude_change = min_amplitude_change,
    min_phase_shift = min_phase_shift,
    min_harmonic_distance = min_harmonic_distance
  )
  bad <- !vapply(numeric_args, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }, logical(1))
  if (any(bad)) stop(names(numeric_args)[which(bad)[1L]], " must be one finite number")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between zero and one")
  if (any(unlist(numeric_args[-1L]) < 0)) {
    stop("effect-size thresholds must be non-negative")
  }

  estimates <- mixedcirc_results(x, "estimates", adjust_method)
  contrasts <- mixedcirc_results(x, "contrasts", adjust_method)
  if (nrow(contrasts) == 0L) return(contrasts)

  key <- paste(estimates$feature, estimates$group, sep = "\r")
  key_1 <- paste(contrasts$feature, contrasts$group_1, sep = "\r")
  key_2 <- paste(contrasts$feature, contrasts$group_2, sep = "\r")
  first <- estimates[match(key_1, key), , drop = FALSE]
  second <- estimates[match(key_2, key), , drop = FALSE]

  contrasts$group_1_amplitude <- first$amplitude
  contrasts$group_2_amplitude <- second$amplitude
  rhythm_column <- if ("q_rhythm" %in% names(estimates)) "q_rhythm" else "p_rhythm"
  differential_column <- if ("q_differential" %in% names(contrasts)) {
    "q_differential"
  } else {
    "p_differential"
  }
  group_1_significance <- first[[rhythm_column]]
  group_2_significance <- second[[rhythm_column]]
  contrasts$group_1_q_rhythm <- group_1_significance
  contrasts$group_2_q_rhythm <- group_2_significance
  contrasts$group_1_rhythmic <- is.finite(group_1_significance) &
    group_1_significance < alpha & is.finite(first$amplitude) &
    first$amplitude >= min_amplitude
  contrasts$group_2_rhythmic <- is.finite(group_2_significance) &
    group_2_significance < alpha & is.finite(second$amplitude) &
    second$amplitude >= min_amplitude

  thresholds_enabled <- any(c(
    min_amplitude_change,
    min_phase_shift,
    min_harmonic_distance
  ) > 0)
  effect_relevant <- rep(TRUE, nrow(contrasts))
  if (thresholds_enabled) {
    effect_relevant <-
      (is.finite(contrasts$delta_amplitude) &
         abs(contrasts$delta_amplitude) >= min_amplitude_change &
         min_amplitude_change > 0) |
      (is.finite(contrasts$delta_peak_time) &
         abs(contrasts$delta_peak_time) >= min_phase_shift &
         min_phase_shift > 0) |
      (is.finite(contrasts$harmonic_distance) &
         contrasts$harmonic_distance >= min_harmonic_distance &
         min_harmonic_distance > 0)
  }
  contrasts$effect_relevant <- effect_relevant
  differential_significance <- contrasts[[differential_column]]
  differential <- is.finite(differential_significance) &
    differential_significance < alpha & effect_relevant

  class <- rep("uncertain", nrow(contrasts))
  both <- contrasts$group_1_rhythmic & contrasts$group_2_rhythmic
  neither <- !contrasts$group_1_rhythmic & !contrasts$group_2_rhythmic
  class[neither & !differential] <- "arrhythmic"
  class[both & !differential] <- "no_detected_change"
  class[both & differential] <- "altered_rhythm"
  class[contrasts$group_1_rhythmic & !contrasts$group_2_rhythmic & differential] <-
    paste0("loss_in_", contrasts$group_2[
      contrasts$group_1_rhythmic & !contrasts$group_2_rhythmic & differential
    ])
  class[!contrasts$group_1_rhythmic & contrasts$group_2_rhythmic & differential] <-
    paste0("gain_in_", contrasts$group_2[
      !contrasts$group_1_rhythmic & contrasts$group_2_rhythmic & differential
    ])
  contrasts$rhythm_class <- class
  contrasts
}

#' @rdname mixedcirc_results
#' @export
methods::setMethod("as.data.frame", "mixedcirc_fit", function(x,
                                                               row.names = NULL,
                                                               optional = FALSE,
                                                               ...) {
  mixedcirc_results(x, type = "wide")
})

#' @rdname mixedcirc_results
#' @export
methods::setMethod("as.data.frame", "mixedcirc_fit_list", function(x,
                                                                    row.names = NULL,
                                                                    optional = FALSE,
                                                                    ...) {
  mixedcirc_results(x, type = "wide")
})
