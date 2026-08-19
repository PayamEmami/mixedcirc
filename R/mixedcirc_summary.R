# Publication-level analysis summary ---------------------------------------

#' Summarize a mixedcirc analysis
#'
#' Produces a compact analysis-level summary of fit quality, rhythmic features,
#' direct differential-rhythm tests, and biological change classes. This is
#' intended for quality control, reports, and compact result tables.
#'
#' @param x A [mixedcirc_fit-class] or [mixedcirc_fit_list-class] object.
#' @param alpha Significance threshold applied after `adjust_method` correction.
#' @param adjust_method A method accepted by [stats::p.adjust()], or `NULL` to
#'   summarize raw p-values.
#' @param min_amplitude,min_amplitude_change,min_phase_shift,min_harmonic_distance
#'   Effect-size thresholds passed to
#'   [mixedcirc_classify()].
#'
#' @return An object of class `mixedcirc_summary`, containing `overview`,
#'   `groups`, `contrasts`, `classes`, and `diagnostics` tables.
#'
#' @details Group medians are calculated among features that meet the selected
#' rhythmicity and amplitude thresholds. Peak time uses a circular median so
#' that values on opposite sides of the cycle boundary remain close.
#' @export
mixedcirc_summary <- function(x,
                              alpha = 0.05,
                              adjust_method = "BH",
                              min_amplitude = 0,
                              min_amplitude_change = 0,
                              min_phase_shift = 0,
                              min_harmonic_distance = 0) {
  circular_median <- function(values, period) {
    values <- values[is.finite(values)] %% period
    if (!length(values)) return(NA_real_)
    distance_sum <- vapply(values, function(candidate) {
      distances <- ((values - candidate + period / 2) %% period) - period / 2
      sum(abs(distances))
    }, numeric(1))
    values[which.min(distance_sum)]
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be one number between 0 and 1", call. = FALSE)
  }
  if (!is.null(adjust_method) &&
      (!is.character(adjust_method) || length(adjust_method) != 1L ||
       is.na(adjust_method) || !adjust_method %in% stats::p.adjust.methods)) {
    stop("adjust_method must be NULL or a method accepted by p.adjust", call. = FALSE)
  }

  estimates <- mixedcirc_results(x, "estimates", adjust_method = adjust_method)
  diagnostics <- mixedcirc_results(x, "diagnostics", adjust_method = adjust_method)
  contrasts <- mixedcirc_contrasts(x, adjust_method = adjust_method)
  significance_column <- if (is.null(adjust_method)) "p_rhythm" else "q_rhythm"
  differential_column <- if (is.null(adjust_method)) {
    "p_differential"
  } else {
    "q_differential"
  }
  global_differential_column <- if (is.null(adjust_method)) {
    "p_differential_global"
  } else {
    "q_differential_global"
  }

  group_rows <- lapply(unique(estimates$group), function(group) {
    current <- estimates[estimates$group == group, , drop = FALSE]
    significant <- is.finite(current[[significance_column]]) &
      current[[significance_column]] < alpha &
      is.finite(current$amplitude) & current$amplitude >= min_amplitude
    data.frame(
      group = group,
      tested_features = sum(is.finite(current$p_rhythm)),
      rhythmic_features = sum(significant),
      rhythmic_proportion = if (sum(is.finite(current$p_rhythm)) > 0L) {
        sum(significant) / sum(is.finite(current$p_rhythm))
      } else {
        NA_real_
      },
      median_amplitude = if (any(significant)) {
        stats::median(current$amplitude[significant], na.rm = TRUE)
      } else {
        NA_real_
      },
      median_peak_time = circular_median(
        current$peak_time[significant],
        current$period[1L]
      ),
      stringsAsFactors = FALSE
    )
  })
  group_summary <- do.call(rbind, group_rows)

  if (nrow(contrasts)) {
    contrast_rows <- lapply(unique(contrasts$contrast), function(contrast) {
      current <- contrasts[contrasts$contrast == contrast, , drop = FALSE]
      significant <- is.finite(current[[differential_column]]) &
        current[[differential_column]] < alpha
      data.frame(
        contrast = contrast,
        tested_features = sum(is.finite(current$p_differential)),
        differential_features = sum(significant),
        differential_proportion = if (sum(is.finite(current$p_differential)) > 0L) {
          sum(significant) / sum(is.finite(current$p_differential))
        } else {
          NA_real_
        },
        median_harmonic_distance = stats::median(
          current$harmonic_distance, na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    })
    contrast_summary <- do.call(rbind, contrast_rows)
    classes <- mixedcirc_classify(
      x,
      alpha = alpha,
      min_amplitude = min_amplitude,
      min_amplitude_change = min_amplitude_change,
      min_phase_shift = min_phase_shift,
      min_harmonic_distance = min_harmonic_distance,
      adjust_method = adjust_method
    )
    class_summary <- as.data.frame(table(
      contrast = classes$contrast,
      rhythm_class = classes$rhythm_class,
      useNA = "ifany"
    ), stringsAsFactors = FALSE)
    class_summary <- class_summary[class_summary$Freq > 0L, , drop = FALSE]
  } else {
    contrast_summary <- data.frame(
      contrast = character(), tested_features = integer(),
      differential_features = integer(), differential_proportion = numeric(),
      median_harmonic_distance = numeric(), stringsAsFactors = FALSE
    )
    class_summary <- data.frame(
      contrast = character(), rhythm_class = character(), Freq = integer(),
      stringsAsFactors = FALSE
    )
  }

  first_estimate <- !duplicated(estimates$feature)
  omnibus_significant <- if (global_differential_column %in% names(estimates)) {
    values <- estimates[[global_differential_column]][first_estimate]
    sum(is.finite(values) & values < alpha)
  } else {
    0L
  }
  overview <- data.frame(
    features = length(unique(estimates$feature)),
    groups = length(unique(estimates$group)),
    period = if (length(unique(estimates$period)) == 1L) unique(estimates$period) else NA_real_,
    analysis_type = paste(unique(estimates$analysis_type), collapse = ", "),
    successful_fits = sum(diagnostics$fit_status == "ok"),
    failed_fits = sum(diagnostics$fit_status != "ok"),
    converged_fits = sum(diagnostics$converged %in% TRUE),
    singular_fits = sum(diagnostics$singular %in% TRUE),
    omnibus_differential_features = omnibus_significant,
    alpha = alpha,
    adjustment = if (is.null(adjust_method)) "none" else adjust_method,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      overview = overview,
      groups = group_summary,
      contrasts = contrast_summary,
      classes = class_summary,
      diagnostics = diagnostics
    ),
    class = "mixedcirc_summary"
  )
}

#' Print a mixedcirc analysis summary
#'
#' @param x An object returned by [mixedcirc_summary()].
#' @param ... Additional arguments; currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.mixedcirc_summary <- function(x, ...) {
  cat("mixedcirc analysis summary\n")
  cat("--------------------------\n")
  cat("Features:", x$overview$features, "\n")
  cat("Groups:", x$overview$groups, "\n")
  cat("Analysis:", x$overview$analysis_type, "\n")
  cat(
    "Successful fits:", x$overview$successful_fits,
    "of", x$overview$features, "\n"
  )
  if (x$overview$groups > 1L) {
    cat(
      "Omnibus differential features:",
      x$overview$omnibus_differential_features,
      "at", x$overview$adjustment, "<", x$overview$alpha, "\n"
    )
  }
  cat("\nGroup rhythmicity:\n")
  print(x$groups, row.names = FALSE)
  if (nrow(x$contrasts)) {
    cat("\nDirect differential-rhythm tests:\n")
    print(x$contrasts, row.names = FALSE)
  }
  invisible(x)
}

#' @param object A [mixedcirc_fit-class] object.
#' @param ... Additional arguments passed to [mixedcirc_summary()].
#' @rdname mixedcirc_summary
#' @export
methods::setMethod("summary", "mixedcirc_fit", function(object, ...) {
  mixedcirc_summary(object, ...)
})

#' @param object A [mixedcirc_fit_list-class] object.
#' @param ... Additional arguments passed to [mixedcirc_summary()].
#' @rdname mixedcirc_summary
#' @export
methods::setMethod("summary", "mixedcirc_fit_list", function(object, ...) {
  mixedcirc_summary(object, ...)
})
