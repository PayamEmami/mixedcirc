# Publication-oriented summaries and plots ----------------------------------

.mc_validate_features <- function(available, features = NULL) {
  if (is.null(features)) return(available)
  if (!is.character(features) || anyNA(features) || any(!nzchar(features))) {
    stop("features must be NULL or a character vector", call. = FALSE)
  }
  missing_features <- setdiff(features, available)
  if (length(missing_features)) {
    stop(
      "Unknown features: ", paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }
  unique(features)
}

.mc_scale_within <- function(value, groups, method) {
  split_index <- split(seq_along(value), groups, drop = TRUE)
  out <- value
  for (index in split_index) {
    current <- value[index]
    if (method == "raw") {
      out[index] <- current
    } else if (method == "center") {
      out[index] <- current - mean(current, na.rm = TRUE)
    } else {
      spread <- stats::sd(current, na.rm = TRUE)
      if (!is.finite(spread) || spread == 0) spread <- 1
      out[index] <- (current - mean(current, na.rm = TRUE)) / spread
    }
  }
  out
}

#' Extract subject-specific circadian rhythm estimates
#'
#' Combines population harmonic coefficients with fitted subject-level random
#' effects. This makes heterogeneity in mesor, amplitude and peak time directly
#' available for longitudinal analyses.
#'
#' @param x A [mixedcirc_fit-class] or [mixedcirc_fit_list-class] object with
#'   retained `lmer` models.
#' @param features Optional feature names to include.
#'
#' @return A data frame with one row per feature and subject.
#' @export
mixedcirc_subject_rhythms <- function(x, features = NULL) {
  objects <- .mc_fit_items(x)
  available <- vapply(seq_along(objects), function(index) {
    .mc_feature_name(objects[[index]], paste0("feature_", index))
  }, character(1))
  features <- .mc_validate_features(available, features)
  objects <- objects[match(features, available)]

  rows <- lapply(seq_along(objects), function(index) {
    object <- objects[[index]]
    fit <- object@fit
    feature <- features[index]
    if (!inherits(fit, "merMod")) return(NULL)

    id_var <- object@params$id
    if (is.null(id_var) || !id_var %in% names(object@exp_design)) return(NULL)
    random_effects <- lme4::ranef(fit)
    if (!id_var %in% names(random_effects)) return(NULL)
    subject_effects <- as.data.frame(random_effects[[id_var]])
    subject_effects$subject <- rownames(random_effects[[id_var]])
    rownames(subject_effects) <- NULL

    design <- object@exp_design
    subject_group <- if (isTRUE(object@params$multiple_groups)) {
      group_var <- object@params$group
      mapping <- unique(data.frame(
        subject = as.character(design[[id_var]]),
        group = as.character(design[[group_var]]),
        stringsAsFactors = FALSE
      ))
      duplicated_subject <- duplicated(mapping$subject) |
        duplicated(mapping$subject, fromLast = TRUE)
      if (any(duplicated_subject)) {
        stop("Each subject must belong to one group", call. = FALSE)
      }
      mapping$group[match(subject_effects$subject, mapping$subject)]
    } else {
      rep("all", nrow(subject_effects))
    }

    estimates <- .mc_estimate_rows(object, feature)
    population <- estimates[match(subject_group, estimates$group), , drop = FALSE]
    column_value <- function(candidates) {
      matched <- candidates[candidates %in% names(subject_effects)]
      if (!length(matched)) return(rep(0, nrow(subject_effects)))
      as.numeric(subject_effects[[matched[1L]]])
    }
    random_intercept <- column_value("(Intercept)")
    if (isTRUE(object@params$RRBS)) {
      random_beta_cos <- column_value(c(
        "scaler:inphase", "inphase:scaler"
      ))
      random_beta_sin <- column_value(c(
        "scaler:outphase", "outphase:scaler"
      ))
      subject_mesor <- population$mesor
    } else {
      random_beta_cos <- column_value("inphase")
      random_beta_sin <- column_value("outphase")
      subject_mesor <- population$mesor + random_intercept
    }

    beta_cos <- population$beta_cos + random_beta_cos
    beta_sin <- population$beta_sin + random_beta_sin
    phase <- mapply(calc_phase, beta_cos, beta_sin)
    data.frame(
      feature = feature,
      subject = subject_effects$subject,
      group = subject_group,
      population_mesor = population$mesor,
      population_amplitude = population$amplitude,
      random_mesor = random_intercept,
      random_beta_cos = random_beta_cos,
      random_beta_sin = random_beta_sin,
      subject_mesor = subject_mesor,
      subject_amplitude = sqrt(beta_cos^2 + beta_sin^2),
      subject_phase = as.numeric(phase),
      subject_peak_time = vapply(
        phase,
        mixedcirc_detectpeaks,
        numeric(1),
        period = population$period[1L]
      ),
      period = population$period,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  rows <- rows[vapply(rows, function(value) !is.null(value), logical(1))]
  if (!length(rows)) {
    stop(
      "No retained mixed-effects models with subject random effects were found",
      call. = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot a phase-ordered fitted rhythm heatmap
#'
#' @param x A fitted mixedcirc object with retained models.
#' @param features Optional feature names to include.
#' @param groups Optional group levels to include.
#' @param scale Row scaling: `"zscore"`, `"center"`, or `"raw"`.
#' @param order_by Order rows by fitted `"peak"` or input `"feature"` order.
#' @param show_feature_labels Logical controlling y-axis feature labels. The
#'   default (`NULL`) shows them for at most 50 features and hides them for
#'   larger panels.
#' @param npoints Number of prediction points per cycle.
#' @param fixed_values,expand_vars,numeric_method,factor_method Passed to
#'   [mixedcirc_interpolate()].
#'
#' @return A ggplot object.
#' @export
mixedcirc_plot_heatmap <- function(x,
                                   features = NULL,
                                   groups = NULL,
                                   scale = c("zscore", "center", "raw"),
                                   order_by = c("peak", "feature"),
                                   show_feature_labels = NULL,
                                   npoints = 96,
                                   fixed_values = list(),
                                   expand_vars = NULL,
                                   numeric_method = "mean",
                                   factor_method = "average_observed") {
  scale <- match.arg(scale)
  order_by <- match.arg(order_by)
  available <- vapply(seq_along(.mc_fit_items(x)), function(index) {
    .mc_feature_name(.mc_fit_items(x)[[index]], paste0("feature_", index))
  }, character(1))
  features <- .mc_validate_features(available, features)
  if (is.null(show_feature_labels)) {
    show_feature_labels <- length(features) <= 50L
  } else if (!is.logical(show_feature_labels) ||
             length(show_feature_labels) != 1L || is.na(show_feature_labels)) {
    stop("show_feature_labels must be NULL, TRUE, or FALSE", call. = FALSE)
  }
  selected <- if (methods::is(x, "mixedcirc_fit")) x else {
    methods::new("mixedcirc_fit_list", results = x@results[match(features, available)])
  }
  predictions <- mixedcirc_interpolate(
    selected,
    npoints = npoints,
    fixed_values = fixed_values,
    expand_vars = expand_vars,
    numeric_method = numeric_method,
    factor_method = factor_method
  )
  if (!"feature" %in% names(predictions)) predictions$feature <- features[1L]
  spec <- .get_plot_spec(if (methods::is(selected, "mixedcirc_fit")) selected else selected[[1L]])
  time_var <- spec$time_var
  group_var <- spec$group_var
  predictions$.group <- if (!is.null(group_var) && group_var %in% names(predictions)) {
    as.character(predictions[[group_var]])
  } else {
    "all"
  }
  group_order <- unique(predictions$.group)
  if (!is.null(groups)) {
    if (!is.character(groups) || anyNA(groups)) stop("groups must be a character vector")
    missing_groups <- setdiff(groups, unique(predictions$.group))
    if (length(missing_groups)) stop("Unknown groups: ", paste(missing_groups, collapse = ", "))
    predictions <- predictions[predictions$.group %in% groups, , drop = FALSE]
  }
  predictions$.value <- .mc_scale_within(
    predictions$Y.hat,
    interaction(predictions$feature, predictions$.group, drop = TRUE),
    scale
  )

  feature_order <- features
  if (order_by == "peak") {
    estimates <- mixedcirc_results(selected, "estimates", adjust_method = NULL)
    estimates <- estimates[estimates$group %in% unique(predictions$.group), , drop = FALSE]
    first_peak <- estimates[!duplicated(estimates$feature), c("feature", "peak_time")]
    feature_order <- first_peak$feature[order(first_peak$peak_time, na.last = TRUE)]
  }
  predictions$feature <- factor(predictions$feature, levels = rev(feature_order))
  predictions$.group <- factor(predictions$.group, levels = group_order)

  plot <- ggplot2::ggplot(
    predictions,
    ggplot2::aes(x = .data[[time_var]], y = .data$feature, fill = .data$.value)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = switch(scale, zscore = "Row z-score", center = "Centered", raw = "Fitted")
    ) +
    ggplot2::labs(x = time_var, y = NULL, title = "Phase-ordered circadian profiles") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = if (isTRUE(show_feature_labels)) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.y = if (isTRUE(show_feature_labels)) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      }
    )
  if (length(unique(predictions$.group)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(~.group, ncol = 1)
  }
  plot
}

#' Plot differential-rhythm effect sizes and significance
#'
#' @param x A grouped fitted mixedcirc object.
#' @param contrast Optional contrast label such as `"control - treatment"`.
#' @param effect Effect size shown on the x-axis.
#' @param alpha Adjusted-significance threshold.
#' @param label_n Number of the most significant features to label per plot.
#' @param ... Effect-size thresholds passed to [mixedcirc_classify()].
#'
#' @return A ggplot object.
#' @export
mixedcirc_plot_differential <- function(x,
                                        contrast = NULL,
                                        effect = c(
                                          "delta_amplitude", "delta_peak_time",
                                          "delta_mesor", "harmonic_distance"
                                        ),
                                        alpha = 0.05,
                                        label_n = 10,
                                        ...) {
  effect <- match.arg(effect)
  if (!is.numeric(label_n) || length(label_n) != 1L || label_n < 0 || label_n %% 1 != 0) {
    stop("label_n must be a non-negative integer", call. = FALSE)
  }
  dat <- mixedcirc_classify(x, alpha = alpha, ...)
  if (!nrow(dat)) stop("x contains no between-group contrasts", call. = FALSE)
  if (!is.null(contrast)) {
    if (!is.character(contrast) || length(contrast) != 1L || !contrast %in% dat$contrast) {
      stop("contrast must identify one available contrast", call. = FALSE)
    }
    dat <- dat[dat$contrast == contrast, , drop = FALSE]
  }
  dat$.significance <- ifelse(
    is.finite(dat$q_differential),
    -log10(pmax(dat$q_differential, .Machine$double.xmin)),
    NA_real_
  )
  label_data <- dat[order(dat$q_differential, na.last = NA), , drop = FALSE]
  if (nrow(label_data) > label_n) label_data <- label_data[seq_len(label_n), , drop = FALSE]
  threshold <- -log10(alpha)
  plot <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = .data[[effect]], y = .data$.significance, color = .data$rhythm_class)
  ) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2, color = "grey45") +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(label = .data$feature),
      check_overlap = TRUE,
      vjust = -0.5,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = gsub("_", " ", effect),
      y = expression(-log[10]("differential rhythm FDR")),
      color = "Rhythm class",
      title = "Differential circadian rhythms"
    ) +
    ggplot2::scale_color_discrete(labels = function(value) gsub("_", " ", value)) +
    ggplot2::theme_classic(base_size = 11)
  if (length(unique(dat$contrast)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(~contrast, scales = "free_x")
  }
  plot
}

#' Plot amplitude and peak time in polar coordinates
#'
#' @param x A fitted mixedcirc object.
#' @param features Optional feature names to include.
#' @param groups Optional group levels to include.
#' @param alpha FDR threshold used to distinguish rhythmic estimates.
#' @param label_n Number of largest-amplitude rhythmic features to label.
#'
#' @return A ggplot object.
#' @export
mixedcirc_plot_polar <- function(x,
                                 features = NULL,
                                 groups = NULL,
                                 alpha = 0.05,
                                 label_n = 8) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be one number between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(label_n) || length(label_n) != 1L || !is.finite(label_n) ||
      label_n < 0 || label_n %% 1 != 0) {
    stop("label_n must be a non-negative integer", call. = FALSE)
  }
  dat <- mixedcirc_results(x, "estimates")
  features <- .mc_validate_features(unique(dat$feature), features)
  dat <- dat[dat$feature %in% features, , drop = FALSE]
  if (!is.null(groups)) {
    if (!all(groups %in% dat$group)) stop("groups contains an unknown group")
    dat <- dat[dat$group %in% groups, , drop = FALSE]
  }
  rhythmic <- is.finite(dat$q_rhythm) & dat$q_rhythm < alpha
  dat$rhythmic <- factor(ifelse(rhythmic, "Yes", "No"), levels = c("Yes", "No"))
  labels <- dat[rhythmic & is.finite(dat$amplitude), , drop = FALSE]
  labels <- labels[order(labels$amplitude, decreasing = TRUE), , drop = FALSE]
  labels <- labels[!duplicated(labels$feature), , drop = FALSE]
  if (nrow(labels) > label_n) labels <- labels[seq_len(label_n), , drop = FALSE]
  period <- unique(dat$period)
  if (length(period) != 1L) stop("x contains more than one fitted period")
  breaks <- seq(0, period, length.out = 9L)[-9L]

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$peak_time, y = .data$amplitude,
      color = .data$group, alpha = .data$rhythmic
    )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(label = .data$feature),
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    ggplot2::scale_x_continuous(limits = c(0, period), breaks = breaks) +
    ggplot2::scale_alpha_manual(values = c(Yes = 0.9, No = 0.2)) +
    ggplot2::labs(
      x = paste0("Peak time (", period, "-unit cycle)"),
      y = "Amplitude", color = "Group", alpha = "Rhythmic",
      title = "Circadian phase and amplitude landscape"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot bootstrap confidence intervals
#'
#' @param x A [mixedcirc_boot-class] object.
#' @param parameters Optional result row names to include.
#'
#' Circular phase intervals that cross the `-pi`/`pi` boundary are drawn on a
#' continuous local scale around their point estimate rather than across almost
#' the entire horizontal axis.
#'
#' @return A ggplot object.
#' @export
mixedcirc_plot_bootstrap <- function(x, parameters = NULL) {
  if (!methods::is(x, "mixedcirc_boot")) {
    stop("x must be a mixedcirc_boot object", call. = FALSE)
  }
  dat <- as.data.frame(x)
  dat$parameter <- rownames(dat)
  required <- c("boot.lower", "boot.upper", "boot.estimate")
  if (!all(required %in% names(dat))) stop("bootstrap result lacks confidence intervals")
  if (!is.null(parameters)) {
    if (!is.character(parameters) || !all(parameters %in% dat$parameter)) {
      stop("parameters contains an unknown result name", call. = FALSE)
    }
    dat <- dat[dat$parameter %in% parameters, , drop = FALSE]
  }
  dat$.plot_lower <- dat$boot.lower
  dat$.plot_upper <- dat$boot.upper
  circular <- grepl("^phase", dat$parameter)
  if (any(circular)) {
    lower_deviation <- .mc_wrap_radians(
      dat$boot.lower[circular] - dat$boot.estimate[circular]
    )
    upper_deviation <- .mc_wrap_radians(
      dat$boot.upper[circular] - dat$boot.estimate[circular]
    )
    dat$.plot_lower[circular] <- dat$boot.estimate[circular] +
      pmin(lower_deviation, upper_deviation)
    dat$.plot_upper[circular] <- dat$boot.estimate[circular] +
      pmax(lower_deviation, upper_deviation)
  }
  dat$parameter <- factor(dat$parameter, levels = rev(dat$parameter))
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$boot.estimate, y = .data$parameter)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$.plot_lower, xend = .data$.plot_upper,
        y = .data$parameter, yend = .data$parameter
      ),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2, color = "#2166AC") +
    ggplot2::labs(
      x = paste0(round(100 * x@metadata$conf_level), "% bootstrap interval"),
      y = NULL,
      title = "Bootstrap uncertainty"
    ) +
    ggplot2::theme_classic(base_size = 11)
}

#' Plot the temporal sampling design
#'
#' @param x A fitted mixedcirc object.
#' @param by_subject Logical. If `FALSE` (default), shows the number of
#'   observations at each time and group combination. If `TRUE`, shows one row
#'   per subject, making the repeated-measures structure of the design
#'   directly visible.
#'
#' @return A ggplot object showing the temporal sampling design.
#' @export
mixedcirc_plot_design <- function(x, by_subject = FALSE) {
  if (!is.logical(by_subject) || length(by_subject) != 1L || is.na(by_subject)) {
    stop("by_subject must be TRUE or FALSE", call. = FALSE)
  }
  objects <- .mc_fit_items(x)
  if (!length(objects)) stop("x contains no fits", call. = FALSE)
  object <- objects[[1L]]
  dat <- object@exp_design
  spec <- .get_plot_spec(object)
  if (isTRUE(spec$RRBS) && "scaler" %in% names(dat)) {
    dat <- dat[dat$scaler == 1, , drop = FALSE]
  }
  dat$.group <- if (!is.null(spec$group_var)) as.character(dat[[spec$group_var]]) else "all"

  if (by_subject) {
    if (is.null(spec$id_var) || !spec$id_var %in% names(dat)) {
      stop("by_subject = TRUE requires a subject/id column identifiable in the fitted object", call. = FALSE)
    }
    dat$.subject <- as.character(dat[[spec$id_var]])
    row_order <- order(dat$.group, dat$.subject)
    subject_order <- unique(dat$.subject[row_order])
    dat$.subject <- factor(dat$.subject, levels = subject_order)
    plot_obj <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = .data[[spec$time_var]], y = .data$.subject, color = .data$.group)
    ) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::labs(
        x = spec$time_var, y = "Subject", color = "Group",
        title = "Temporal sampling design"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(legend.position = "right")
    return(plot_obj)
  }

  counts <- stats::aggregate(
    rep(1L, nrow(dat)),
    list(time = dat[[spec$time_var]], group = dat$.group),
    length
  )
  names(counts)[3L] <- "n"
  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data$time, y = .data$group)
  ) +
    ggplot2::geom_line(ggplot2::aes(group = .data$group), linewidth = 0.4, color = "grey70") +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$n, color = .data$group),
      alpha = 0.85
    ) +
    ggplot2::scale_size_area(max_size = 8, breaks = unique(sort(counts$n))) +
    ggplot2::labs(
      x = spec$time_var, y = NULL, size = "Observations", color = "Group",
      title = "Temporal sampling design"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "right")
}

#' Plot model residual diagnostics across features
#'
#' @param x A fitted mixedcirc object.
#' @param features Optional feature names to include.
#' @param type `"residual_fitted"` or `"qq"`.
#'
#' @return A faceted ggplot object.
#' @export
mixedcirc_plot_diagnostics <- function(x,
                                       features = NULL,
                                       type = c("residual_fitted", "qq")) {
  type <- match.arg(type)
  objects <- .mc_fit_items(x)
  available <- vapply(seq_along(objects), function(index) {
    .mc_feature_name(objects[[index]], paste0("feature_", index))
  }, character(1))
  features <- .mc_validate_features(available, features)
  objects <- objects[match(features, available)]
  dat <- do.call(rbind, lapply(seq_along(objects), function(index) {
    current <- objects[[index]]@exp_design
    if (!all(c(".fitted", ".residual") %in% names(current))) return(NULL)
    data.frame(
      feature = features[index],
      fitted = as.numeric(current$.fitted),
      residual = as.numeric(current$.residual),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(dat)) {
    stop("No stored residual diagnostics are available", call. = FALSE)
  }
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (!nrow(dat)) stop("No stored residual diagnostics are available", call. = FALSE)

  if (type == "residual_fitted") {
    return(
      ggplot2::ggplot(dat, ggplot2::aes(x = .data$fitted, y = .data$residual)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        ggplot2::geom_point(alpha = 0.45, size = 1) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5) +
        ggplot2::facet_wrap(~feature, scales = "free") +
        ggplot2::labs(x = "Fitted value", y = "Residual", title = "Residual diagnostics") +
        ggplot2::theme_classic(base_size = 11)
    )
  }

  qq <- do.call(rbind, lapply(split(dat$residual, dat$feature), function(value) {
    points <- stats::qqnorm(value, plot.it = FALSE)
    data.frame(theoretical = points$x, observed = points$y)
  }))
  qq$feature <- rep(names(split(dat$residual, dat$feature)), vapply(
    split(dat$residual, dat$feature), length, integer(1)
  ))
  ggplot2::ggplot(qq, ggplot2::aes(x = .data$theoretical, y = .data$observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = 2) +
    ggplot2::geom_point(alpha = 0.55, size = 1) +
    ggplot2::facet_wrap(~feature, scales = "free") +
    ggplot2::labs(x = "Theoretical normal quantile", y = "Residual quantile", title = "Residual Q-Q plots") +
    ggplot2::theme_classic(base_size = 11)
}

#' Plot subject-specific fitted rhythms
#'
#' @param x A fitted mixedcirc object with subject random effects.
#' @param features Optional feature names to include.
#' @param npoints Number of points per reconstructed curve.
#'
#' @return A faceted ggplot object.
#' @export
mixedcirc_plot_subjects <- function(x, features = NULL, npoints = 100) {
  if (!is.numeric(npoints) || length(npoints) != 1L || npoints < 2 || npoints %% 1 != 0) {
    stop("npoints must be an integer of at least 2", call. = FALSE)
  }
  subject <- mixedcirc_subject_rhythms(x, features = features)
  rows <- lapply(seq_len(nrow(subject)), function(index) {
    time <- seq(0, subject$period[index], length.out = npoints)
    data.frame(
      feature = subject$feature[index],
      subject = subject$subject[index],
      group = subject$group[index],
      time = time,
      value = subject$subject_mesor[index] + subject$subject_amplitude[index] *
        cos(2 * pi * time / subject$period[index] + subject$subject_phase[index]),
      stringsAsFactors = FALSE
    )
  })
  dat <- do.call(rbind, rows)
  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$time, y = .data$value,
      group = interaction(.data$feature, .data$subject), color = .data$group
    )
  ) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.45) +
    ggplot2::facet_wrap(~feature, scales = "free_y") +
    ggplot2::labs(
      x = "Circadian time", y = "Subject-specific fitted value",
      color = "Group", title = "Individual rhythm heterogeneity"
    ) +
    ggplot2::theme_classic(base_size = 11)
}
