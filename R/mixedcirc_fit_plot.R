# =========================================================
# Helper functions for mixedcirc plotting
# =========================================================

.get_model_frame <- function(fit) {
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)

  if (!is.null(mf)) {
    return(mf)
  }

  if (inherits(fit, "merMod") && !is.null(fit@frame)) {
    return(fit@frame)
  }

  if (inherits(fit, "lm") && !is.null(fit$model)) {
    return(fit$model)
  }

  stop("Could not extract model frame from `fit`.")
}

.get_fit_classes <- function(fit) {
  list(
    is_lm  = inherits(fit, "lm") && !inherits(fit, c("merMod", "lme")),
    is_mer = inherits(fit, "merMod"),
    is_lme = inherits(fit, "lme")
  )
}

.get_prediction_formula_vars <- function(fit) {
  if (inherits(fit, "merMod")) {
    ff <- lme4::nobars(stats::formula(fit))
  } else {
    ff <- stats::formula(fit)
  }

  tt <- stats::terms(ff)
  unique(all.vars(stats::delete.response(tt)))
}

.make_same_class_as_model <- function(newdata, mf) {
  common_vars <- intersect(names(newdata), names(mf))

  for (v in common_vars) {
    if (is.factor(mf[[v]])) {
      newdata[[v]] <- factor(newdata[[v]], levels = levels(mf[[v]]))
    } else if (is.logical(mf[[v]])) {
      newdata[[v]] <- as.logical(newdata[[v]])
    } else if (is.character(mf[[v]])) {
      newdata[[v]] <- as.character(newdata[[v]])
    } else if (is.numeric(mf[[v]])) {
      newdata[[v]] <- as.numeric(newdata[[v]])
    }
  }

  newdata
}

.is_factor_like <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

.first_non_missing <- function(x) {
  x2 <- x[!is.na(x)]
  if (length(x2) == 0) return(NA)
  x2[1]
}

.get_call_arg_name <- function(object, arg_name) {
  if (!"call" %in% methods::slotNames(object)) {
    return(NULL)
  }

  cl <- object@call
  if (is.null(cl) || is.null(cl[[arg_name]])) {
    return(NULL)
  }

  out <- tryCatch(as.character(cl[[arg_name]]), error = function(e) NULL)
  if (length(out) == 0) return(NULL)
  out[1]
}

.get_plot_spec <- function(object) {
  exp_design <- object@exp_design
  params <- if ("params" %in% methods::slotNames(object)) object@params else list()

  time_var <- .get_call_arg_name(object, "time")
  group_var_requested <- .get_call_arg_name(object, "group")
  id_var <- .get_call_arg_name(object, "id")
  replicate_var <- .get_call_arg_name(object, "replicate_id")

  multiple_groups <- isTRUE(params$multiple_groups)
  rrbs <- isTRUE(params$RRBS)

  if (is.null(time_var) || !time_var %in% names(exp_design)) {
    if ("time" %in% names(exp_design)) {
      time_var <- "time"
    } else {
      stop("Could not identify the time column in exp_design.")
    }
  }

  group_var <- NULL
  if (multiple_groups) {
    if (!is.null(group_var_requested) && group_var_requested %in% names(exp_design)) {
      group_var <- group_var_requested
    } else if ("group" %in% names(exp_design)) {
      group_var <- "group"
    } else {
      stop("Could not identify the group column in exp_design.")
    }
  }

  if (is.null(id_var) || !id_var %in% names(exp_design)) {
    if ("id" %in% names(exp_design)) {
      id_var <- "id"
    } else if ("rep" %in% names(exp_design)) {
      id_var <- "rep"
    } else {
      id_var <- NULL
    }
  }

  if (is.null(replicate_var) || !replicate_var %in% names(exp_design)) {
    if ("replicate_id" %in% names(exp_design)) {
      replicate_var <- "replicate_id"
    } else {
      replicate_var <- NULL
    }
  }

  list(
    time_var = time_var,
    group_var = group_var,
    id_var = id_var,
    replicate_var = replicate_var,
    multiple_groups = multiple_groups,
    RRBS = rrbs
  )
}

.get_factor_weights <- function(x, method = c("average_observed", "average_equal", "reference")) {
  method <- match.arg(method)

  if (is.factor(x)) {
    levs <- levels(x)
    x_chr <- as.character(x)
  } else if (is.logical(x)) {
    levs <- sort(unique(x))
    x_chr <- as.character(x)
  } else {
    levs <- sort(unique(as.character(x)))
    x_chr <- as.character(x)
  }

  levs <- levs[!is.na(levs)]

  if (length(levs) == 0) {
    return(data.frame(level = NA, weight = 1))
  }

  if (method == "reference") {
    return(data.frame(level = levs[1], weight = 1, stringsAsFactors = FALSE))
  }

  if (method == "average_equal") {
    w <- rep(1 / length(levs), length(levs))
    return(data.frame(level = levs, weight = w, stringsAsFactors = FALSE))
  }

  tab <- table(factor(x_chr, levels = levs))
  w <- as.numeric(tab) / sum(tab)
  data.frame(level = levs, weight = w, stringsAsFactors = FALSE)
}

.build_covariate_grid <- function(newdata,
                                  fit,
                                  fixed_values = list(),
                                  expand_vars = NULL,
                                  numeric_method = c("mean", "median"),
                                  factor_method = c("average_observed", "average_equal", "reference"),
                                  exclude_vars = c("measure", ".fitted", ".residual", ".weight",
                                                   ".used_in_fit", "inphase", "outphase",
                                                   "id", "replicate_id")) {

  numeric_method <- match.arg(numeric_method)
  factor_method  <- match.arg(factor_method)

  mf <- .get_model_frame(fit)
  model_vars <- .get_prediction_formula_vars(fit)

  vars_to_fill <- setdiff(model_vars, c(names(newdata), exclude_vars))
  if (length(vars_to_fill) == 0) {
    newdata$.avg_weight <- 1
    return(.make_same_class_as_model(newdata, mf))
  }

  fill_list <- list()
  hidden_avg_vars <- character(0)

  for (v in vars_to_fill) {
    if (!v %in% names(mf)) next
    x <- mf[[v]]

    if (v %in% names(fixed_values)) {
      fill_list[[v]] <- fixed_values[[v]]
      next
    }

    if (!is.null(expand_vars) && v %in% expand_vars) {
      if (is.factor(x)) {
        fill_list[[v]] <- levels(x)
      } else if (is.logical(x)) {
        fill_list[[v]] <- sort(unique(x))
      } else if (is.character(x)) {
        fill_list[[v]] <- sort(unique(x))
      } else if (is.numeric(x)) {
        fill_list[[v]] <- sort(unique(x))
      } else {
        fill_list[[v]] <- unique(x)
      }
      next
    }

    if (is.numeric(x)) {
      fill_list[[v]] <- if (numeric_method == "mean") {
        mean(x, na.rm = TRUE)
      } else {
        stats::median(x, na.rm = TRUE)
      }
    } else if (.is_factor_like(x)) {
      wdf <- .get_factor_weights(x, method = factor_method)
      fill_list[[v]] <- wdf$level
      fill_list[[paste0(".w__", v)]] <- wdf$weight
      if (factor_method != "reference") {
        hidden_avg_vars <- c(hidden_avg_vars, v)
      }
    } else {
      fill_list[[v]] <- .first_non_missing(x)
    }
  }

  if (length(fill_list) > 0) {
    extra_grid <- do.call(expand.grid, c(fill_list, stringsAsFactors = FALSE))
    newdata <- merge(newdata, extra_grid, by = NULL, sort = FALSE)
  }

  weight_cols <- grep("^\\.w__", names(newdata), value = TRUE)
  if (length(weight_cols) == 0) {
    newdata$.avg_weight <- 1
  } else {
    newdata$.avg_weight <- apply(newdata[, weight_cols, drop = FALSE], 1, prod)
    newdata[weight_cols] <- NULL
  }

  attr(newdata, "hidden_avg_vars") <- unique(hidden_avg_vars)
  .make_same_class_as_model(newdata, mf)
}

.add_expanded_label <- function(dat, expand_vars = NULL) {
  if (is.null(expand_vars) || length(expand_vars) == 0) {
    dat$.expanded_label <- "default"
    return(dat)
  }

  expand_vars <- intersect(expand_vars, names(dat))

  if (length(expand_vars) == 0) {
    dat$.expanded_label <- "default"
    return(dat)
  }

  dat$.expanded_label <- apply(
    dat[, expand_vars, drop = FALSE],
    1,
    function(z) paste(paste0(expand_vars, "=", z), collapse = ", ")
  )

  dat
}

.collapse_rrbs_difference <- function(dat,
                                      value_col = "value",
                                      scaler_col = "scaler") {

  if (!scaler_col %in% names(dat)) {
    stop("`", scaler_col, "` not found in data.")
  }
  if (!value_col %in% names(dat)) {
    stop("`", value_col, "` not found in data.")
  }

  dat1 <- dat[dat[[scaler_col]] == 1, , drop = FALSE]
  dat0 <- dat[dat[[scaler_col]] == 0, , drop = FALSE]

  if (nrow(dat1) == 0 || nrow(dat0) == 0) {
    stop("For RRBS plotting, both scaler == 1 and scaler == 0 must exist.")
  }

  key_cols <- setdiff(names(dat), c(value_col, scaler_col))

  dat1$.pair_id <- ave(seq_len(nrow(dat1)),
                       interaction(dat1[, key_cols, drop = FALSE], drop = TRUE),
                       FUN = seq_along)
  dat0$.pair_id <- ave(seq_len(nrow(dat0)),
                       interaction(dat0[, key_cols, drop = FALSE], drop = TRUE),
                       FUN = seq_along)

  merge_keys <- c(key_cols, ".pair_id")

  merged <- merge(
    dat1[, c(merge_keys, value_col), drop = FALSE],
    dat0[, c(merge_keys, value_col), drop = FALSE],
    by = merge_keys,
    suffixes = c(".1", ".0"),
    sort = FALSE
  )

  out <- merged[, key_cols, drop = FALSE]
  out[[value_col]] <- merged[[paste0(value_col, ".1")]] - merged[[paste0(value_col, ".0")]]

  out
}

.average_hidden_predictions <- function(pred_df, hidden_avg_vars = NULL, value_col = "Y.hat") {
  if (is.null(hidden_avg_vars) || length(hidden_avg_vars) == 0) {
    return(pred_df)
  }

  hidden_avg_vars <- intersect(hidden_avg_vars, names(pred_df))
  if (length(hidden_avg_vars) == 0) {
    return(pred_df)
  }

  if (!".avg_weight" %in% names(pred_df)) {
    pred_df$.avg_weight <- 1
  }

  keep_cols <- setdiff(names(pred_df), c(hidden_avg_vars, value_col, ".avg_weight"))
  if (length(keep_cols) == 0) {
    stop("No grouping columns left after averaging hidden nuisance variables.")
  }

  split_key <- interaction(pred_df[, keep_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  idx_list <- split(seq_len(nrow(pred_df)), split_key)

  out_list <- lapply(idx_list, function(ii) {
    tmp <- pred_df[ii, , drop = FALSE]
    w <- tmp$.avg_weight
    if (all(is.na(w)) || sum(w, na.rm = TRUE) == 0) {
      w <- rep(1, nrow(tmp))
    }
    val <- stats::weighted.mean(tmp[[value_col]], w = w, na.rm = TRUE)
    row <- tmp[1, keep_cols, drop = FALSE]
    row[[value_col]] <- val
    row
  })

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

.get_used_exp_design <- function(object) {
  exp_design <- object@exp_design
  if (!".used_in_fit" %in% names(exp_design)) {
    stop("object@exp_design must contain `.used_in_fit`. Please refit with the updated mixedcirc_detect().")
  }
  exp_design[exp_design$.used_in_fit %in% TRUE, , drop = FALSE]
}

.build_prediction_data <- function(object,
                                   period = 24,
                                   min_time = NULL,
                                   max_time = NULL,
                                   n_points = 200,
                                   fixed_values = list(),
                                   expand_vars = NULL,
                                   numeric_method = c("mean", "median"),
                                   factor_method = c("average_observed", "average_equal", "reference")) {

  numeric_method <- match.arg(numeric_method)
  factor_method  <- match.arg(factor_method)

  fit <- object@fit
  mf <- .get_model_frame(fit)
  exp_design <- .get_used_exp_design(object)
  spec <- .get_plot_spec(object)

  if (nrow(exp_design) == 0) {
    stop("No rows marked as used in fit inside exp_design.")
  }

  time_var <- spec$time_var
  group_var <- spec$group_var

  if (is.null(min_time)) min_time <- min(exp_design[[time_var]], na.rm = TRUE)
  if (is.null(max_time)) max_time <- max(exp_design[[time_var]], na.rm = TRUE)

  timeax <- seq(min_time, max_time, length.out = n_points)

  has_group  <- !is.null(group_var) && group_var %in% names(exp_design) && group_var %in% names(mf)
  has_scaler <- "scaler" %in% names(mf)
  has_rep_id <- !is.null(spec$replicate_var) && spec$replicate_var %in% names(mf)

  group_levels <- if (has_group) unique(as.character(exp_design[[group_var]])) else NA_character_

  pred_list <- vector("list", 0)
  hidden_avg_all <- character(0)

  for (gr in group_levels) {
    base_newdata <- data.frame(
      inphase = cos(2 * pi * timeax / period),
      outphase = sin(2 * pi * timeax / period)
    )
    base_newdata[[time_var]] <- timeax

    if (has_group) {
      base_newdata[[group_var]] <- gr
    }

    if (spec$RRBS && has_scaler) {
      rep_id <- if (has_rep_id) unique(as.character(mf[[spec$replicate_var]]))[1] else NULL

      nd1 <- base_newdata
      nd1$scaler <- 1
      if (!is.null(rep_id)) nd1[[spec$replicate_var]] <- rep_id

      nd0 <- base_newdata
      nd0$scaler <- 0
      if (!is.null(rep_id)) nd0[[spec$replicate_var]] <- rep_id

      nd1 <- .build_covariate_grid(
        newdata = nd1,
        fit = fit,
        fixed_values = fixed_values,
        expand_vars = expand_vars,
        numeric_method = numeric_method,
        factor_method = factor_method,
        exclude_vars = c("measure", ".fitted", ".residual", ".weight", ".used_in_fit",
                         time_var, "inphase", "outphase", spec$id_var, spec$replicate_var)
      )
      nd0 <- .build_covariate_grid(
        newdata = nd0,
        fit = fit,
        fixed_values = fixed_values,
        expand_vars = expand_vars,
        numeric_method = numeric_method,
        factor_method = factor_method,
        exclude_vars = c("measure", ".fitted", ".residual", ".weight", ".used_in_fit",
                         time_var, "inphase", "outphase", spec$id_var, spec$replicate_var)
      )

      hidden_avg_all <- union(hidden_avg_all, attr(nd1, "hidden_avg_vars"))
      hidden_avg_all <- union(hidden_avg_all, attr(nd0, "hidden_avg_vars"))

      pred_list[[length(pred_list) + 1]] <- nd1
      pred_list[[length(pred_list) + 1]] <- nd0
    } else {
      nd <- base_newdata

      if (has_scaler && "scaler" %in% .get_prediction_formula_vars(fit)) {
        nd$scaler <- if ("scaler" %in% names(fixed_values)) fixed_values[["scaler"]] else 1
      }

      nd <- .build_covariate_grid(
        newdata = nd,
        fit = fit,
        fixed_values = fixed_values,
        expand_vars = expand_vars,
        numeric_method = numeric_method,
        factor_method = factor_method,
        exclude_vars = c("measure", ".fitted", ".residual", ".weight", ".used_in_fit",
                         time_var, "inphase", "outphase", spec$id_var, spec$replicate_var)
      )

      hidden_avg_all <- union(hidden_avg_all, attr(nd, "hidden_avg_vars"))
      pred_list[[length(pred_list) + 1]] <- nd
    }
  }

  pred_df <- do.call(rbind, pred_list)
  rownames(pred_df) <- NULL
  attr(pred_df, "hidden_avg_vars") <- hidden_avg_all
  pred_df
}

.predict_circadian_curve <- function(object, pred_df) {
  fit <- object@fit
  cls <- .get_fit_classes(fit)
  hidden_avg_vars <- attr(pred_df, "hidden_avg_vars")

  if (cls$is_mer) {
    pred_df$Y.hat <- stats::predict(
      object = fit,
      newdata = pred_df,
      re.form = NA,
      allow.new.levels = TRUE
    )
  } else if (cls$is_lme) {
    pred_df$Y.hat <- stats::predict(
      object = fit,
      newdata = pred_df,
      level = 0
    )
  } else {
    pred_df$Y.hat <- stats::predict(
      object = fit,
      newdata = pred_df
    )
  }

  if (.get_plot_spec(object)$RRBS && "scaler" %in% names(pred_df)) {
    pred_df <- .collapse_rrbs_difference(pred_df, value_col = "Y.hat", scaler_col = "scaler")
  }

  pred_df <- .average_hidden_predictions(pred_df, hidden_avg_vars = hidden_avg_vars, value_col = "Y.hat")
  pred_df
}

.get_raw_plot_data <- function(object) {
  exp_design <- .get_used_exp_design(object)
  spec <- .get_plot_spec(object)

  if (!"measure" %in% names(exp_design)) {
    stop("exp_design must contain a `measure` column.")
  }

  raw_data <- exp_design

  if (spec$RRBS && "scaler" %in% names(raw_data)) {
    raw_data <- .collapse_rrbs_difference(raw_data, value_col = "measure", scaler_col = "scaler")
  }

  raw_data
}

# =========================================================
# Main plotting function
# =========================================================

#' Plot circadian rhythm fits from mixedcirc
#'
#' Plots fitted circadian curves from a `mixedcirc_fit` object and can
#' automatically fix or average over extra covariates used in the model.
#'
#' @param x A `mixedcirc_fit` object.
#' @param y Not used.
#' @param period Rhythm period. Default 24.
#' @param min_time Minimum time for prediction. If NULL, taken from data.
#' @param max_time Maximum time for prediction. If NULL, taken from data.
#' @param n_points Number of time points used for the fitted curve.
#' @param plot_title Plot title. If NULL, derived from result rownames if possible.
#' @param plot_points Whether to plot raw points. Default TRUE.
#' @param plot_smooth Whether to add dashed smoothing based on raw data. Default FALSE.
#' @param plot_trend Whether to plot fitted curves. Default TRUE.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param fixed_values Named list of covariate values to force in prediction.
#' @param expand_vars Character vector of extra covariates to show separately.
#' @param show_expanded One of "none", "linetype", "facet", "label".
#' @param numeric_method How numeric nuisance covariates are fixed: "mean" or "median".
#' @param factor_method How nuisance factor covariates are handled:
#'   "average_observed", "average_equal", or "reference".
#'
#' @return A ggplot object.
#' @import stats
#' @import methods
#' @import lme4
#' @import ggplot2
#' @import ggpubr
#' @import ggsci
#' @export
mixedcirc_fit_plot <- function(x,
                               y = NULL,
                               period = 24,
                               min_time = NULL,
                               max_time = NULL,
                               n_points = 200,
                               plot_title = NULL,
                               plot_points = TRUE,
                               plot_smooth = FALSE,
                               plot_trend = TRUE,
                               xlab = NULL,
                               ylab = "Expression",
                               fixed_values = list(),
                               expand_vars = NULL,
                               show_expanded = c("none", "linetype", "facet", "label")[1],
                               numeric_method = c("mean", "median")[1],
                               factor_method = c("average_observed", "average_equal", "reference")[1]) {

  show_expanded  <- match.arg(show_expanded, c("none", "linetype", "facet", "label"))
  numeric_method <- match.arg(numeric_method, c("mean", "median"))
  factor_method  <- match.arg(factor_method, c("average_observed", "average_equal", "reference"))

  object <- x
  spec <- .get_plot_spec(object)
  time_var <- spec$time_var
  group_var <- spec$group_var

  if (is.null(xlab)) {
    xlab <- time_var
  }

  if (is.null(plot_title)) {
    plot_title <- tryCatch(rownames(object@results), error = function(e) NULL)
    if (length(plot_title) > 1) plot_title <- plot_title[1]
    if (is.null(plot_title) || length(plot_title) == 0 || is.na(plot_title)) {
      plot_title <- "Variable"
    }
  }

  pred_df <- .build_prediction_data(
    object = object,
    period = period,
    min_time = min_time,
    max_time = max_time,
    n_points = n_points,
    fixed_values = fixed_values,
    expand_vars = expand_vars,
    numeric_method = numeric_method,
    factor_method = factor_method
  )

  pred_df <- .predict_circadian_curve(object, pred_df)
  raw_data <- .get_raw_plot_data(object)

  pred_df <- .add_expanded_label(pred_df, expand_vars = expand_vars)
  raw_data <- .add_expanded_label(raw_data, expand_vars = expand_vars)

  has_group_pred <- !is.null(group_var) && group_var %in% names(pred_df)
  has_group_raw  <- !is.null(group_var) && group_var %in% names(raw_data)
  has_group <- has_group_pred && has_group_raw

  if (has_group) {
    plot_tmp <- ggplot2::ggplot(
      data = raw_data,
      ggplot2::aes_string(x = time_var, y = "measure", group = group_var, color = group_var, shape = group_var)
    )
  } else {
    plot_tmp <- ggplot2::ggplot(
      data = raw_data,
      ggplot2::aes_string(x = time_var, y = "measure")
    )
  }

  if (plot_points) {
    if (has_group) {
      plot_tmp <- plot_tmp +
        ggplot2::geom_point(
          data = raw_data,
          ggplot2::aes_string(x = time_var, y = "measure", group = group_var, color = group_var, shape = group_var),
          position = ggplot2::position_dodge(2),
          alpha = 0.2
        )
    } else {
      plot_tmp <- plot_tmp +
        ggplot2::geom_point(
          data = raw_data,
          ggplot2::aes_string(x = time_var, y = "measure"),
          alpha = 0.2
        )
    }
  }

  if (plot_smooth) {
    if (has_group) {
      if (!is.null(expand_vars) && length(expand_vars) > 0 && show_expanded == "linetype") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_smooth(
            data = raw_data,
            ggplot2::aes_string(
              x = time_var,
              y = "measure",
              group = paste0("interaction(", group_var, ", .expanded_label)"),
              color = group_var,
              linetype = ".expanded_label"
            ),
            se = FALSE,
            alpha = 0.2
          )
      } else {
        plot_tmp <- plot_tmp +
          ggplot2::geom_smooth(
            data = raw_data,
            ggplot2::aes_string(x = time_var, y = "measure", group = group_var, color = group_var),
            linetype = "dashed",
            se = FALSE,
            alpha = 0.2
          )
      }
    } else {
      if (!is.null(expand_vars) && length(expand_vars) > 0 && show_expanded == "linetype") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_smooth(
            data = raw_data,
            ggplot2::aes_string(x = time_var, y = "measure", group = ".expanded_label", linetype = ".expanded_label"),
            se = FALSE,
            alpha = 0.2
          )
      } else {
        plot_tmp <- plot_tmp +
          ggplot2::geom_smooth(
            data = raw_data,
            ggplot2::aes_string(x = time_var, y = "measure"),
            linetype = "dashed",
            se = FALSE,
            alpha = 0.2
          )
      }
    }
  }

  if (plot_trend) {
    if (has_group) {
      if (is.null(expand_vars) || length(expand_vars) == 0 || show_expanded == "none") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat", color = group_var)
          )
      } else if (show_expanded == "linetype") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(
              x = time_var,
              y = "Y.hat",
              color = group_var,
              linetype = ".expanded_label",
              group = paste0("interaction(", group_var, ", .expanded_label)")
            )
          )
      } else if (show_expanded == "label") {
        pred_df$.curve_label <- paste0(group_var, "=", pred_df[[group_var]], ", ", pred_df$.expanded_label)
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat", color = ".curve_label", group = ".curve_label")
          )
      } else if (show_expanded == "facet") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat", color = group_var)
          )
      }
    } else {
      if (is.null(expand_vars) || length(expand_vars) == 0 || show_expanded == "none") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat")
          )
      } else if (show_expanded == "linetype") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat", linetype = ".expanded_label", group = ".expanded_label")
          )
      } else if (show_expanded == "label") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat", color = ".expanded_label", group = ".expanded_label")
          )
      } else if (show_expanded == "facet") {
        plot_tmp <- plot_tmp +
          ggplot2::geom_line(
            data = pred_df,
            ggplot2::aes_string(x = time_var, y = "Y.hat")
          )
      }
    }
  }

  if (!is.null(expand_vars) && length(expand_vars) > 0 && show_expanded == "facet") {
    plot_tmp <- plot_tmp + ggplot2::facet_wrap(~ .expanded_label)
  }

  plot_tmp <- plot_tmp +
    ggpubr::theme_classic2() +
    ggsci::scale_color_aaas() +
    ggplot2::scale_x_continuous(
      labels = sort(unique(raw_data[[time_var]])),
      breaks = sort(unique(raw_data[[time_var]]))
    ) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_y_continuous() +
    ggplot2::theme(axis.text = ggplot2::element_text(color = "black")) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlab)

  if (has_group) {
    if (length(unique(raw_data[[group_var]])) < 2 &&
        (is.null(expand_vars) || length(expand_vars) == 0 || show_expanded == "none")) {
      plot_tmp <- plot_tmp + ggplot2::theme(legend.position = "none")
    }
  } else {
    if (is.null(expand_vars) || length(expand_vars) == 0 || show_expanded == "none") {
      plot_tmp <- plot_tmp + ggplot2::theme(legend.position = "none")
    }
  }

  plot_tmp
}
