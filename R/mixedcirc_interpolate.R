#' Generate interpolated circadian rhythm predictions
#'
#' This function takes one or more fitted `mixedcirc` objects and generates
#' interpolated circadian predictions across a sequence of time points.
#'
#'
#' @param x An object of class `mixedcirc_fit` or `mixedcirc_fit_list`.
#' @param period Rhythm period. When omitted, the period stored by
#'   `mixedcirc_detect()` is used; the fallback is `24` for older objects.
#' @param min_time Minimum time used for interpolation. If `NULL`, taken from the
#'   fitted object's stored experimental design.
#' @param max_time Maximum time used for interpolation. If `NULL`, taken from the
#'   fitted object's stored experimental design.
#' @param npoints Number of interpolated time points to generate. Default `200`.
#' @param fixed_values Named list of covariate values to force in prediction.
#'   For example `list(age = 50, batch = "B1")`.
#' @param expand_vars Character vector of extra covariates to keep separate in
#'   the interpolated output.
#' @param numeric_method How numeric nuisance covariates are fixed:
#'   `"mean"` or `"median"`. Default `"mean"`.
#' @param factor_method How nuisance factor covariates are handled:
#'   `"average_observed"`, `"average_equal"`, or `"reference"`.
#'   Default `"average_observed"`.
#' @param add_expanded_label Logical. If `TRUE`, adds the `.expanded_label`
#'   column to the returned data. Default `TRUE`.
#'
#' @return
#' If `x` is a `mixedcirc_fit`, a `data.frame` containing the interpolated
#' prediction grid together with a `Y.hat` column of fitted values.
#'
#' If `x` is a `mixedcirc_fit_list`, a single combined `data.frame` is returned,
#' with one additional column named `feature` identifying the source feature.
#'
#' @details
#' This function uses the same interpolation and prediction rules as
#' `mixedcirc_fit_plot()`.
#'
#' @examples
#' data("circa_data")
#'
#' results <- mixedcirc_detect(
#'   data_input = circa_data$data_matrix,
#'   meta_input = data.frame(
#'     time = circa_data$time,
#'     group = circa_data$group,
#'     id = circa_data$id
#'   ),
#'   time = "time",
#'   group = "group",
#'   id = "id",
#'   period = 24,
#'   lm_method = "lm",
#'   verbose = TRUE
#' )
#'
#' int_data_one <- mixedcirc_interpolate(results[1])
#' int_data_all <- mixedcirc_interpolate(results)
#'
#' @export
mixedcirc_interpolate <- function(x,
                                  period = 24,
                                  min_time = NULL,
                                  max_time = NULL,
                                  npoints = 200,
                                  fixed_values = list(),
                                  expand_vars = NULL,
                                  numeric_method = c("mean", "median")[1],
                                  factor_method = c("average_observed", "average_equal", "reference")[1],
                                  add_expanded_label = TRUE) {

  if (!(is(x, "mixedcirc_fit") || is(x, "mixedcirc_fit_list"))) {
    stop("x must be an object of class `mixedcirc_fit` or `mixedcirc_fit_list`.")
  }
  if (is(x, "mixedcirc_fit_list") && length(x) == 0L) {
    stop("x must contain at least one fitted feature.")
  }

  first_object <- if (is(x, "mixedcirc_fit")) x else x[[1L]]
  if (missing(period)) {
    if (is(x, "mixedcirc_fit_list")) {
      stored_periods <- vapply(x@results, function(object) {
        value <- object@params$period
        if (is.null(value)) NA_real_ else as.numeric(value)[1L]
      }, numeric(1))
      known_periods <- unique(stored_periods[is.finite(stored_periods)])
      if (length(known_periods) > 1L) {
        stop("x contains fits with different stored periods; supply period explicitly.")
      }
      if (length(known_periods) == 1L) period <- known_periods
    } else if (!is.null(first_object@params$period)) {
      period <- first_object@params$period
    }
  }
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0) {
    stop("period must be a single positive finite number.")
  }
  if (!is.numeric(npoints) || length(npoints) != 1L || !is.finite(npoints) ||
      npoints < 2 || npoints %% 1 != 0) {
    stop("npoints must be an integer of at least 2.")
  }
  npoints <- as.integer(npoints)
  if (!is.null(min_time) &&
      (!is.numeric(min_time) || length(min_time) != 1L || !is.finite(min_time))) {
    stop("min_time must be NULL or one finite number.")
  }
  if (!is.null(max_time) &&
      (!is.numeric(max_time) || length(max_time) != 1L || !is.finite(max_time))) {
    stop("max_time must be NULL or one finite number.")
  }
  if (!is.null(min_time) && !is.null(max_time) && max_time <= min_time) {
    stop("max_time must be greater than min_time.")
  }
  if (!is.list(fixed_values) ||
      (length(fixed_values) > 0L && (is.null(names(fixed_values)) || any(names(fixed_values) == "")))) {
    stop("fixed_values must be a named list.")
  }
  if (!is.null(expand_vars) && !is.character(expand_vars)) {
    stop("expand_vars must be NULL or a character vector.")
  }
  if (!is.logical(add_expanded_label) || length(add_expanded_label) != 1L || is.na(add_expanded_label)) {
    stop("add_expanded_label must be TRUE or FALSE.")
  }

  numeric_method <- match.arg(numeric_method, c("mean", "median"))
  factor_method  <- match.arg(factor_method, c("average_observed", "average_equal", "reference"))

  .get_feature_name <- function(object, fallback = "feature") {
    nm <- tryCatch(rownames(object@results), error = function(e) NULL)
    if (is.null(nm) || length(nm) == 0 || is.na(nm[1]) || nm[1] == "") {
      return(fallback)
    }
    nm[1]
  }

  .interpolate_one <- function(object,
                               period,
                               min_time,
                               max_time,
                               npoints,
                               fixed_values,
                               expand_vars,
                               numeric_method,
                               factor_method,
                               add_expanded_label) {

    pred_df <- .build_prediction_data(
      object = object,
      period = period,
      min_time = min_time,
      max_time = max_time,
      n_points = npoints,
      fixed_values = fixed_values,
      expand_vars = expand_vars,
      numeric_method = numeric_method,
      factor_method = factor_method
    )

    pred_df <- .predict_circadian_curve(
      object = object,
      pred_df = pred_df
    )

    if (isTRUE(add_expanded_label)) {
      pred_df <- .add_expanded_label(pred_df, expand_vars = expand_vars)
    }

    pred_df <- as.data.frame(pred_df)
    rownames(pred_df) <- NULL
    pred_df
  }

  if (is(x, "mixedcirc_fit")) {
    return(
      .interpolate_one(
        object = x,
        period = period,
        min_time = min_time,
        max_time = max_time,
        npoints = npoints,
        fixed_values = fixed_values,
        expand_vars = expand_vars,
        numeric_method = numeric_method,
        factor_method = factor_method,
        add_expanded_label = add_expanded_label
      )
    )
  }

  # mixedcirc_fit_list -> return one combined data.frame
  out_list <- vector("list", length(x))

  for (i in seq_along(x)) {
    object <- x[i]
    feature_name <- .get_feature_name(object, fallback = paste0("feature_", i))

    df_i <- .interpolate_one(
      object = object,
      period = period,
      min_time = min_time,
      max_time = max_time,
      npoints = npoints,
      fixed_values = fixed_values,
      expand_vars = expand_vars,
      numeric_method = numeric_method,
      factor_method = factor_method,
      add_expanded_label = add_expanded_label
    )

    df_i$feature <- feature_name
    out_list[[i]] <- df_i
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
