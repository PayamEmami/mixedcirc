#' Generate interpolated circadian rhythm predictions
#'
#' This function takes one or more fitted `mixedcirc` objects and generates
#' interpolated circadian predictions across a sequence of time points.
#'
#'
#' @param x An object of class `mixedcirc_fit` or `mixedcirc_fit_list`.
#' @param period Rhythm period. Default `24`.
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
#'   verbose = TRUE
#' )
#'
#' int_data_one <- mixedcirc_interpolate(results[1])
#' int_data_all <- mixedcirc_interpolate(results)
#'
#' @import methods
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

  if (!is.numeric(npoints) || length(npoints) != 1 || is.na(npoints) || npoints <= 0) {
    stop("npoints must be a single positive numeric value.")
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
