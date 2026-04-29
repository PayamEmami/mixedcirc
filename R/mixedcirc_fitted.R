#' Extract fitted circadian values
#'
#' This function extracts fitted values from one or more fitted `mixedcirc`
#' models and returns them as a matrix aligned with the stored experimental
#' design.
#'
#' @param fit An object of class `mixedcirc_fit` or `mixedcirc_fit_list`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default `FALSE`.
#'
#' @return
#' A numeric matrix of fitted values. Rows correspond to the rows of
#' `exp_design` stored in the fitted object, and columns correspond to features.
#'
#' @details
#' This function returns the fitted values stored in the `.fitted` column of
#' `exp_design`, which is created by the updated `mixedcirc_detect()`.
#' If some rows were dropped during model fitting, their fitted values remain `NA`.
#'
#' @examples
#' library(mixedcirc)
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
#' fitted_vals <- mixedcirc_fitted(results)
#'
#' @import methods
#' @export
mixedcirc_fitted <- function(fit = NULL, verbose = FALSE) {

  if (verbose) cat("Checking inputs ...\n")

  if (!(is(fit, "mixedcirc_fit") || is(fit, "mixedcirc_fit_list"))) {
    stop("input must be `mixedcirc_fit` or `mixedcirc_fit_list`")
  }

  if (verbose) cat("Extracting fitted values ...\n")

  .get_feature_name <- function(object, fallback = "feature") {
    nm <- tryCatch(rownames(object@results), error = function(e) NULL)
    if (is.null(nm) || length(nm) == 0 || is.na(nm[1]) || nm[1] == "") {
      return(fallback)
    }
    nm[1]
  }

  .extract_single_fitted <- function(object, verbose = FALSE) {
    exp_design <- object@exp_design
    feature_name <- .get_feature_name(object)

    if (verbose) {
      cat("Extracting fitted values for", feature_name, "...\n")
    }

    if (!".fitted" %in% names(exp_design)) {
      stop("exp_design does not contain `.fitted`. Please refit with the updated mixedcirc_detect().")
    }

    out_mat <- matrix(as.numeric(exp_design$.fitted), ncol = 1)
    colnames(out_mat) <- feature_name
    out_mat
  }

  if (is(fit, "mixedcirc_fit")) {
    return(.extract_single_fitted(fit, verbose = verbose))
  }

  out_list <- vector("list", length(fit))

  for (i in seq_along(fit)) {
    out_list[[i]] <- .extract_single_fitted(fit[i], verbose = verbose)
  }

  out <- do.call(cbind, out_list)
  out
}
