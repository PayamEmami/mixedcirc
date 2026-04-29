#' Remove linear trend from fitted circadian values
#'
#' This function removes a linear trend from the fitted values obtained from
#' a `mixedcirc_fit` or `mixedcirc_fit_list` object. The detrending is performed
#' by regressing the fitted values on the original time variable and returning
#' the residuals from that regression.
#'
#' @param fit An object of class `mixedcirc_fit` or `mixedcirc_fit_list`.
#' @param per_group Logical. If `TRUE`, detrending is performed separately within
#'   each group when group information was used in the original model. Default `FALSE`.
#' @param use_mixed_for_detrend Logical. If `TRUE`, mixed-model fits are detrended
#'   using a random intercept for the original subject identifier when possible.
#'   If the mixed detrending model fails, the function falls back to `lm`.
#'   Default `TRUE`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default `FALSE`.
#' @param ... Additional arguments passed to the detrending model (`lm` or `lmer`).
#'
#' @return
#' A numeric matrix of detrended fitted values. Rows correspond to rows in the
#' original `exp_design` stored inside the fitted object, and columns correspond
#' to features.
#'
#' @details
#' The function first tries to use the `.fitted` column stored in `object@exp_design`
#' by the updated `mixedcirc_detect()` function. If that column is unavailable or
#' unusable, fitted values are extracted directly from the fitted model and aligned
#' back to `exp_design`.
#'
#' It then fits a regression of the form
#' \preformatted{
#' fitted_value ~ time
#' }
#' and returns the residuals from that regression.
#'
#' If `per_group = TRUE`, the detrending is performed separately within each group.
#'
#' For mixed-model fits, a detrending model with a random intercept,
#' \preformatted{
#' fitted_value ~ time + (1 | id)
#' }
#' is attempted when the original `id` variable is available and
#' `use_mixed_for_detrend = TRUE`. If this fails, the function falls back to
#' ordinary weighted linear regression.
#'
#' Stored `.weight` values are only used when they exist and contain at least
#' some non-missing values. If `.weight` is entirely `NA` (for example because
#' the original model was unweighted), weights are ignored automatically.
#'
#' @import stats
#' @import methods
#' @import lme4
#' @export
mixedcirc_detrend <- function(fit = NULL,
                              per_group = FALSE,
                              use_mixed_for_detrend = TRUE,
                              verbose = FALSE,
                              ...) {

  if (verbose) cat("Checking inputs ...\n")

  if (!(is(fit, "mixedcirc_fit") || is(fit, "mixedcirc_fit_list"))) {
    stop("input must be `mixedcirc_fit` or `mixedcirc_fit_list`")
  }

  # ---------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------

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

  .get_spec <- function(object) {
    exp_design <- object@exp_design
    params <- if ("params" %in% methods::slotNames(object)) object@params else list()

    time_var <- .get_call_arg_name(object, "time")
    group_var_requested <- .get_call_arg_name(object, "group")
    id_var <- .get_call_arg_name(object, "id")

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
        warning("Grouped analysis detected, but group column could not be identified. Detrending without grouping.")
        group_var <- NULL
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

    list(
      time_var = time_var,
      group_var = group_var,
      id_var = id_var,
      multiple_groups = multiple_groups,
      RRBS = rrbs
    )
  }

  .get_feature_name <- function(object, fallback = "feature") {
    nm <- tryCatch(rownames(object@results), error = function(e) NULL)
    if (is.null(nm) || length(nm) == 0 || is.na(nm[1]) || nm[1] == "") {
      return(fallback)
    }
    nm[1]
  }

  .get_stored_weights <- function(object) {
    exp_design <- object@exp_design

    if (!".weight" %in% names(exp_design)) {
      return(NULL)
    }

    w <- suppressWarnings(as.numeric(exp_design$.weight))

    if (length(w) != nrow(exp_design)) {
      return(NULL)
    }

    if (all(is.na(w))) {
      return(NULL)
    }

    w
  }

  .get_model_frame <- function(model_obj) {
    mf <- tryCatch(stats::model.frame(model_obj), error = function(e) NULL)

    if (!is.null(mf)) {
      return(mf)
    }

    if (inherits(model_obj, "merMod") && !is.null(model_obj@frame)) {
      return(model_obj@frame)
    }

    if (inherits(model_obj, "lm") && !is.null(model_obj$model)) {
      return(model_obj$model)
    }

    NULL
  }

  .extract_aligned_fitted <- function(object) {
    exp_design <- object@exp_design
    fit_obj <- object@fit

    if (".fitted" %in% names(exp_design)) {
      fitted_stored <- suppressWarnings(as.numeric(exp_design$.fitted))
      if (length(fitted_stored) == nrow(exp_design) && any(!is.na(fitted_stored))) {
        return(fitted_stored)
      }
    }

    fitted_raw <- tryCatch(as.numeric(stats::fitted(fit_obj)), error = function(e) NULL)
    if (is.null(fitted_raw)) {
      stop("Could not extract fitted values from the model.")
    }

    out <- rep(NA_real_, nrow(exp_design))
    mf <- .get_model_frame(fit_obj)

    if (is.null(mf)) {
      if (length(fitted_raw) == nrow(exp_design)) {
        return(fitted_raw)
      }
      stop("Could not align fitted values back to exp_design because the model frame could not be extracted.")
    }

    mf_rows <- rownames(mf)
    exp_rows <- rownames(exp_design)

    if (!is.null(mf_rows) && !is.null(exp_rows)) {
      pos <- match(mf_rows, exp_rows)
      if (all(!is.na(pos))) {
        out[pos] <- fitted_raw
        return(out)
      }
    }

    if (".used_in_fit" %in% names(exp_design)) {
      used_idx <- which(!is.na(exp_design$.used_in_fit) & exp_design$.used_in_fit)
      if (length(used_idx) == length(fitted_raw)) {
        out[used_idx] <- fitted_raw
        return(out)
      }
    }

    if (length(fitted_raw) == nrow(exp_design)) {
      return(fitted_raw)
    }

    stop("Could not align fitted values back to exp_design.")
  }

  .safe_detrend_subset <- function(fitted_v,
                                   time_v,
                                   id_v = NULL,
                                   weights_v = NULL,
                                   use_mixed = TRUE,
                                   ...) {

    keep <- !is.na(fitted_v) & !is.na(time_v)

    if (!is.null(weights_v)) {
      if (all(is.na(weights_v))) {
        weights_v <- NULL
      } else {
        keep_w <- !is.na(weights_v)
        if (sum(keep & keep_w) >= 2) {
          keep <- keep & keep_w
        } else {
          weights_v <- NULL
        }
      }
    }

    out <- rep(NA_real_, length(fitted_v))

    if (sum(keep) < 2) {
      out[keep] <- fitted_v[keep]
      return(out)
    }

    fitted_sub <- fitted_v[keep]
    time_sub <- time_v[keep]
    weights_sub <- if (!is.null(weights_v)) weights_v[keep] else NULL

    detrended_sub <- NULL

    if (use_mixed && !is.null(id_v)) {
      id_keep <- id_v[keep]
      mixed_keep <- !is.na(id_keep)

      if (sum(mixed_keep) >= 2 && length(unique(id_keep[mixed_keep])) > 1) {
        dat_sub <- data.frame(
          measure = fitted_sub[mixed_keep],
          time = time_sub[mixed_keep],
          id_tmp = as.factor(id_keep[mixed_keep])
        )

        weights_mixed <- if (!is.null(weights_sub)) weights_sub[mixed_keep] else NULL

        fit_try <- tryCatch(
          {
            if (is.null(weights_mixed)) {
              lme4::lmer(measure ~ time + (1 | id_tmp), data = dat_sub, ...)
            } else {
              lme4::lmer(measure ~ time + (1 | id_tmp), data = dat_sub, weights = weights_mixed, ...)
            }
          },
          error = function(e) NULL
        )

        if (!is.null(fit_try)) {
          tmp <- rep(NA_real_, length(fitted_sub))
          tmp[mixed_keep] <- stats::residuals(fit_try)
          detrended_sub <- tmp
        }
      }
    }

    if (is.null(detrended_sub)) {
      fit_lm <- if (is.null(weights_sub)) {
        stats::lm(fitted_sub ~ time_sub, ...)
      } else {
        stats::lm(fitted_sub ~ time_sub, weights = weights_sub, ...)
      }
      detrended_sub <- stats::residuals(fit_lm)
    }

    out[keep] <- detrended_sub
    out
  }

  .detrend_single_object <- function(object,
                                     per_group = FALSE,
                                     use_mixed_for_detrend = TRUE,
                                     verbose = FALSE,
                                     ...) {

    spec <- .get_spec(object)
    fit_obj <- object@fit
    exp_design <- object@exp_design
    feature_name <- .get_feature_name(object)

    time_var <- spec$time_var
    group_var <- spec$group_var
    id_var <- spec$id_var

    if (verbose) {
      cat("Performing detrending on", feature_name, "...\n")
    }

    fitted_v <- .extract_aligned_fitted(object)
    time_v <- exp_design[[time_var]]
    id_v <- if (!is.null(id_var) && id_var %in% names(exp_design)) exp_design[[id_var]] else NULL
    weights_v <- .get_stored_weights(object)

    if (!is.null(weights_v) && length(weights_v) != nrow(exp_design)) {
      warning("Stored weights length does not match exp_design for feature ", feature_name, ". Ignoring weights.")
      weights_v <- NULL
    }

    group_assign <- rep("all", nrow(exp_design))
    if (per_group) {
      if (!is.null(group_var) && group_var %in% names(exp_design)) {
        group_assign <- as.character(exp_design[[group_var]])
      } else {
        warning("per_group = TRUE but no valid group column was found. Performing detrending without grouping.")
      }
    }

    detrended <- rep(NA_real_, length(fitted_v))

    for (gr in unique(group_assign)) {
      idx <- which(group_assign == gr)

      if (verbose && per_group && gr != "all") {
        cat("  Group:", gr, "\n")
      }

      detrended[idx] <- .safe_detrend_subset(
        fitted_v = fitted_v[idx],
        time_v = time_v[idx],
        id_v = if (!is.null(id_v)) id_v[idx] else NULL,
        weights_v = if (!is.null(weights_v)) weights_v[idx] else NULL,
        use_mixed = use_mixed_for_detrend && inherits(fit_obj, "merMod"),
        ...
      )
    }

    out <- matrix(detrended, ncol = 1)
    colnames(out) <- feature_name
    out
  }

  # ---------------------------------------------------------
  # Main logic
  # ---------------------------------------------------------

  if (verbose) cat("Performing detrending ...\n")
  if (verbose && per_group) cat("Performing detrending separately within groups ...\n")

  if (is(fit, "mixedcirc_fit")) {
    return(
      .detrend_single_object(
        object = fit,
        per_group = per_group,
        use_mixed_for_detrend = use_mixed_for_detrend,
        verbose = verbose,
        ...
      )
    )
  }

  out_list <- vector("list", length(fit))

  for (i in seq_along(fit)) {
    out_list[[i]] <- .detrend_single_object(
      object = fit[i],
      per_group = per_group,
      use_mixed_for_detrend = use_mixed_for_detrend,
      verbose = verbose,
      ...
    )
  }

  out <- do.call(cbind, out_list)
  out
}
