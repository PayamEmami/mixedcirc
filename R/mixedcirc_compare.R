#' Add lagged or lead versions of one circadian signal
#'
#' Creates lagged and/or lead versions of the measured or fitted values from a
#' single `mixedcirc_fit` object.
#'
#' Shifting can be performed either in the raw row order or after sorting by time
#' within user-defined groups such as subject, group, or any other metadata
#' columns stored in the fitted object's `exp_design`.
#'
#' Positive values in `lags` create lagged variables:
#' \preformatted{
#' x_lag1[t] = x[t-1]
#' }
#'
#' Negative values in `lags` create lead variables:
#' \preformatted{
#' x_lead1[t] = x[t+1]
#' }
#'
#' For RRBS fits, the doubled row structure is preserved. The shift is applied at
#' the biological observation level, while keeping the paired `scaler = 1` and
#' `scaler = 0` rows aligned. In other words, both rows belonging to the same
#' biological observation are shifted together.
#'
#' @param x A `mixedcirc_fit` object.
#' @param lags Integer vector of lags / leads to create.
#' @param lag_by Optional character vector of column names in `x@exp_design`
#'   defining the grouping structure used before shifting. If `NULL`, the
#'   function uses the original `id` column when available, otherwise the
#'   original `group` column when available, otherwise no grouping.
#'   To force no grouping, use `lag_by = character(0)`.
#' @param lag_mode One of `"within_sorted"` or `"raw_order"`.
#' @param type One of `"original"` or `"fitted"`.
#' @param value_name Name of the output column containing the unshifted signal.
#' @param keep_order Logical. If `TRUE`, the returned data frame is restored to
#'   the original row order after shift creation.
#'
#' @return
#' A data.frame containing the experimental design, the selected signal column,
#' and one new column per requested lag/lead.
#'
#' @import methods
#' @export
mixedcirc_add_lag <- function(x,
                              lags = 1,
                              lag_by = NULL,
                              lag_mode = c("within_sorted", "raw_order")[1],
                              type = c("original", "fitted")[1],
                              value_name = "measure",
                              keep_order = TRUE) {

  if (!is(x, "mixedcirc_fit")) {
    stop("x must be a mixedcirc_fit object")
  }

  type <- match.arg(type, c("original", "fitted"))
  lag_mode <- match.arg(lag_mode, c("within_sorted", "raw_order"))

  if (length(lags) == 0 || any(is.na(lags)) || any(lags %% 1 != 0)) {
    stop("lags must be a vector of integers")
  }

  lags <- sort(unique(as.integer(lags)))

  if (!is.character(value_name) || length(value_name) != 1 || is.na(value_name) || value_name == "") {
    stop("value_name must be a single non-empty character string")
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
      RRBS = rrbs
    )
  }

  .get_signal_data <- function(object, type, value_name) {
    dat <- object@exp_design

    if (!".used_in_fit" %in% names(dat)) {
      stop("object@exp_design must contain .used_in_fit. Please refit with the updated mixedcirc_detect().")
    }

    dat <- dat[isTRUE(dat$.used_in_fit) | dat$.used_in_fit %in% TRUE, , drop = FALSE]

    if (type == "original") {
      if (!"measure" %in% names(dat)) {
        stop("object@exp_design must contain measure.")
      }
      dat[[value_name]] <- dat$measure
    } else {
      if (!".fitted" %in% names(dat)) {
        stop("object@exp_design must contain .fitted for type = 'fitted'.")
      }
      dat[[value_name]] <- dat$.fitted
    }

    dat
  }

  .shift_vector <- function(vals, k) {
    n <- length(vals)

    if (k > 0) {
      if (n <= k) return(rep(NA_real_, n))
      return(c(rep(NA_real_, k), vals[seq_len(n - k)]))
    }

    if (k < 0) {
      kk <- abs(k)
      if (n <= kk) return(rep(NA_real_, n))
      return(c(vals[(kk + 1):n], rep(NA_real_, kk)))
    }

    vals
  }

  .shift_name <- function(value_name, k) {
    if (k > 0) {
      paste0(value_name, "_lag", k)
    } else if (k < 0) {
      paste0(value_name, "_lead", abs(k))
    } else {
      paste0(value_name, "_lag0")
    }
  }

  dat <- .get_signal_data(x, type = type, value_name = value_name)
  spec <- .get_spec(x)
  time_var <- spec$time_var

  if (is.null(lag_by)) {
    if (!is.null(spec$id_var) && spec$id_var %in% names(dat)) {
      lag_by <- spec$id_var
    } else if (!is.null(spec$group_var) && spec$group_var %in% names(dat)) {
      lag_by <- spec$group_var
    } else {
      lag_by <- character(0)
    }
  }

  if (!is.character(lag_by)) {
    stop("lag_by must be NULL or a character vector of column names")
  }

  lag_by <- unique(lag_by)

  if (length(lag_by) > 0 && !all(lag_by %in% names(dat))) {
    missing_vars <- lag_by[!lag_by %in% names(dat)]
    stop("The following lag_by columns were not found in the data: ",
         paste(missing_vars, collapse = ", "))
  }

  if (lag_mode == "within_sorted" && !time_var %in% names(dat)) {
    stop("The time column '", time_var, "' was not found in the extracted data.")
  }

  dat$.orig_order <- seq_len(nrow(dat))

  if (!spec$RRBS) {
    if (lag_mode == "within_sorted") {
      if (length(lag_by) > 0) {
        ord <- do.call(order, c(dat[, lag_by, drop = FALSE], list(dat[[time_var]], dat$.orig_order)))
      } else {
        ord <- order(dat[[time_var]], dat$.orig_order)
      }
      dat2 <- dat[ord, , drop = FALSE]
    } else {
      dat2 <- dat
    }

    if (length(lag_by) > 0) {
      split_index <- split(
        seq_len(nrow(dat2)),
        interaction(dat2[, lag_by, drop = FALSE], drop = TRUE, lex.order = TRUE)
      )
    } else {
      split_index <- list(all = seq_len(nrow(dat2)))
    }

    xvec <- dat2[[value_name]]

    for (k in lags) {
      out_col <- rep(NA_real_, nrow(dat2))

      for (ii in split_index) {
        vals <- xvec[ii]
        out_col[ii] <- .shift_vector(vals, k)
      }

      dat2[[.shift_name(value_name, k)]] <- out_col
    }

    if (keep_order && lag_mode == "within_sorted") {
      dat2 <- dat2[order(dat2$.orig_order), , drop = FALSE]
    }

    dat2$.orig_order <- NULL
    rownames(dat2) <- NULL
    return(dat2)
  }

  if (!"scaler" %in% names(dat)) {
    stop("RRBS object detected, but `scaler` column is missing from exp_design.")
  }

  base_group_cols <- lag_by
  obs_key_cols <- setdiff(names(dat), c(value_name, "scaler", ".orig_order"))

  dat$.row_id_internal <- seq_len(nrow(dat))
  dat$.pair_id_internal <- ave(
    seq_len(nrow(dat)),
    interaction(dat[, obs_key_cols, drop = FALSE], drop = TRUE, lex.order = TRUE),
    FUN = seq_along
  )

  scaler_levels_present <- sort(unique(dat$scaler))
  if (!all(c(0, 1) %in% scaler_levels_present)) {
    warning("RRBS data do not contain both scaler levels 0 and 1 in the extracted rows.")
  }

  dat1 <- dat[dat$scaler == 1, , drop = FALSE]
  dat0 <- dat[dat$scaler == 0, , drop = FALSE]

  merge_keys <- c(obs_key_cols, ".pair_id_internal")

  paired <- merge(
    dat1[, c(merge_keys, value_name, ".row_id_internal"), drop = FALSE],
    dat0[, c(merge_keys, value_name, ".row_id_internal"), drop = FALSE],
    by = merge_keys,
    suffixes = c(".1", ".0"),
    sort = FALSE
  )

  if (nrow(paired) == 0) {
    stop("Could not reconstruct paired RRBS rows for lagging.")
  }

  obs_dat <- paired[, merge_keys, drop = FALSE]
  obs_dat[[paste0(value_name, ".1")]] <- paired[[paste0(value_name, ".1")]]
  obs_dat[[paste0(value_name, ".0")]] <- paired[[paste0(value_name, ".0")]]
  obs_dat$.row_id_internal.1 <- paired$.row_id_internal.1
  obs_dat$.row_id_internal.0 <- paired$.row_id_internal.0

  if (lag_mode == "within_sorted") {
    if (length(base_group_cols) > 0) {
      ord_obs <- do.call(order, c(obs_dat[, base_group_cols, drop = FALSE], list(obs_dat[[time_var]], obs_dat$.pair_id_internal)))
    } else {
      ord_obs <- order(obs_dat[[time_var]], obs_dat$.pair_id_internal)
    }
    obs_dat2 <- obs_dat[ord_obs, , drop = FALSE]
  } else {
    obs_dat2 <- obs_dat
  }

  if (length(base_group_cols) > 0) {
    split_obs <- split(
      seq_len(nrow(obs_dat2)),
      interaction(obs_dat2[, base_group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    )
  } else {
    split_obs <- list(all = seq_len(nrow(obs_dat2)))
  }

  for (k in lags) {
    shifted_1 <- rep(NA_real_, nrow(obs_dat2))
    shifted_0 <- rep(NA_real_, nrow(obs_dat2))

    for (ii in split_obs) {
      vals1 <- obs_dat2[[paste0(value_name, ".1")]][ii]
      vals0 <- obs_dat2[[paste0(value_name, ".0")]][ii]

      shifted_1[ii] <- .shift_vector(vals1, k)
      shifted_0[ii] <- .shift_vector(vals0, k)
    }

    new_nm <- .shift_name(value_name, k)
    obs_dat2[[paste0(new_nm, ".1")]] <- shifted_1
    obs_dat2[[paste0(new_nm, ".0")]] <- shifted_0
  }

  dat2 <- dat
  for (k in lags) {
    new_nm <- .shift_name(value_name, k)
    dat2[[new_nm]] <- NA_real_

    map1 <- setNames(obs_dat2[[paste0(new_nm, ".1")]], obs_dat2$.row_id_internal.1)
    map0 <- setNames(obs_dat2[[paste0(new_nm, ".0")]], obs_dat2$.row_id_internal.0)

    idx1 <- dat2$scaler == 1
    idx0 <- dat2$scaler == 0

    dat2[[new_nm]][idx1] <- unname(map1[as.character(dat2$.row_id_internal[idx1])])
    dat2[[new_nm]][idx0] <- unname(map0[as.character(dat2$.row_id_internal[idx0])])
  }

  if (keep_order) {
    dat2 <- dat2[order(dat2$.orig_order), , drop = FALSE]
  }

  dat2$.orig_order <- NULL
  dat2$.row_id_internal <- NULL
  dat2$.pair_id_internal <- NULL
  rownames(dat2) <- NULL
  dat2
}


#' Compare two circadian signals using lagged association or correlation models
#'
#' Compares two fitted `mixedcirc_fit` objects after creating lagged / lead
#' versions of `x`.
#'
#' Three analysis modes are supported:
#' \describe{
#'   \item{`"association"`}{Directional lagged regression. The fitted model is
#'   `measure.y ~ shifted measure.x`, or `measure.y ~ shifted measure.x + (1|id)`
#'   when an `id` column is available after matching.}
#'   \item{`"pearson"`}{Pearson correlation between `measure.y` and the shifted
#'   `measure.x`.}
#'   \item{`"rmcorr"`}{Repeated-measures correlation using `rmcorr::rmcorr()`.
#'   This requires an `id` column to still be present after matching.}
#' }
#'
#' The function uses the stored `exp_design` from each fitted object, including
#' the feature-specific `measure`, `.used_in_fit`, and `.fitted` columns created
#' by the updated `mixedcirc_detect()`.
#'
#' @param x A `mixedcirc_fit` object providing the lagged predictor signal.
#' @param y A `mixedcirc_fit` object providing the response signal.
#' @param lags Integer vector of lags / leads to evaluate.
#' @param analysis One of `"association"`, `"pearson"`, or `"rmcorr"`.
#'   Default `"rmcorr"`.
#' @param type One of `"original"` or `"fitted"`.
#' @param separate_groups Logical. If `TRUE` and a group column is still present
#'   after matching, analyses are run separately within each group.
#' @param lag_by Optional character vector of metadata columns defining the lagging
#'   units for `x`.
#' @param lag_mode One of `"within_sorted"` or `"raw_order"`.
#' @param match_by Optional character vector of metadata columns used to match
#'   rows between `x` and `y`.
#'
#'   If `NULL`, the function chooses matching columns automatically from the
#'   shared metadata in the stored experimental designs.
#'
#'   If `match_by` is provided, it is used exactly as supplied for matching, and
#'   no additional columns such as time, group, id, or replicate_id are added
#'   automatically.
#'
#'   This means that if the user wants grouping or repeated-measures structure
#'   to be preserved after matching, those variables must also be included in
#'   `match_by` when needed.
#'
#'   For example:
#'   - if `separate_groups = TRUE`, the grouping variable should usually be part
#'   of `match_by`
#'   - if `analysis = "association"` and a mixed model with `(1|id)` is desired,
#'     the subject identifier should usually be part of `match_by`
#'   - if `analysis = "rmcorr"`, the subject identifier must be part of
#'     `match_by`
#'
#' @param lmer.df Denominator degrees-of-freedom method passed to `emmeans` for
#'   mixed models in `analysis = "association"`.
#' @param return_models Logical. If `TRUE`, fitted objects / test objects are returned.
#' @param return_data Logical. If `TRUE`, the matched analysis data are returned.
#' @param verbose Logical. If `TRUE`, progress messages are printed.
#' @param ... Additional arguments passed to `lm()` or `lmer()` in
#'   `analysis = "association"`.
#'
#' @return
#' A list with at least:
#' \describe{
#'   \item{results}{A data.frame with one row per lag (and per group when relevant).}
#' }
#'
#' If `return_models = TRUE`, the output also contains:
#' \describe{
#'   \item{models}{A named list of fitted model objects or test objects.}
#' }
#'
#' If `return_data = TRUE`, the output also contains:
#' \describe{
#'   \item{data}{The matched analysis data used for the analyses.}
#' }
#'
#' @details
#' `analysis = "association"` is directional, not symmetric. It answers:
#' \preformatted{
#' does shifted x predict y?
#' }
#' Therefore, in general, fitting
#' \preformatted{
#' y ~ shifted x
#' }
#' is not equivalent to fitting
#' \preformatted{
#' x ~ shifted y
#' }
#'
#' If the reverse direction is of interest, swap `x` and `y` and run the
#' function again.
#'
#' Please note that the association analysis might fail for identical time series.
#'
#' `analysis = "pearson"` and `analysis = "rmcorr"` are correlation-style
#' summaries rather than directional regression models.
#'
#' For RRBS fits, the signal is collapsed before analysis to:
#' \preformatted{
#' scaler = 1 - scaler = 0
#' }
#'
#' @import stats
#' @import methods
#' @import emmeans
#' @import performance
#' @import MuMIn
#' @import rmcorr
#' @import lme4
#' @export
mixedcirc_compare <- function(x,
                              y,
                              lags = -1:1,
                              analysis = c("association", "pearson", "rmcorr")[3],
                              type = c("original", "fitted")[1],
                              separate_groups = TRUE,
                              lag_by = NULL,
                              lag_mode = c("within_sorted", "raw_order")[1],
                              match_by = NULL,
                              lmer.df = c("Satterthwaite", "Kenward-Roger")[1],
                              return_models = FALSE,
                              return_data = FALSE,
                              verbose = FALSE,
                              ...) {

  if (missing(x)) stop("x must be provided")
  if (missing(y)) stop("y must be provided")
  if (!is(x, "mixedcirc_fit")) stop("x must be a mixedcirc_fit object")
  if (!is(y, "mixedcirc_fit")) stop("y must be a mixedcirc_fit object")

  analysis <- match.arg(analysis, c("association", "pearson", "rmcorr"))
  type <- match.arg(type, c("original", "fitted"))
  lag_mode <- match.arg(lag_mode, c("within_sorted", "raw_order"))
  lmer.df <- match.arg(lmer.df, c("Satterthwaite", "Kenward-Roger"))

  if (length(lags) == 0 || any(is.na(lags)) || any(lags %% 1 != 0)) {
    stop("lags must be a vector of integers")
  }
  lags <- sort(unique(as.integer(lags)))

  .get_call_arg_name <- function(object, arg_name) {
    if (!"call" %in% methods::slotNames(object)) return(NULL)
    cl <- object@call
    if (is.null(cl) || is.null(cl[[arg_name]])) return(NULL)
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
      RRBS = rrbs
    )
  }

  .collapse_rrbs_difference <- function(dat,
                                        value_col = "value",
                                        scaler_col = "scaler") {

    dat1 <- dat[dat[[scaler_col]] == 1, , drop = FALSE]
    dat0 <- dat[dat[[scaler_col]] == 0, , drop = FALSE]

    if (nrow(dat1) == 0 || nrow(dat0) == 0) {
      stop("For RRBS modeling, both scaler == 1 and scaler == 0 must exist.")
    }

    key_cols <- setdiff(names(dat), c(value_col, scaler_col))

    dat1$.pair_id <- ave(
      seq_len(nrow(dat1)),
      interaction(dat1[, key_cols, drop = FALSE], drop = TRUE),
      FUN = seq_along
    )

    dat0$.pair_id <- ave(
      seq_len(nrow(dat0)),
      interaction(dat0[, key_cols, drop = FALSE], drop = TRUE),
      FUN = seq_along
    )

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

  .get_analysis_data <- function(object, type, value_name) {
    dat <- object@exp_design
    spec <- .get_spec(object)

    if (!".used_in_fit" %in% names(dat)) {
      stop("object@exp_design must contain .used_in_fit. Please refit with the updated mixedcirc_detect().")
    }

    dat <- dat[dat$.used_in_fit %in% TRUE, , drop = FALSE]

    if (type == "original") {
      if (!"measure" %in% names(dat)) {
        stop("object@exp_design must contain measure.")
      }
      dat[[value_name]] <- dat$measure
    } else {
      if (!".fitted" %in% names(dat)) {
        stop("object@exp_design must contain .fitted for type = 'fitted'.")
      }
      dat[[value_name]] <- dat$.fitted
    }

    if (spec$RRBS && "scaler" %in% names(dat)) {
      dat <- .collapse_rrbs_difference(dat, value_col = value_name, scaler_col = "scaler")
    }

    dat
  }

  .default_match_by <- function(x_dat, y_dat, spec_x, spec_y) {
    shared <- intersect(names(x_dat), names(y_dat))
    shared <- setdiff(shared, c("measure.x", "measure.y"))

    preferred <- c(
      spec_x$time_var,
      spec_x$group_var,
      spec_x$id_var,
      spec_x$replicate_var,
      spec_y$time_var,
      spec_y$group_var,
      spec_y$id_var,
      spec_y$replicate_var
    )
    preferred <- unique(preferred[!is.null(preferred) & !is.na(preferred)])
    preferred <- intersect(preferred, shared)

    if (length(preferred) > 0) {
      return(preferred)
    }

    shared <- setdiff(shared, c("inphase", "outphase", "scaler", ".used_in_fit",
                                ".fitted", ".residual", ".weight"))
    if (length(shared) == 0) {
      stop("Could not determine matching columns automatically. Please provide match_by.")
    }
    shared
  }

  .get_effect_from_emmeans <- function(model_obj, pred_name, lmer.df) {
    em_obj <- emmeans::emtrends(
      model_obj,
      ~ 1,
      var = pred_name,
      lmer.df = lmer.df
    )

    em_df <- as.data.frame(em_obj)
    trend_col <- grep("\\.trend$", names(em_df), value = TRUE)[1]
    test_df <- as.data.frame(emmeans::test(em_obj, null = 0))
    stat_col <- grep("(ratio|statistic|F\\.ratio|t\\.ratio)$", names(test_df), value = TRUE)[1]

    list(
      estimate = em_df[[trend_col]][1],
      std.error = em_df$SE[1],
      df = if ("df" %in% names(em_df)) em_df$df[1] else NA_real_,
      statistic = if (!is.na(stat_col)) test_df[[stat_col]][1] else NA_real_,
      p.value = if ("p.value" %in% names(test_df)) test_df$p.value[1] else NA_real_
    )
  }

  .get_r2_lm <- function(model_obj) {
    s <- summary(model_obj)
    c(
      r.squared = unname(s$r.squared),
      r.squared.marginal = unname(s$r.squared),
      r.squared.conditional = unname(s$r.squared)
    )
  }

  .get_r2_lmer <- function(model_obj) {
    if (requireNamespace("performance", quietly = TRUE)) {
      rr <- performance::r2_nakagawa(model_obj)
      return(c(
        r.squared = unname(rr$R2_marginal),
        r.squared.marginal = unname(rr$R2_marginal),
        r.squared.conditional = unname(rr$R2_conditional)
      ))
    }

    if (requireNamespace("MuMIn", quietly = TRUE)) {
      rr <- MuMIn::r.squaredGLMM(model_obj)
      return(c(
        r.squared = unname(rr[1, 1]),
        r.squared.marginal = unname(rr[1, 1]),
        r.squared.conditional = unname(rr[1, 2])
      ))
    }

    fixed_eff <- lme4::fixef(model_obj)
    fittedvar <- stats::var(as.vector(stats::model.matrix(model_obj) %*% fixed_eff), na.rm = TRUE)
    vc <- as.data.frame(lme4::VarCorr(model_obj))
    residvar <- vc$vcov[vc$grp == "Residual"][1]
    randvar <- sum(vc$vcov[vc$grp != "Residual" & is.na(vc$var2)], na.rm = TRUE)

    marginal <- fittedvar / (fittedvar + randvar + residvar)
    conditional <- (fittedvar + randvar) / (fittedvar + randvar + residvar)

    c(
      r.squared = marginal,
      r.squared.marginal = marginal,
      r.squared.conditional = conditional
    )
  }

  .fit_one_association <- function(dat_sub, pred_name, id_var, lmer.df, ...) {
    dat_sub <- dat_sub[!is.na(dat_sub[[pred_name]]) & !is.na(dat_sub$measure.y), , drop = FALSE]

    if (nrow(dat_sub) < 3) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = NA_character_, note = "Too few complete observations"))
    }

    use_lmer <- FALSE
    if (!is.null(id_var) && id_var %in% names(dat_sub)) {
      id_nonmis <- dat_sub[[id_var]][!is.na(dat_sub[[id_var]])]
      if (length(unique(id_nonmis)) > 1) {
        use_lmer <- TRUE
      }
    }

    if (use_lmer) {
      fml <- stats::as.formula(paste0("measure.y ~ ", pred_name, " + (1 | ", id_var, ")"))
      model_obj <- tryCatch(
        lme4::lmer(fml, data = dat_sub, ...),
        error = function(e) NULL
      )

      if (!is.null(model_obj)) {
        eff <- .get_effect_from_emmeans(model_obj, pred_name, lmer.df = lmer.df)
        r2 <- .get_r2_lmer(model_obj)
        return(list(
          model = model_obj,
          stats = c(eff, as.list(r2)),
          n.used = stats::nobs(model_obj),
          method = "lmer",
          note = ""
        ))
      }
    }

    fml <- stats::as.formula(paste0("measure.y ~ ", pred_name))
    model_obj <- tryCatch(
      stats::lm(fml, data = dat_sub, ...),
      error = function(e) NULL
    )

    if (is.null(model_obj)) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = NA_character_, note = "Model fitting failed"))
    }

    eff <- .get_effect_from_emmeans(model_obj, pred_name, lmer.df = lmer.df)
    r2 <- .get_r2_lm(model_obj)

    list(
      model = model_obj,
      stats = c(eff, as.list(r2)),
      n.used = stats::nobs(model_obj),
      method = "lm",
      note = if (use_lmer) "Mixed model failed; used lm instead" else ""
    )
  }

  .fit_one_pearson <- function(dat_sub, pred_name) {
    dat_sub <- dat_sub[!is.na(dat_sub[[pred_name]]) & !is.na(dat_sub$measure.y), , drop = FALSE]

    if (nrow(dat_sub) < 3) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "pearson", note = "Too few complete observations"))
    }

    tst <- tryCatch(
      stats::cor.test(dat_sub[[pred_name]], dat_sub$measure.y, method = "pearson"),
      error = function(e) NULL
    )

    if (is.null(tst)) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "pearson", note = "Correlation test failed"))
    }

    r <- unname(tst$estimate[1])

    list(
      model = tst,
      stats = list(
        estimate = r,
        std.error = NA_real_,
        df = if (!is.null(tst$parameter)) unname(tst$parameter[1]) else NA_real_,
        statistic = if (!is.null(tst$statistic)) unname(tst$statistic[1]) else NA_real_,
        p.value = tst$p.value,
        r.squared = r^2,
        r.squared.marginal = r^2,
        r.squared.conditional = r^2
      ),
      n.used = sum(stats::complete.cases(dat_sub[, c(pred_name, "measure.y")])),
      method = "pearson",
      note = ""
    )
  }

  .fit_one_rmcorr <- function(dat_sub, pred_name, id_var) {
    dat_sub <- dat_sub[!is.na(dat_sub[[pred_name]]) &
                         !is.na(dat_sub$measure.y) &
                         !is.na(dat_sub[[id_var]]), , drop = FALSE]

    if (!requireNamespace("rmcorr", quietly = TRUE)) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "rmcorr", note = "Package 'rmcorr' is not installed"))
    }

    if (nrow(dat_sub) < 3) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "rmcorr", note = "Too few complete observations"))
    }

    if (is.null(id_var) || !id_var %in% names(dat_sub)) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "rmcorr", note = "id column is required for rmcorr"))
    }

    if (length(unique(dat_sub[[id_var]])) < 3) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "rmcorr", note = "rmcorr requires at least 3 subjects"))
    }

    rm_obj <- tryCatch({
      eval(
        substitute(
          rmcorr::rmcorr(participant = ID, measure1 = X, measure2 = Y, dataset = DATA),
          list(
            ID = as.name(id_var),
            X = as.name(pred_name),
            Y = as.name("measure.y"),
            DATA = dat_sub
          )
        )
      )
    }, error = function(e) e)

    if (inherits(rm_obj, "error")) {
      return(list(model = NULL, stats = NULL, n.used = nrow(dat_sub), method = "rmcorr", note = conditionMessage(rm_obj)))
    }

    r <- unname(rm_obj$r)

    list(
      model = rm_obj,
      stats = list(
        estimate = r,
        std.error = NA_real_,
        df = if (!is.null(rm_obj$df)) unname(rm_obj$df) else NA_real_,
        statistic = NA_real_,
        p.value = if (!is.null(rm_obj$p)) rm_obj$p else NA_real_,
        r.squared = r^2,
        r.squared.marginal = r^2,
        r.squared.conditional = r^2
      ),
      n.used = nrow(dat_sub),
      method = "rmcorr",
      note = ""
    )
  }

  .shift_name <- function(value_name, k) {
    if (k > 0) {
      paste0(value_name, "_lag", k)
    } else if (k < 0) {
      paste0(value_name, "_lead", abs(k))
    } else {
      paste0(value_name, "_lag0")
    }
  }

  spec_x <- .get_spec(x)
  spec_y <- .get_spec(y)

  x_signal <- .get_analysis_data(x, type = type, value_name = "measure.x")
  y_signal <- .get_analysis_data(y, type = type, value_name = "measure.y")

  if (separate_groups && !is.null(spec_x$group_var)) {
    if (is.null(lag_by)) {
      lag_by_eff <- spec_x$group_var
      if (!is.null(spec_x$id_var) && spec_x$id_var %in% names(x_signal)) {
        lag_by_eff <- unique(c(lag_by_eff, spec_x$id_var))
      }
    } else {
      lag_by_eff <- unique(c(spec_x$group_var, lag_by))
    }
  } else if (is.null(lag_by)) {
    if (!is.null(spec_x$id_var) && spec_x$id_var %in% names(x_signal)) {
      lag_by_eff <- spec_x$id_var
    } else {
      lag_by_eff <- character(0)
    }
  } else {
    lag_by_eff <- unique(lag_by)
  }

  x_lagged <- mixedcirc_add_lag(
    x = x,
    lags = lags,
    lag_by = lag_by_eff,
    lag_mode = lag_mode,
    type = type,
    value_name = "measure.x",
    keep_order = TRUE
  )

  if (is.null(match_by)) {
    match_by_eff <- .default_match_by(x_lagged, y_signal, spec_x, spec_y)
  } else {
    if (!is.character(match_by)) stop("match_by must be NULL or a character vector")
    match_by_eff <- unique(match_by)
  }

  missing_match_x <- setdiff(match_by_eff, names(x_lagged))
  missing_match_y <- setdiff(match_by_eff, names(y_signal))
  if (length(missing_match_x) > 0) {
    stop("The following match_by columns were not found in x analysis data: ",
         paste(missing_match_x, collapse = ", "))
  }
  if (length(missing_match_y) > 0) {
    stop("The following match_by columns were not found in y analysis data: ",
         paste(missing_match_y, collapse = ", "))
  }

  # IMPORTANT:
  # match_by is used exactly as supplied. No extra columns are added automatically.
  x_keep <- unique(c(
    match_by_eff,
    "measure.x",
    vapply(lags, function(k) .shift_name("measure.x", k), character(1))
  ))
  x_keep <- intersect(x_keep, names(x_lagged))

  y_keep <- unique(c(match_by_eff, "measure.y"))
  y_keep <- intersect(y_keep, names(y_signal))

  data_for_model <- merge(
    x_lagged[, x_keep, drop = FALSE],
    y_signal[, y_keep, drop = FALSE],
    by = match_by_eff,
    sort = FALSE
  )

  if (nrow(data_for_model) == 0) {
    stop("No matched observations were found between x and y.")
  }

  # group/id are only available here if the user preserved them via match_by
  group_var <- spec_x$group_var
  id_var <- spec_x$id_var

  if (!is.null(group_var) && !group_var %in% names(data_for_model)) {
    if (separate_groups && verbose) {
      message("Group column not present after matching; running pooled analysis.")
    }
    group_var <- NULL
  }

  if (!is.null(id_var) && !id_var %in% names(data_for_model)) {
    id_var <- NULL
  }

  if (is.null(group_var) || !separate_groups) {
    group_levels <- "all"
    data_for_model$.analysis_group <- "all"
    group_var_out <- NULL
  } else {
    data_for_model$.analysis_group <- as.character(data_for_model[[group_var]])
    group_levels <- unique(data_for_model$.analysis_group)
    group_var_out <- group_var
  }

  if (verbose) cat("Running lagged", analysis, "models ...\n")

  results_list <- list()
  models_list <- list()
  counter <- 1L

  for (gr in group_levels) {
    dat_g <- data_for_model[data_for_model$.analysis_group == gr, , drop = FALSE]

    if (verbose) cat("Processing group:", gr, "\n")

    for (k in lags) {
      pred_name <- .shift_name("measure.x", k)

      fit_res <- switch(
        analysis,
        association = .fit_one_association(
          dat_sub = dat_g,
          pred_name = pred_name,
          id_var = id_var,
          lmer.df = lmer.df,
          ...
        ),
        pearson = .fit_one_pearson(
          dat_sub = dat_g,
          pred_name = pred_name
        ),
        rmcorr = .fit_one_rmcorr(
          dat_sub = dat_g,
          pred_name = pred_name,
          id_var = id_var
        )
      )

      row_out <- data.frame(
        group = if (!is.null(group_var_out)) gr else NA_character_,
        lag = k,
        predictor = pred_name,
        analysis = analysis,
        method = fit_res$method,
        estimate_type = if (analysis == "association") "slope" else "correlation",
        n.used = fit_res$n.used,
        estimate = NA_real_,
        std.error = NA_real_,
        df = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        r.squared = NA_real_,
        r.squared.marginal = NA_real_,
        r.squared.conditional = NA_real_,
        note = fit_res$note,
        stringsAsFactors = FALSE
      )

      if (!is.null(fit_res$stats)) {
        row_out$estimate <- fit_res$stats$estimate
        row_out$std.error <- fit_res$stats$std.error
        row_out$df <- fit_res$stats$df
        row_out$statistic <- fit_res$stats$statistic
        row_out$p.value <- fit_res$stats$p.value
        row_out$r.squared <- fit_res$stats$r.squared
        row_out$r.squared.marginal <- fit_res$stats$r.squared.marginal
        row_out$r.squared.conditional <- fit_res$stats$r.squared.conditional
      }

      results_list[[counter]] <- row_out

      if (return_models) {
        nm <- paste0("group_", gr, "__lag_", k)
        models_list[[nm]] <- fit_res$model
      }

      counter <- counter + 1L
    }
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  out <- list(results = results_df)

  if (return_models) out$models <- models_list
  if (return_data) out$data <- data_for_model

  out
}
