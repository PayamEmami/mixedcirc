#' Remove linear trend and nuisance covariate effects from the data
#'
#' Fits one regression model per variable and returns the residuals after removing
#' a non-circadian linear trend. By default, the removed trend is the linear effect
#' of `time`. Additional nuisance covariates can be supplied through
#' `formula_extra`.
#'
#' The function uses the same metadata-driven interface as `mixedcirc_detect()`.
#' For RRBS / BS-seq style data, the model respects the paired methylated /
#' unmethylated structure using `scaler` and `replicate_id`.
#'
#' @param data_input A numeric matrix, data.frame, or `DGEList` containing the input data.
#'   Rows correspond to observations and columns correspond to variables when `data_input`
#'   is a matrix or data.frame. If a `DGEList` is provided, rows correspond to variables and
#'   columns correspond to samples.
#' @param meta_input A data.frame or matrix containing the sample-level metadata used to build
#'   the detrending design. This must include the columns referenced by `time`, `group`, `id`,
#'   `replicate_id`, and any variables appearing in `formula_extra`.
#' @param time A character string giving the name of the column in `meta_input` that contains
#'   time.
#' @param group Optional character string giving the name of a grouping column in `meta_input`.
#'   This is not removed by default. If `remove_trend_separate_groups = TRUE`, detrending is
#'   performed separately within each group.
#' @param id Optional character string giving the name of the subject identifier column in
#'   `meta_input`. This is used for random effects when `lm_method = "lme"`. If `NULL`, the
#'   function switches to ordinary linear regression (`lm`).
#' @param lm_method Regression framework to use. Must be one of `"lm"` or `"lme"`.
#'   `"lm"` fits an ordinary linear model. `"lme"` fits a mixed-effects model using
#'   `lme4::lmer()`. Default is `"lme"`.
#' @param formula_extra Optional formula containing additional nuisance fixed effects to remove.
#'   The variables in this formula must be column names in `meta_input`. For RRBS data, these
#'   extra terms are automatically interacted with `scaler`.
#' @param random_effect_type Type of random-effect structure to use when `lm_method = "lme"`.
#'   Must be one of `"none"` or `"mesor"`. `"none"` adds no random effect and `"mesor"` adds
#'   a random intercept `(1 | id)`. Default is `"mesor"`.
#' @param obs_weights Optional matrix of observational weights with the same dimensions as
#'   `data_input` when `data_input` is a matrix or data.frame. Each column corresponds to one
#'   variable and each row to one observation. For RRBS input, `obs_weights` is required unless
#'   weights are estimated internally from a `DGEList` or by setting
#'   `force_weight_estimation = TRUE`.
#' @param RRBS Logical. If `TRUE`, the data are assumed to represent RRBS / BS-seq style
#'   methylation data, with each biological observation represented by two consecutive rows:
#'   methylated count followed by unmethylated count. In this mode, `replicate_id` must be
#'   provided. Default is `FALSE`.
#' @param replicate_id Character string giving the name of the column in `meta_input` that
#'   identifies unique technical / replicate observations for RRBS data.
#' @param remove_trend_separate_groups Logical. If `TRUE` and group information is provided,
#'   detrending is performed separately within each group. Default is `TRUE`.
#' @param force_weight_estimation Logical. If `TRUE`, weights are estimated internally using
#'   `voom()` or `voomWithDreamWeights()` even when `data_input` is not a `DGEList`.
#'   Default is `FALSE`.
#' @param ncores Number of cores to use for parallel processing. Default is `1`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default is `FALSE`.
#' @param ... Additional arguments passed to the underlying model-fitting function.
#'
#' @return
#' A numeric matrix of detrended values. Rows correspond to observations and columns correspond
#' to variables. The values are the residuals after removing the fitted trend model.
#'
#' @details
#' For non-RRBS data, the default detrending model is:
#' \preformatted{
#' measure ~ 1 + time
#' }
#'
#' If `formula_extra` is supplied, its terms are appended to the model. These terms are treated
#' as nuisance covariates to remove.
#'
#' For RRBS data, the default detrending model is:
#' \preformatted{
#' measure ~ 0 + replicate_id + scaler + time:scaler
#' }
#'
#' If `formula_extra` is supplied for RRBS data, the extra terms are automatically interacted
#' with `scaler`, so that the nuisance structure is removed on the methylation scale.
#'
#' If `remove_trend_separate_groups = TRUE`, the detrending model is fit separately within each
#' group. Otherwise, a single common detrending model is fit across all groups. Group itself is
#' not removed unless the user explicitly includes it in `formula_extra`.
#'
#' If weights are not supplied directly, they may be estimated from count-like data using
#' `limma::voom()` for ordinary linear models or `voomWithDreamWeights()`
#' for mixed-effects models.
#'
#' @examples
#' data("circa_data")
#'
#' detrended <- mixedcirc_remove_trend(
#'   data_input = circa_data$data_matrix,
#'   meta_input = data.frame(
#'     time = circa_data$time,
#'     group = circa_data$group,
#'     id = circa_data$id
#'   ),
#'   time = "time",
#'   group = "group",
#'   id = "id",
#'   verbose = TRUE
#' )
#'
#' @import stats
#' @import methods
#' @import doFuture
#' @import future
#' @import foreach
#' @import lme4
#' @import limma
#' @import parallel
#' @import BiocParallel
#' @export
mixedcirc_remove_trend <- function(
    data_input = NULL,
    meta_input = NULL,
    time = NULL,
    group = NULL,
    id = NULL,
    lm_method = c("lm", "lme")[2],
    formula_extra = NULL,
    random_effect_type = c("none", "mesor")[2],
    obs_weights = NULL,
    RRBS = FALSE,
    replicate_id = NULL,
    remove_trend_separate_groups = TRUE,
    force_weight_estimation = FALSE,
    ncores = 1,
    verbose = FALSE,
    ...
) {

  registerDoFuture()
  if (ncores > 1) {
    plan(multisession, workers = ncores)
  } else {
    plan(sequential)
  }

  on.exit({
    future:::ClusterRegistry("stop")
  }, add = TRUE)

  if (verbose) cat("Checking inputs ...\n")

  # ---------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------
  if (is.null(data_input)) {
    stop("data_input must be a data frame, matrix, or DGEList")
  }

  if (is.null(meta_input)) {
    stop("meta_input must be a data frame or matrix")
  }

  if (!is.matrix(data_input) && !is.data.frame(data_input) && !is(data_input, "DGEList")) {
    stop("data_input must be a data frame, matrix, or DGEList")
  }

  if (!is.matrix(meta_input) && !is.data.frame(meta_input)) {
    stop("meta_input must be a data frame or matrix")
  }

  meta_input <- as.data.frame(meta_input)

  if (is.null(time) || !time %in% colnames(meta_input)) {
    stop("time must be a column name in meta_input")
  }

  if (!is(data_input, "DGEList")) {
    if (RRBS) {
      if (nrow(meta_input) * 2 != nrow(data_input)) {
        stop("For RRBS, meta_input must have one row per biological sample, so nrow(data_input) must equal 2 * nrow(meta_input)")
      }
    } else {
      if (nrow(meta_input) != nrow(data_input)) {
        stop("The number of rows of meta_input must equal the number of rows of data_input")
      }
    }
  } else {
    if (nrow(meta_input) != ncol(data_input)) {
      stop("The number of rows of meta_input must equal the number of columns of data_input for DGEList input")
    }
  }

  multiple_groups <- FALSE
  if (!is.null(group)) {
    if (!group %in% colnames(meta_input)) {
      stop("group must be the name of a column in meta_input")
    }
    multiple_groups <- TRUE
  } else {
    meta_input$dummy_group <- "all"
    group <- "dummy_group"
  }

  if (!is.null(id)) {
    if (!id %in% colnames(meta_input)) {
      stop("id must be the name of a column in meta_input")
    }
  } else {
    warning("No id information was provided! Switching to normal linear regression!")
    lm_method <- "lm"
    meta_input$dummy_id <- "all"
    id <- "dummy_id"
  }

  lm_method <- match.arg(lm_method, c("lm", "lme"))
  random_effect_type <- match.arg(random_effect_type, c("none", "mesor"))

  if (RRBS) {
    if (is.null(replicate_id) || !replicate_id %in% colnames(meta_input)) {
      stop("For RRBS data, replicate_id must be set and must be a column name in meta_input")
    }

    if (is.null(obs_weights) && !is(data_input, "DGEList") && !force_weight_estimation) {
      stop("For RRBS data, obs_weights must be provided unless weights are estimated internally from a DGEList or by setting force_weight_estimation = TRUE")
    }
  }

  if (!is.null(obs_weights) && !is(data_input, "DGEList")) {
    if (!all(dim(obs_weights) == dim(data_input))) {
      stop("obs_weights must be the same size as data_input")
    }
  }

  if (!is.null(formula_extra)) {
    if (!inherits(formula_extra, "formula")) {
      stop("formula_extra must be a formula")
    }

    missing_formula_extra_terms <- setdiff(all.vars(formula_extra), colnames(meta_input))
    if (length(missing_formula_extra_terms) > 0) {
      stop(
        "formula_extra terms must be column names in meta_input: ",
        paste(missing_formula_extra_terms, collapse = ", ")
      )
    }
  }

  if (!is(data_input, "DGEList")) {
    eset <- data_input[, , drop = FALSE]
  } else {
    eset <- NA
  }

  # ---------------------------------------------------------
  # Build experiment design
  # ---------------------------------------------------------
  if (verbose) cat("Building experiment ...\n")

  exp_design <- meta_input

  if (RRBS) {
    exp_design <- meta_input[rep(seq_len(nrow(meta_input)), each = 2), , drop = FALSE]
    exp_design$scaler <- rep(c(1, 0), nrow(meta_input))
  }

  exp_design[, group] <- factor(exp_design[, group])

  # ---------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------
  build_random_term <- function(random_effect_type, id) {
    if (random_effect_type == "none") return(NULL)
    if (random_effect_type == "mesor") return(paste0("(1 | ", id, ")"))
    stop("Unsupported random_effect_type")
  }

  append_extra_terms <- function(base_formula, formula_extra, RRBS = FALSE) {
    if (is.null(formula_extra)) {
      return(base_formula)
    }

    rhs1 <- paste(deparse(base_formula[[3]]), collapse = "")

    if (RRBS) {
      formula_extra_mod <- update(formula_extra, . ~ scaler:(.))
      rhs2 <- paste(deparse(formula_extra_mod[[3]]), collapse = "")
    } else {
      rhs2 <- paste(deparse(formula_extra[[2]]), collapse = "")
    }

    as.formula(
      paste(all.vars(base_formula)[1], "~", paste(rhs1, rhs2, sep = " + "))
    )
  }

  add_random_effect <- function(formula_fixed, rand_term) {
    if (is.null(rand_term)) {
      return(formula_fixed)
    }

    rhs <- paste(paste(deparse(formula_fixed[[3]]), collapse = ""), rand_term, sep = " + ")
    as.formula(paste(all.vars(formula_fixed)[1], "~", rhs))
  }

  get_model_matrix_formula <- function(formula_fixed) {
    as.formula(paste("~", paste(deparse(formula_fixed[[3]]), collapse = "")))
  }

  safe_fit_model <- function(formula_use, data_use, weights_use, lm_method, ...) {
    if (lm_method == "lm") {
      if (is.null(weights_use)) {
        stats::lm(formula_use, data = data_use, ...)
      } else {
        stats::lm(formula_use, data = data_use, weights = weights_use, ...)
      }
    } else {
      if (is.null(weights_use)) {
        lme4::lmer(formula_use, data = data_use, ...)
      } else {
        lme4::lmer(formula_use, data = data_use, weights = weights_use, ...)
      }
    }
  }

  fill_residuals_back <- function(model_obj, data_rows) {
    out <- rep(NA_real_, length(data_rows))
    mf <- stats::model.frame(model_obj)

    mf_rows <- rownames(mf)
    data_rows <- as.character(data_rows)

    if (is.null(mf_rows) || is.null(data_rows)) {
      if (length(stats::residuals(model_obj)) != length(out)) {
        stop("Could not align model residuals back to the supplied data rows.")
      }
      out <- stats::residuals(model_obj)
      return(out)
    }

    pos <- match(mf_rows, data_rows)

    if (any(is.na(pos))) {
      stop("Could not align model-frame rows back to the supplied data rows.")
    }

    out[pos] <- stats::residuals(model_obj)
    out
  }

  # ---------------------------------------------------------
  # Build detrending formula
  # ---------------------------------------------------------
  if (verbose) cat("Building detrending formula ...\n")

  if (RRBS) {
    formula_input <- as.formula(
      paste0(
        "measure ~ 0 + ", replicate_id,
        " + scaler + ", time, ":scaler"
      )
    )
  } else {
    formula_input <- as.formula(
      paste0("measure ~ 1 + ", time)
    )
  }

  formula_fixed <- append_extra_terms(
    base_formula = formula_input,
    formula_extra = formula_extra,
    RRBS = RRBS
  )

  formula_with_random <- formula_fixed
  if (lm_method == "lme") {
    rand_term <- build_random_term(
      random_effect_type = random_effect_type,
      id = id
    )
    formula_with_random <- add_random_effect(formula_fixed, rand_term)
  }

  if (lm_method == "lm" && random_effect_type != "none") {
    warning("random_effect_type was ignored because lm_method = 'lm'")
  }

  if (verbose) {
    cat("Final detrending formula:\n")
    print(formula_with_random)
  }

  model_matrix_formula <- get_model_matrix_formula(formula_fixed)
  model_matrix <- stats::model.matrix(model_matrix_formula, exp_design)

  # ---------------------------------------------------------
  # Weight estimation
  # ---------------------------------------------------------
  if (is(data_input, "DGEList") || force_weight_estimation) {
    if (verbose) cat("Estimating weights for the input variables ...\n")

    if (lm_method == "lm") {
      if (verbose) cat("lm_method is lm, voom will be used ...\n")
      voom_res <- limma::voom(data_input, model_matrix, plot = FALSE)
    } else {
      voom_res <- voomWithDreamWeights(
        data_input,
        formula = formula_with_random,
        data = exp_design,
        BPPARAM = BiocParallel::SnowParam(workers = ncores)
      )
    }

    obs_weights <- t(voom_res$weights)
    eset <- t(voom_res$E)
  }

  feature_names <- colnames(eset)
  if (length(feature_names) == 0) {
    feature_names <- seq_len(ncol(eset))
  }

  chunks <- parallel::splitIndices(ncol(eset), min(ncol(eset), ncores))
  if (verbose) cat("Splitting data into", length(chunks), "chunks ...\n")

  # ---------------------------------------------------------
  # Main fitting loop
  # ---------------------------------------------------------
  res <- foreach(chIndx = chunks) %dopar% {

    outputs_fn <- foreach(i = chIndx) %do% {

      if (verbose) cat("Processing variable", i, "...\n")

      data_grouped <- cbind(
        measure = as.numeric(eset[, i]),
        exp_design
      )
      data_grouped <- as.data.frame(data_grouped)

      if (is.null(rownames(data_grouped)) ||
          anyNA(rownames(data_grouped)) ||
          any(rownames(data_grouped) == "") ||
          anyDuplicated(rownames(data_grouped))) {
        rownames(data_grouped) <- paste0("row_", seq_len(nrow(data_grouped)))
      }

      weights_i <- NULL
      if (!is.null(obs_weights)) {
        weights_i <- obs_weights[, i]
      }

      residual_out <- rep(NA_real_, nrow(data_grouped))

      if (multiple_groups && remove_trend_separate_groups) {
        group_levels <- unique(as.character(data_grouped[, group]))

        for (gr in group_levels) {
          idx <- which(as.character(data_grouped[, group]) == gr)
          data_sub <- data_grouped[idx, , drop = FALSE]
          weights_sub <- if (!is.null(weights_i)) weights_i[idx] else NULL

          if (is.null(rownames(data_sub)) ||
              anyNA(rownames(data_sub)) ||
              any(rownames(data_sub) == "") ||
              anyDuplicated(rownames(data_sub))) {
            rownames(data_sub) <- rownames(data_grouped)[idx]
          }

          model_ln <- suppressMessages(
            safe_fit_model(
              formula_use = formula_with_random,
              data_use = data_sub,
              weights_use = weights_sub,
              lm_method = lm_method,
              ...
            )
          )

          residual_out[idx] <- fill_residuals_back(model_ln, rownames(data_sub))
        }
      } else {
        model_ln <- suppressMessages(
          safe_fit_model(
            formula_use = formula_with_random,
            data_use = data_grouped,
            weights_use = weights_i,
            lm_method = lm_method,
            ...
          )
        )

        residual_out <- fill_residuals_back(model_ln, rownames(data_grouped))
      }

      residual_out
    }

    outputs_fn
  }

  results_out <- do.call("cbind", lapply(rapply(res, enquote, how = "unlist"), eval))
  colnames(results_out) <- feature_names
  rownames(results_out) <- rownames(eset)

  results_out
}
