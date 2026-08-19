# Feature-wise models and mixed-model voom-style weights ---------------------
# This package-local implementation keeps the RNA-seq and repeated-measures
# workflows independent of an external variance-decomposition package.

.mixedcirc_vp_control <- lme4::lmerControl(
  calc.derivs = FALSE,
  check.rankX = "stop.deficient",
  check.conv.singular = list(action = "ignore", tol = 1e-4)
)

methods::setClass(
  "VarParFitList",
  slots = c(method = "character"),
  contains = "list"
)

methods::setClass(
  "varPartResults",
  slots = c(type = "character", method = "character"),
  contains = "data.frame"
)

.mc_scalar_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  value
}

# lme4 evaluates its `weights` argument non-standardly in the formula
# environment.  When a caller has a global object with that exact name, some
# versions can replace the binding while constructing the model frame.  Public
# entry points preserve the caller's binding explicitly.
.mc_capture_global_weights <- function() {
  present <- exists("weights", envir = .GlobalEnv, inherits = FALSE)
  list(
    present = present,
    value = if (present) get("weights", envir = .GlobalEnv, inherits = FALSE)
  )
}

.mc_restore_global_weights <- function(state) {
  if (state$present) {
    caller_globals <- globalenv()
    caller_globals[["weights"]] <- state$value
  } else if (exists("weights", envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = "weights", envir = .GlobalEnv)
  }
  invisible(NULL)
}

.mc_is_mixed_formula <- function(formula) {
  formula <- stats::as.formula(formula)
  length(reformulas::findbars(formula)) > 0L
}

.mc_response_formula <- function(formula) {
  formula <- stats::as.formula(formula)
  if (length(formula) != 2L) {
    stop(
      "formula must be one-sided because expression rows supply the response",
      call. = FALSE
    )
  }
  result <- stats::as.formula(
    call("~", as.name(".mixedcirc_response"), formula[[2L]]),
    env = environment(formula)
  )
  environment(result) <- environment(formula)
  result
}

.mc_expression_input <- function(exprObj, weightsMatrix, useWeights) {
  useWeights <- .mc_scalar_logical(useWeights, "useWeights")
  inherited_weights <- NULL

  if (methods::is(exprObj, "EList")) {
    expression <- as.matrix(exprObj$E)
    inherited_weights <- exprObj$weights
  } else if (methods::is(exprObj, "ExpressionSet")) {
    expression <- as.matrix(Biobase::exprs(exprObj))
  } else if (is.data.frame(exprObj)) {
    expression <- as.matrix(exprObj)
  } else if (is.matrix(exprObj) || methods::is(exprObj, "sparseMatrix")) {
    expression <- as.matrix(exprObj)
  } else {
    stop(
      "exprObj must be a matrix, data frame, sparse matrix, ExpressionSet, or EList",
      call. = FALSE
    )
  }

  if (!is.numeric(expression)) {
    stop("exprObj must contain numeric values", call. = FALSE)
  }
  if (!nrow(expression) || !ncol(expression)) {
    stop("exprObj must have at least one row and one column", call. = FALSE)
  }

  weights <- if (!is.null(weightsMatrix)) weightsMatrix else inherited_weights
  if (!useWeights) weights <- NULL

  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    if (length(dim(weights)) != 2L || any(dim(expression) != dim(weights))) {
      stop("exprObj and weightsMatrix must have the same dimensions", call. = FALSE)
    }
    if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("weights must contain finite positive values", call. = FALSE)
    }
    dimnames(weights) <- dimnames(expression)
  }

  list(expression = expression, weights = weights)
}

.mc_prepare_model_input <- function(exprObj,
                                    formula,
                                    data,
                                    weightsMatrix,
                                    useWeights,
                                    allow_response_na = FALSE) {
  input <- .mc_expression_input(exprObj, weightsMatrix, useWeights)
  expression <- input$expression
  weights <- input$weights
  formula <- stats::as.formula(formula)
  data <- as.data.frame(data)

  variables <- unique(all.vars(formula))
  missing_variables <- setdiff(variables, colnames(data))
  if (length(missing_variables)) {
    stop(
      "Variable in formula not found in data: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }
  if (ncol(expression) != nrow(data)) {
    stop(
      "the number of samples in exprObj (columns) must equal the rows in data",
      call. = FALSE
    )
  }

  complete_metadata <- if (length(variables)) {
    stats::complete.cases(data[, variables, drop = FALSE])
  } else {
    rep(TRUE, nrow(data))
  }
  if (any(!complete_metadata)) {
    warning(
      sum(!complete_metadata),
      " sample(s) with missing model metadata were dropped",
      call. = FALSE,
      immediate. = TRUE
    )
    data <- droplevels(data[complete_metadata, , drop = FALSE])
    expression <- expression[, complete_metadata, drop = FALSE]
    if (!is.null(weights)) weights <- weights[, complete_metadata, drop = FALSE]
  }

  if (nrow(data) < 2L) {
    stop("at least two complete samples are required", call. = FALSE)
  }
  if (!allow_response_na && any(!is.finite(expression))) {
    stop("exprObj contains NA, NaN, or infinite values", call. = FALSE)
  }
  if (allow_response_na && any(is.infinite(expression))) {
    stop("exprObj contains infinite values", call. = FALSE)
  }

  row_variances <- apply(expression, 1L, stats::var, na.rm = allow_response_na)
  invalid_variance <- which(!is.finite(row_variances) | row_variances <= 0)
  if (length(invalid_variance)) {
    feature <- rownames(expression)[invalid_variance[[1L]]]
    if (is.null(feature) || is.na(feature) || !nzchar(feature)) {
      feature <- invalid_variance[[1L]]
    }
    stop("response feature ", feature, " has zero or undefined variance", call. = FALSE)
  }

  if (
    !is.null(colnames(expression)) &&
      !is.null(rownames(data)) &&
      !identical(colnames(expression), rownames(data))
  ) {
    warning(
      "sample names in exprObj do not match row names in data; positional order is used",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  list(
    expression = expression,
    weights = weights,
    formula = formula,
    fit_formula = .mc_response_formula(formula),
    data = data
  )
}

.mc_weighted_var <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  if (length(x) != length(weights)) {
    stop("x and weights must have the same length", call. = FALSE)
  }
  keep <- is.finite(x) & is.finite(weights)
  x <- x[keep]
  weights <- weights[keep]
  if (length(x) < 2L || any(weights < 0) || sum(weights) <= 0) return(NA_real_)
  weights <- weights * length(weights) / sum(weights)
  sum(weights * (x - stats::weighted.mean(x, weights))^2) /
    (sum(weights) - 1)
}

.mc_collinearity_score <- function(fit) {
  covariance <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
  if (is.null(covariance) || !length(covariance) || any(!is.finite(covariance))) {
    return(1)
  }
  correlation <- stats::cov2cor(covariance)
  if (nrow(correlation) < 2L) return(0)
  max(abs(correlation[lower.tri(correlation)]))
}

.mc_check_model <- function(fit,
                            showWarnings = TRUE,
                            colinearityCutoff = 0.999,
                            checkIntercept = TRUE,
                            collinearity_score = .mc_collinearity_score) {
  if (!isTRUE(showWarnings)) return(invisible(fit))

  model_terms <- tryCatch(stats::terms(fit), error = function(error) NULL)
  if (isTRUE(checkIntercept) &&
      !is.null(model_terms) &&
      identical(attr(model_terms, "intercept"), 0L)) {
    stop("No intercept term was specified in the formula", call. = FALSE)
  }

  coefficients <- tryCatch(stats::coef(fit), error = function(e) numeric())
  if (anyNA(coefficients)) {
    stop("the fitted design matrix is rank deficient", call. = FALSE)
  }
  if (inherits(fit, "lm")) {
    score <- collinearity_score(fit)
    if (is.finite(score) && score > colinearityCutoff) {
      stop(
        "model collinearity score ", format(score, digits = 4L),
        " exceeds ", colinearityCutoff,
        call. = FALSE
      )
    }
  } else if (inherits(fit, "merMod")) {
    covariance <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
    if (!is.null(covariance) && nrow(covariance) > 1L) {
      condition_number <- tryCatch(
        kappa(stats::cov2cor(covariance)),
        error = function(e) NA_real_
      )
      if (is.finite(condition_number) && condition_number > 1e8) {
        stop(
          "fixed-effect covariance condition number is too high: ",
          format(condition_number, digits = 3L),
          call. = FALSE
        )
      }
    }
  }
  invisible(fit)
}

.mc_fit_one_model <- function(response,
                              weights,
                              fit_formula,
                              data,
                              mixed,
                              REML,
                              control,
                              showWarnings,
                              colinearityCutoff,
                              fxn,
                              check_model,
                              collinearity_score,
                              dots) {
  model_data <- data
  model_data$.mixedcirc_response <- as.numeric(response)
  args <- c(
    list(
      formula = fit_formula,
      data = model_data,
      na.action = stats::na.exclude
    ),
    dots
  )
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    weights <- weights / mean(weights, na.rm = TRUE)
    args$weights <- weights
  }
  if (mixed) {
    args$REML <- REML
    args$control <- control
    fit <- if (showWarnings) {
      do.call(lme4::lmer, args)
    } else {
      suppressWarnings(do.call(lme4::lmer, args))
    }
  } else {
    fit <- if (showWarnings) {
      do.call(stats::lm, args)
    } else {
      suppressWarnings(do.call(stats::lm, args))
    }
  }
  check_model(
    fit,
    showWarnings = showWarnings,
    colinearityCutoff = colinearityCutoff,
    collinearity_score = collinearity_score
  )
  fxn(fit)
}

.mc_model_chunks <- function(number_of_features, BPPARAM) {
  if (!methods::is(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must be a BiocParallelParam object", call. = FALSE)
  }
  workers <- tryCatch(BiocParallel::bpworkers(BPPARAM), error = function(e) 1L)
  workers <- max(1L, as.integer(workers[[1L]]))
  number_of_chunks <- min(number_of_features, max(1L, workers * 4L))
  parallel::splitIndices(number_of_features, number_of_chunks)
}

.mc_fit_models <- function(exprObj,
                           formula,
                           data,
                           REML,
                           useWeights,
                           weightsMatrix,
                           showWarnings,
                           fxn,
                           control,
                           quiet,
                           BPPARAM,
                           colinearityCutoff,
                           allow_response_na,
                           dots) {
  REML <- .mc_scalar_logical(REML, "REML")
  showWarnings <- .mc_scalar_logical(showWarnings, "showWarnings")
  quiet <- .mc_scalar_logical(quiet, "quiet")
  if (!is.function(fxn)) stop("fxn must be a function", call. = FALSE)
  if (is.null(control)) control <- .mixedcirc_vp_control

  prepared <- .mc_prepare_model_input(
    exprObj = exprObj,
    formula = formula,
    data = data,
    weightsMatrix = weightsMatrix,
    useWeights = useWeights,
    allow_response_na = allow_response_na
  )
  mixed <- .mc_is_mixed_formula(prepared$formula)
  chunks <- .mc_model_chunks(nrow(prepared$expression), BPPARAM)
  if (!quiet) {
    message(
      "Fitting ", nrow(prepared$expression), " feature(s) in ",
      length(chunks), " chunk(s)"
    )
  }

  fit_chunk <- function(indices,
                        expression,
                        weights,
                        fit_formula,
                        model_data,
                        mixed,
                        REML,
                        control,
                        showWarnings,
                        colinearityCutoff,
                        fxn,
                        fit_one_model,
                        check_model,
                        collinearity_score,
                        dots) {
    lapply(indices, function(index) {
      feature_name <- rownames(expression)[index]
      if (is.null(feature_name) || is.na(feature_name) || !nzchar(feature_name)) {
        feature_name <- as.character(index)
      }
      tryCatch(
        fit_one_model(
          response = expression[index, ],
          weights = if (is.null(weights)) NULL else weights[index, ],
          fit_formula = fit_formula,
          data = model_data,
          mixed = mixed,
          REML = REML,
          control = control,
          showWarnings = showWarnings,
          colinearityCutoff = colinearityCutoff,
          fxn = fxn,
          check_model = check_model,
          collinearity_score = collinearity_score,
          dots = dots
        ),
        error = function(error) {
          stop("feature ", feature_name, ": ", conditionMessage(error), call. = FALSE)
        }
      )
    })
  }

  chunk_results <- BiocParallel::bplapply(
    chunks,
    fit_chunk,
    expression = prepared$expression,
    weights = prepared$weights,
    fit_formula = prepared$fit_formula,
    model_data = prepared$data,
    mixed = mixed,
    REML = REML,
    control = control,
    showWarnings = showWarnings,
    colinearityCutoff = colinearityCutoff,
    fxn = fxn,
    fit_one_model = .mc_fit_one_model,
    check_model = .mc_check_model,
    collinearity_score = .mc_collinearity_score,
    dots = dots,
    BPPARAM = BPPARAM
  )
  results <- unlist(chunk_results, recursive = FALSE, use.names = FALSE)
  names(results) <- rownames(prepared$expression)

  methods::new(
    "VarParFitList",
    results,
    method = if (mixed) "lmer" else "lm"
  )
}

#' Fit one linear or mixed model per feature
#'
#' Fits a separate model to every row of a feature-by-sample matrix. The
#' implementation is local to mixedcirc and supports ordinary and mixed models,
#' observation-level weights, and `BiocParallel` backends.
#'
#' @param exprObj A numeric feature-by-sample matrix, data frame, sparse matrix,
#'   `Biobase::ExpressionSet`, or `limma::EList`.
#' @param formula A one-sided model formula, for example
#'   `~ age + (1 | subject)`.
#' @param data Sample metadata with one row per column of `exprObj`.
#' @param REML Logical; use restricted maximum likelihood for mixed models.
#' @param useWeights Logical; use weights supplied by an `EList` or
#'   `weightsMatrix`.
#' @param weightsMatrix Optional positive feature-by-sample weight matrix.
#' @param showWarnings Logical; retain fitting warnings and run model checks.
#' @param fxn Function applied to each fitted model. The identity function
#'   returns the models themselves.
#' @param control Optional `lme4::lmerControl()` object.
#' @param quiet Logical; suppress progress messages.
#' @param BPPARAM A `BiocParallelParam` object.
#' @param ... Additional arguments passed to `lm()` or `lmer()`.
#'
#' @return A `VarParFitList` with one entry per feature.
#' @aliases fitVarPartModel,ANY-method
#'
#' @examples
#' set.seed(1)
#' metadata <- data.frame(subject = factor(rep(1:5, each = 3)))
#' expression <- matrix(rnorm(3 * nrow(metadata)), nrow = 3)
#' fits <- fitVarPartModel(
#'   expression, ~ (1 | subject), metadata,
#'   BPPARAM = BiocParallel::SerialParam(), quiet = TRUE
#' )
#' length(fits)
#'
#' @export
fitVarPartModel <- function(exprObj,
                            formula,
                            data,
                            REML = FALSE,
                            useWeights = TRUE,
                            weightsMatrix = NULL,
                            showWarnings = TRUE,
                            fxn = identity,
                            control = NULL,
                            quiet = FALSE,
                            BPPARAM = BiocParallel::SerialParam(),
                            ...) {
  global_weights <- .mc_capture_global_weights()
  on.exit(.mc_restore_global_weights(global_weights), add = TRUE)
  .mc_fit_models(
    exprObj = exprObj,
    formula = formula,
    data = data,
    REML = REML,
    useWeights = useWeights,
    weightsMatrix = weightsMatrix,
    showWarnings = showWarnings,
    fxn = fxn,
    control = control,
    quiet = quiet,
    BPPARAM = BPPARAM,
    colinearityCutoff = 0.999,
    allow_response_na = FALSE,
    dots = list(...)
  )
}

# Internal fitting helper used by the chunked RRBS workflow.
variance_fit <- function(exprObj,
                         formula,
                         data,
                         REML = FALSE,
                         useWeights = TRUE,
                         weightsMatrix = NULL,
                         showWarnings = TRUE,
                         fxn = identity,
                         colinearityCutoff = 0.999,
                         control = NULL,
                         quiet = FALSE,
                         BPPARAM = BiocParallel::SerialParam(),
                         ignore_na = FALSE,
                         ...) {
  global_weights <- .mc_capture_global_weights()
  on.exit(.mc_restore_global_weights(global_weights), add = TRUE)
  ignore_na <- .mc_scalar_logical(ignore_na, "ignore_na")
  .mc_fit_models(
    exprObj = exprObj,
    formula = formula,
    data = data,
    REML = REML,
    useWeights = useWeights,
    weightsMatrix = weightsMatrix,
    showWarnings = showWarnings,
    fxn = fxn,
    control = control,
    quiet = quiet,
    BPPARAM = BPPARAM,
    colinearityCutoff = colinearityCutoff,
    allow_response_na = ignore_na,
    dots = list(...)
  )
}

calcVarPart <- function(fit, returnFractions = TRUE, showWarnings = TRUE, ...) {
  returnFractions <- .mc_scalar_logical(returnFractions, "returnFractions")
  .mc_check_model(fit, showWarnings = showWarnings)

  if (inherits(fit, "lm")) {
    design <- stats::model.matrix(fit)
    term_labels <- attr(stats::terms(fit), "term.labels")
    assignment <- fit$assign
    coefficients <- stats::coef(fit)
    weights <- stats::weights(fit)
    if (is.null(weights)) weights <- rep(1, nrow(design))

    sums_of_squares <- numeric()
    if (length(term_labels)) {
      sums_of_squares <- vapply(seq_along(term_labels), function(term_index) {
        columns <- which(assignment == term_index)
        contribution <- if (length(columns)) {
          as.vector(design[, columns, drop = FALSE] %*% coefficients[columns])
        } else {
          rep(0, nrow(design))
        }
        .mc_weighted_var(contribution, weights) * (length(contribution) - 1)
      }, numeric(1))
      names(sums_of_squares) <- term_labels
    }
    sums_of_squares <- c(
      sums_of_squares,
      Residuals = stats::sigma(fit)^2 * stats::df.residual(fit)
    )
    result <- if (returnFractions) {
      sums_of_squares / sum(sums_of_squares)
    } else {
      sums_of_squares / nrow(design)
    }
    return(result)
  }

  if (!inherits(fit, "lmerMod")) {
    stop("calcVarPart supports lm and lmerMod fits", call. = FALSE)
  }

  variance_correlation <- lme4::VarCorr(fit)
  grouping_names <- sort(names(variance_correlation))
  components <- unlist(lapply(grouping_names, function(grouping) {
    standard_deviations <- attr(variance_correlation[[grouping]], "stddev")
    values <- standard_deviations^2
    coefficient_names <- names(standard_deviations)
    names(values) <- ifelse(
      coefficient_names == "(Intercept)",
      grouping,
      paste(grouping, coefficient_names, sep = ".")
    )
    values
  }), use.names = TRUE)

  design <- lme4::getME(fit, "X")
  fixed_coefficients <- lme4::fixef(fit)
  fixed_columns <- which(colnames(design) != "(Intercept)")
  if (length(fixed_columns)) {
    fixed_contributions <- sweep(
      design[, fixed_columns, drop = FALSE],
      2L,
      fixed_coefficients[fixed_columns],
      `*`
    )
    weights <- stats::weights(fit)
    if (is.null(weights)) weights <- rep(1, nrow(design))
    marginal_variances <- apply(fixed_contributions, 2L, .mc_weighted_var, weights = weights)
    scaling <- (nrow(design) - 1) / nrow(design)
    marginal_variances <- marginal_variances * scaling
    total_fixed_variance <- .mc_weighted_var(
      rowSums(fixed_contributions), weights
    ) * scaling
    if (sum(marginal_variances) > 0 && is.finite(total_fixed_variance)) {
      fixed_components <- marginal_variances / sum(marginal_variances) *
        total_fixed_variance
    } else {
      fixed_components <- setNames(rep(0, length(marginal_variances)),
                                   names(marginal_variances))
    }
    components <- c(components, fixed_components)
  }
  components <- c(components, Residuals = stats::sigma(fit)^2)
  if (returnFractions) components <- components / sum(components)
  components
}

extractVarPart <- function(modelList, ...) {
  if (!length(modelList)) stop("modelList is empty", call. = FALSE)
  fractions <- lapply(modelList, calcVarPart, ...)
  component_names <- unique(unlist(lapply(fractions, names), use.names = FALSE))
  result <- matrix(
    0,
    nrow = length(fractions),
    ncol = length(component_names),
    dimnames = list(names(modelList), component_names)
  )
  for (index in seq_along(fractions)) {
    result[index, names(fractions[[index]])] <- fractions[[index]]
  }
  result <- as.data.frame(result, check.names = FALSE)
  methods::new(
    "varPartResults",
    result,
    type = if (inherits(modelList[[1L]], "lm")) "anova" else "linear mixed model",
    method = "Variance explained (fraction)"
  )
}

#' Estimate variance fractions for every feature
#'
#' Fits the requested model to each feature and reports the fraction attributed
#' to each fixed term, random-effect component, and residual variation.
#' For mixed models with several fixed-effect columns, their combined fitted
#' variance is allocated in proportion to the marginal variance of each
#' column's contribution. The result is a descriptive package-local
#' decomposition rather than a call to the `variancePartition` package.
#'
#' @inheritParams fitVarPartModel
#'
#' @return A data-frame-like `varPartResults` object. Every row sums to one,
#'   subject to numerical tolerance.
#' @aliases fitExtractVarPartModel,ANY-method
#'
#' @examples
#' set.seed(2)
#' metadata <- data.frame(
#'   subject = factor(rep(1:6, each = 3)),
#'   x = rnorm(18)
#' )
#' subject_effect <- rnorm(6)[metadata$subject]
#' expression <- rbind(
#'   feature_1 = metadata$x + subject_effect + rnorm(18, sd = 0.5),
#'   feature_2 = subject_effect + rnorm(18)
#' )
#' fractions <- fitExtractVarPartModel(
#'   expression, ~ x + (1 | subject), metadata,
#'   BPPARAM = BiocParallel::SerialParam(), quiet = TRUE
#' )
#' rowSums(as.data.frame(fractions))
#'
#' @export
fitExtractVarPartModel <- function(exprObj,
                                   formula,
                                   data,
                                   REML = FALSE,
                                   useWeights = TRUE,
                                   weightsMatrix = NULL,
                                   showWarnings = TRUE,
                                   control = NULL,
                                   quiet = FALSE,
                                   BPPARAM = BiocParallel::SerialParam(),
                                   ...) {
  models <- fitVarPartModel(
    exprObj = exprObj,
    formula = formula,
    data = data,
    REML = REML,
    useWeights = useWeights,
    weightsMatrix = weightsMatrix,
    showWarnings = showWarnings,
    fxn = identity,
    control = control,
    quiet = quiet,
    BPPARAM = BPPARAM,
    ...
  )
  extractVarPart(models, showWarnings = FALSE)
}

.mc_prepare_count_input <- function(counts, data, formula, lib.size, weights) {
  out <- list()
  if (methods::is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(lib.size)) {
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    }
    count_matrix <- as.matrix(counts$counts)
  } else if (methods::is(counts, "ExpressionSet")) {
    if (length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
    if (length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
    count_matrix <- as.matrix(Biobase::exprs(counts))
  } else if (!is.null(dim(counts))) {
    count_matrix <- as.matrix(counts)
  } else {
    stop("counts must be a count matrix, DGEList, or ExpressionSet", call. = FALSE)
  }

  if (!is.numeric(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
    stop("counts must contain finite non-negative numeric values", call. = FALSE)
  }
  if (nrow(count_matrix) < 2L) {
    stop("at least two features are required to fit a mean-variance trend", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (ncol(count_matrix) != nrow(data)) {
    stop("the columns of counts must equal the rows in data", call. = FALSE)
  }

  formula <- stats::as.formula(formula)
  variables <- unique(all.vars(formula))
  missing_variables <- setdiff(variables, colnames(data))
  if (length(missing_variables)) {
    stop(
      "Variable in formula not found in data: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }
  complete_metadata <- if (length(variables)) {
    stats::complete.cases(data[, variables, drop = FALSE])
  } else {
    rep(TRUE, nrow(data))
  }
  if (any(!complete_metadata)) {
    warning(
      sum(!complete_metadata),
      " sample(s) with missing model metadata were dropped",
      call. = FALSE
    )
    count_matrix <- count_matrix[, complete_metadata, drop = FALSE]
    data <- droplevels(data[complete_metadata, , drop = FALSE])
    if (!is.null(lib.size)) lib.size <- lib.size[complete_metadata]
    if (!is.null(weights)) {
      weights <- if (is.matrix(weights)) {
        weights[, complete_metadata, drop = FALSE]
      } else {
        weights[complete_metadata]
      }
    }
    if (!is.null(out$targets) && nrow(out$targets) == length(complete_metadata)) {
      out$targets <- out$targets[complete_metadata, , drop = FALSE]
    }
  }
  if (is.null(lib.size)) lib.size <- colSums(count_matrix)
  if (
    !is.numeric(lib.size) || length(lib.size) != ncol(count_matrix) ||
      any(!is.finite(lib.size)) || any(lib.size < 0)
  ) {
    stop("lib.size must contain one finite non-negative value per sample", call. = FALSE)
  }
  if (any(lib.size == 0)) {
    stop("library sizes must be positive", call. = FALSE)
  }

  list(
    counts = count_matrix,
    data = data,
    formula = formula,
    lib.size = lib.size,
    weights = weights,
    out = out
  )
}

.mc_prior_counts <- function(counts, lib.size, prior.count, scaledByLib) {
  scaledByLib <- .mc_scalar_logical(scaledByLib, "scaledByLib")
  if (!is.numeric(prior.count) || length(prior.count) != 1L ||
      !is.finite(prior.count) || prior.count < 0) {
    stop("prior.count must be one finite non-negative number", call. = FALSE)
  }
  scale <- if (scaledByLib) lib.size / mean(lib.size) else rep(1, length(lib.size))
  counts + matrix(
    prior.count * scale,
    nrow = nrow(counts),
    ncol = ncol(counts),
    byrow = TRUE
  )
}

.mc_voom_trend <- function(x, y, span) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) return(NULL)

  if (is.character(span)) {
    if (length(span) != 1L || span != "auto") {
      stop("span must be numeric or 'auto'", call. = FALSE)
    }
    if (length(unique(x)) >= 4L) {
      smooth <- tryCatch(
        stats::smooth.spline(x, y, cv = TRUE),
        error = function(e) NULL
      )
      if (!is.null(smooth)) {
        prediction <- stats::predict(smooth, sort(unique(x)))
        return(list(x = prediction$x, y = prediction$y))
      }
    }
    span <- 0.5
  }
  if (!is.numeric(span) || length(span) != 1L || !is.finite(span) ||
      span <= 0 || span > 1) {
    stop("span must be in (0, 1] or equal to 'auto'", call. = FALSE)
  }
  stats::lowess(x, y, f = span)
}

#' Transform RNA-seq counts for mixed-model analysis
#'
#' Converts counts to log2 counts-per-million and estimates observation-level
#' precision weights from the fitted mean-variance trend. Unlike `limma::voom`,
#' the model formula may contain random effects.
#'
#' @param counts A non-negative feature-by-sample count matrix,
#'   `Biobase::ExpressionSet`, or `edgeR::DGEList`.
#' @param formula A one-sided fixed- or mixed-effects formula.
#' @param data Sample metadata with one row per count-matrix column.
#' @param lib.size Optional effective library sizes.
#' @param normalize.method Method passed to `limma::normalizeBetweenArrays()`.
#' @param span LOWESS span in `(0, 1]`, or `"auto"`.
#' @param weights Optional feature-by-sample or sample-level prior weights.
#' @param prior.count Prior count used for the log transformation.
#' @param prior.count.for.weights Prior count used when deriving weights from
#'   augmented counts.
#' @param plot Logical; draw the fitted mean-variance trend.
#' @param save.plot Logical; retain trend coordinates in the output.
#' @param rescaleWeightsAfter Logical; multiply the estimated precision weights
#'   by the normalized prior weights.
#' @param scaledByLib Logical; scale the prior count by library size.
#' @param priorWeightsAsCounts Logical; derive prior weights from augmented
#'   counts when no weights are supplied.
#' @param quiet Logical; suppress model-fitting messages.
#' @param BPPARAM A `BiocParallelParam` object.
#' @param ... Additional fitting arguments.
#'
#' @return A `limma::EList` containing logCPM values and precision weights.
#'
#' @examples
#' set.seed(3)
#' metadata <- data.frame(subject = factor(rep(1:4, each = 3)))
#' counts <- matrix(rpois(20 * nrow(metadata), lambda = 30), nrow = 20)
#' vobj <- voomWithDreamWeights(
#'   counts, ~ (1 | subject), metadata,
#'   BPPARAM = BiocParallel::SerialParam(), quiet = TRUE
#' )
#' dim(vobj$weights)
#'
#' @seealso [limma::voom()]
#' @export
voomWithDreamWeights <- function(counts,
                                 formula,
                                 data,
                                 lib.size = NULL,
                                 normalize.method = "none",
                                 span = 0.5,
                                 weights = NULL,
                                 prior.count = 0.5,
                                 prior.count.for.weights = prior.count,
                                 plot = FALSE,
                                 save.plot = FALSE,
                                 rescaleWeightsAfter = FALSE,
                                 scaledByLib = FALSE,
                                 priorWeightsAsCounts = FALSE,
                                 quiet = FALSE,
                                 BPPARAM = BiocParallel::SerialParam(),
                                 ...) {
  global_weights <- .mc_capture_global_weights()
  on.exit(.mc_restore_global_weights(global_weights), add = TRUE)
  quiet <- .mc_scalar_logical(quiet, "quiet")
  plot <- .mc_scalar_logical(plot, "plot")
  save.plot <- .mc_scalar_logical(save.plot, "save.plot")
  rescaleWeightsAfter <- .mc_scalar_logical(
    rescaleWeightsAfter, "rescaleWeightsAfter"
  )
  priorWeightsAsCounts <- .mc_scalar_logical(
    priorWeightsAsCounts, "priorWeightsAsCounts"
  )

  prepared <- .mc_prepare_count_input(counts, data, formula, lib.size, weights)
  counts <- prepared$counts
  data <- prepared$data
  formula <- prepared$formula
  lib.size <- prepared$lib.size
  weights <- prepared$weights
  out <- prepared$out

  augmented <- .mc_prior_counts(counts, lib.size, prior.count, scaledByLib)
  if (priorWeightsAsCounts && (is.null(weights) || all(weights == 1))) {
    if (!is.numeric(prior.count.for.weights) ||
        length(prior.count.for.weights) != 1L ||
        !is.finite(prior.count.for.weights) || prior.count.for.weights < 0) {
      stop("prior.count.for.weights must be finite and non-negative", call. = FALSE)
    }
    weights <- augmented + prior.count.for.weights
  }

  log_cpm <- t(log2(t(augmented) / (lib.size + 1) * 1e6))
  log_cpm <- limma::normalizeBetweenArrays(log_cpm, method = normalize.method)
  if (is.null(weights)) weights <- rep(1, ncol(log_cpm))
  if (is.matrix(weights)) {
    if (length(dim(weights)) != 2L || any(dim(weights) != dim(log_cpm))) {
      stop("matrix weights must have the same dimensions as counts", call. = FALSE)
    }
    weights_matrix <- weights
  } else {
    if (length(weights) != ncol(log_cpm)) {
      stop("sample weights must have length ncol(counts)", call. = FALSE)
    }
    weights_matrix <- limma::asMatrixWeights(weights, dim(log_cpm))
  }
  if (any(!is.finite(weights_matrix)) || any(weights_matrix <= 0)) {
    stop("weights must contain finite positive values", call. = FALSE)
  }
  weights_matrix <- weights_matrix / rowMeans(weights_matrix)
  dimnames(weights_matrix) <- dimnames(log_cpm)

  design <- NULL
  if (.mc_is_mixed_formula(formula)) {
    fitted_summary <- fitVarPartModel(
      log_cpm,
      formula,
      data,
      weightsMatrix = weights_matrix,
      showWarnings = FALSE,
      fxn = function(fit) {
        list(
          sigma = stats::sigma(fit),
          average = mean(stats::model.response(stats::model.frame(fit)), na.rm = TRUE),
          fitted = stats::fitted(fit)
        )
      },
      quiet = quiet,
      BPPARAM = BPPARAM,
      ...
    )
    sigma_values <- vapply(fitted_summary, `[[`, numeric(1), "sigma")
    average_values <- vapply(fitted_summary, `[[`, numeric(1), "average")
    fitted_values <- do.call(rbind, lapply(fitted_summary, `[[`, "fitted"))
    rownames(fitted_values) <- rownames(log_cpm)
    residual_df <- rep(2, length(sigma_values))
  } else {
    design <- stats::model.matrix(formula, data)
    fit <- limma::lmFit(log_cpm, design, weights = weights_matrix, ...)
    sigma_values <- fit$sigma
    average_values <- if (is.null(fit$Amean)) rowMeans(log_cpm) else fit$Amean
    residual_df <- fit$df.residual
    if (fit$rank < ncol(design)) {
      estimable <- fit$pivot[seq_len(fit$rank)]
      fitted_values <- fit$coefficients[, estimable, drop = FALSE] %*%
        t(fit$design[, estimable, drop = FALSE])
    } else {
      fitted_values <- fit$coefficients %*% t(fit$design)
    }
  }

  replicated <- sum(residual_df > 0)
  if (replicated < 2L) {
    warning("insufficient replication to estimate a mean-variance trend; using unit weights")
    out$E <- log_cpm
    out$weights <- matrix(1, nrow(log_cpm), ncol(log_cpm), dimnames = dimnames(log_cpm))
    if (!is.null(design)) out$design <- design
    out$targets <- if (is.null(out$targets)) {
      data.frame(lib.size = lib.size)
    } else {
      transform(out$targets, lib.size = lib.size)
    }
    return(methods::new("EList", out))
  }

  trend_x <- average_values + mean(log2(lib.size + 1)) - log2(1e6)
  trend_y <- sqrt(pmax(sigma_values, 0))
  nonzero <- rowSums(counts) > 0
  trend <- .mc_voom_trend(trend_x[nonzero], trend_y[nonzero], span)
  if (is.null(trend)) {
    warning("too few nonzero features for a mean-variance trend; using unit weights")
    precision <- matrix(1, nrow(log_cpm), ncol(log_cpm), dimnames = dimnames(log_cpm))
  } else {
    trend_function <- stats::approxfun(trend, rule = 2)
    fitted_cpm <- 2^fitted_values
    fitted_count <- 1e-6 * t(t(fitted_cpm) * (lib.size + 1))
    fitted_log_count <- log2(pmax(fitted_count, .Machine$double.xmin))
    predicted_root_sd <- trend_function(fitted_log_count)
    predicted_root_sd <- pmax(predicted_root_sd, sqrt(.Machine$double.eps))
    precision <- 1 / predicted_root_sd^4
    dim(precision) <- dim(fitted_log_count)
    dimnames(precision) <- dimnames(log_cpm)
  }

  out$E <- log_cpm
  out$weights <- precision
  if (!is.null(design)) out$design <- design
  out$targets <- if (is.null(out$targets)) {
    data.frame(lib.size = lib.size)
  } else {
    transform(out$targets, lib.size = lib.size)
  }
  if (rescaleWeightsAfter) {
    out$weights <- out$weights * weights_matrix
    out$targets$sample.weights <- colMeans(weights_matrix)
  }
  if (save.plot && !is.null(trend)) {
    out$voom.xy <- list(
      x = trend_x[nonzero],
      y = trend_y[nonzero],
      xlab = "log2(count size + 0.5)",
      ylab = "sqrt(standard deviation)"
    )
    out$voom.line <- trend
  }
  if (plot && !is.null(trend)) {
    graphics::plot(
      trend_x[nonzero], trend_y[nonzero],
      xlab = "log2(count size + 0.5)",
      ylab = "sqrt(standard deviation)",
      pch = 16, cex = 0.25
    )
    graphics::title("voom: mean-variance trend")
    graphics::lines(trend, col = "red")
  }
  if (max(out$weights) > 1e8) {
    warning(
      "very large precision weights were estimated; inspect the mean-variance trend",
      call. = FALSE
    )
  }
  methods::new("EList", out)
}
