#' Experimental integration of circadian datasets using SGCCA
#'
#' Performs multi-block data integration using sparse generalized canonical
#' correlation analysis (SGCCA). The biological datasets are
#' integrated together with one or two design blocks representing circadian
#' rhythm and, optionally, group structure.
#'
#' @param data_input A named list containing at least two numeric matrices or
#'   data frames. Each element
#'   is one dataset, with rows corresponding to samples and columns to variables.
#'   All datasets must have the same number of rows.
#' @param meta_input A data.frame or matrix containing sample-level metadata.
#'   Must include the columns referenced by `time`, `group`, and `id` when used.
#' @param time A character string giving the name of the column in `meta_input`
#'   that contains time.
#' @param group Optional character string giving the name of the grouping column
#'   in `meta_input`. If supplied, a group-design block is added.
#' @param id Optional character string giving the name of the subject identifier
#'   column in `meta_input`. Used only when `decompose = TRUE`.
#' @param period Numeric scalar specifying the assumed period. Default `24`.
#' @param ncomp Number of SGCCA components to estimate. Must be a positive integer.
#'   Default `2`.
#' @param variables Optional named list controlling the number of variables kept
#'   per component in each dataset. Each supplied element must be an integer
#'   vector of length `ncomp`; its name must identify a block in `data_input`.
#'   Blocks omitted from this list keep all variables.
#' @param decompose Logical. If `TRUE`, within-subject variation is extracted for
#'   each dataset before integration.
#'   Requires `id`. Default `FALSE`.
#' @param center Logical. If `TRUE`, columns of each dataset are centered.
#'   Default `TRUE`.
#' @param scale Logical. If `TRUE`, columns of each dataset are scaled to unit
#'   variance. Default `FALSE`.
#' @param merge Logical. If `TRUE`, group and rhythm design variables are merged
#'   into a single design block. Default `FALSE`.
#' @param no_correlation Logical. If `TRUE` and group information is provided in
#'   separate blocks, the corresponding entries of the SGCCA design matrix are
#'   set to zero, so the group and rhythm blocks are not directly coupled.
#'   Default `FALSE`.
#' @param no_interaction Logical. If `TRUE`, no group-specific rhythm interaction
#'   is encoded. In grouped analysis, rhythm is then represented by common
#'   `inphase` and `outphase` terms and group is represented separately.
#'   Default `FALSE`.
#' @param max.iter Maximum number of SGCCA iterations. Default `100000`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default `FALSE`.
#'
#' @return
#' An object of class `mixedcirc_integration` containing:
#' \describe{
#'   \item{partial}{A named list of partial block scores, one matrix per dataset.}
#'   \item{average}{A matrix of average sample scores across biological datasets.}
#'   \item{loadings}{A named list of loading matrices, one per dataset.}
#'   \item{model}{The fitted model object.}
#' }
#'
#' @details
#' The function augments the biological datasets with design information derived
#' from the sample metadata:
#'
#' \deqn{inphase = \cos(2\pi t / period)}
#' \deqn{outphase = \sin(2\pi t / period)}
#'
#' If group information is provided and `no_interaction = FALSE`, the rhythm
#' block uses group-specific circadian terms. If `no_interaction = TRUE`, the
#' rhythm block is shared across groups and the group structure is represented
#' separately.
#'
#' @section Experimental status:
#' This function is experimental. Because the SGCCA model is supervised by the
#' time and optional group variables, its scores are selected to follow those
#' variables. Testing the learned scores against the same time or group variables
#' in the same samples does not provide independent confirmatory evidence. Use the
#' scores for exploratory analysis, or validate them on independent samples or
#' with a resampling procedure that repeats the complete integration workflow.
#'
#' @examples
#' library(mixedcirc)
#' data("circa_data")
#'
#' data_input <- list(
#'   a = circa_data$data_matrix[, 1:3],
#'   b = circa_data$data_matrix[, 4:7]
#' )
#'
#' results <- mixedcirc_integrate(
#'   data_input = data_input,
#'   meta_input = data.frame(
#'     time = circa_data$time,
#'     group = circa_data$group,
#'     id = circa_data$id
#'   ),
#'   time = "time",
#'   group = "group",
#'   id = "id"
#' )
#'
#' score_matrix <- mixedcirc_getscore(results, type = "average")
#'
#' @export
mixedcirc_integrate <- function(data_input = NULL,
                                meta_input = NULL,
                                time = NULL,
                                group = NULL,
                                id = NULL,
                                period = 24,
                                ncomp = 2,
                                variables = NULL,
                                decompose = FALSE,
                                center = TRUE,
                                scale = FALSE,
                                merge = FALSE,
                                no_correlation = FALSE,
                                no_interaction = FALSE,
                                max.iter = 100000,
                                verbose = FALSE) {

  if (isTRUE(verbose)) cat("Checking inputs ...\n")

  logical_args <- list(
    decompose = decompose,
    center = center,
    scale = scale,
    merge = merge,
    no_correlation = no_correlation,
    no_interaction = no_interaction,
    verbose = verbose
  )
  bad_logical <- !vapply(
    logical_args,
    function(z) is.logical(z) && length(z) == 1L && !is.na(z),
    logical(1)
  )
  if (any(bad_logical)) {
    stop(names(logical_args)[which(bad_logical)[1L]], " must be TRUE or FALSE")
  }

  if (is.null(data_input) || !is.list(data_input) || length(data_input) < 2L) {
    stop("data_input must be a named list containing at least two biological data blocks")
  }

  if (is.null(names(data_input)) || anyNA(names(data_input)) ||
      any(names(data_input) == "") || anyDuplicated(names(data_input))) {
    stop("data_input must have unique, non-empty block names")
  }

  ok_blocks <- vapply(data_input, function(x) is.matrix(x) || is.data.frame(x), logical(1))
  if (any(!ok_blocks)) {
    stop(
      "All elements of data_input must be matrices or data frames. Problematic blocks: ",
      paste(names(data_input)[!ok_blocks], collapse = ", ")
    )
  }

  non_numeric_blocks <- !vapply(data_input, function(x) {
    if (is.data.frame(x)) all(vapply(x, is.numeric, logical(1))) else is.numeric(x)
  }, logical(1))
  if (any(non_numeric_blocks)) {
    stop(
      "All data blocks must contain only numeric columns. Problematic blocks: ",
      paste(names(data_input)[non_numeric_blocks], collapse = ", ")
    )
  }

  data_input <- Map(function(x, nm) {
    x <- as.matrix(x)
    if (ncol(x) == 0L) stop("Data block '", nm, "' has no variables")
    if (any(!is.finite(x))) {
      stop("Data block '", nm, "' contains missing or non-finite values")
    }
    if (is.null(colnames(x))) colnames(x) <- paste0("feature", seq_len(ncol(x)))
    if (anyDuplicated(colnames(x))) {
      warning("Duplicated variable names in data block '", nm, "' were made unique")
      colnames(x) <- make.unique(colnames(x))
    }
    x
  }, data_input, names(data_input))

  nrows_unique <- unique(vapply(data_input, nrow, integer(1)))
  if (length(nrows_unique) != 1) {
    stop("All elements of data_input must have the same number of rows")
  }
  n_samples <- nrows_unique

  if (is.null(meta_input) || !(is.data.frame(meta_input) || is.matrix(meta_input))) {
    stop("meta_input must be a data.frame or matrix")
  }
  meta_input <- as.data.frame(meta_input)

  if (nrow(meta_input) != n_samples) {
    stop("meta_input must have the same number of rows as each dataset in data_input")
  }

  if (!is.character(time) || length(time) != 1L || is.na(time) ||
      time == "" || !time %in% colnames(meta_input)) {
    stop("time must be the name of a column in meta_input")
  }
  if (!is.numeric(meta_input[[time]]) || any(!is.finite(meta_input[[time]]))) {
    stop("the time column must contain only finite numeric values")
  }

  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0) {
    stop("period must be a single positive finite numeric value")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1L || !is.finite(ncomp) ||
      ncomp <= 0 || ncomp %% 1 != 0) {
    stop("ncomp must be a positive integer")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp > min(vapply(data_input, ncol, integer(1)))) {
    stop("ncomp cannot exceed the number of variables in the smallest biological block")
  }
  if (!is.numeric(max.iter) || length(max.iter) != 1L || !is.finite(max.iter) ||
      max.iter < 1 || max.iter %% 1 != 0) {
    stop("max.iter must be a positive integer")
  }
  max.iter <- as.integer(max.iter)

  multiple_groups <- FALSE
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || is.na(group) ||
        group == "" || !group %in% colnames(meta_input)) {
      stop("group must be the name of a column in meta_input")
    }
    if (anyNA(meta_input[[group]])) {
      stop("the group column must not contain missing values")
    }
    multiple_groups <- length(unique(meta_input[[group]])) > 1L
    if (!multiple_groups) {
      warning("The supplied group column has fewer than two groups; fitting an ungrouped integration model")
      group <- NULL
    }
  }

  if (decompose) {
    if (!is.character(id) || length(id) != 1L || is.na(id) ||
        id == "" || !id %in% colnames(meta_input)) {
      stop("If decompose = TRUE, id must be provided as a column name in meta_input")
    }
    if (anyNA(meta_input[[id]]) || length(unique(meta_input[[id]])) < 2L) {
      stop("decompose = TRUE requires at least two non-missing subject identifiers")
    }
  }

  if (!is.null(variables)) {
    if (!is.list(variables) || is.null(names(variables)) ||
        anyNA(names(variables)) || any(names(variables) == "") ||
        anyDuplicated(names(variables))) {
      stop("variables must be a named list")
    }

    if (!all(names(variables) %in% names(data_input))) {
      stop("variables must have names matching data_input")
    }

    bad_vars <- !vapply(variables, function(z) {
      is.numeric(z) && length(z) == ncomp && all(is.finite(z)) &&
        all(z >= 1) && all(z %% 1 == 0)
    }, logical(1))
    if (any(bad_vars)) {
      stop("Each element of variables must contain ncomp positive integers")
    }
    too_large <- vapply(names(variables), function(nm) {
      any(variables[[nm]] > ncol(data_input[[nm]]))
    }, logical(1))
    if (any(too_large)) {
      stop(
        "variables cannot request more variables than a block contains. Problematic blocks: ",
        paste(names(variables)[too_large], collapse = ", ")
      )
    }
  }

  if (verbose) cat("Building experiment design ...\n")

  exp_design <- meta_input
  exp_design$inphase <- cos(2 * pi * exp_design[[time]] / period)
  exp_design$outphase <- sin(2 * pi * exp_design[[time]] / period)

  if (multiple_groups) {
    exp_design[[group]] <- droplevels(factor(exp_design[[group]]))
  }

  if (decompose) {
    if (verbose) cat("Performing within-subject decomposition ...\n")
    id_df <- as.data.frame(meta_input[[id]])
    names(id_df) <- id
    data_input <- lapply(data_input, function(x) mixOmics::withinVariation(x, id_df))
  }

  if (verbose) cat("Centering / scaling datasets ...\n")
  zero_variance <- vapply(data_input, function(x) {
    any(apply(x, 2L, stats::sd) <= sqrt(.Machine$double.eps))
  }, logical(1))
  if (any(zero_variance)) {
    stop(
      "Biological blocks contain constant variables after preprocessing: ",
      paste(names(data_input)[zero_variance], collapse = ", ")
    )
  }
  data_input <- lapply(
    data_input,
    function(x) {
      scale_values <- if (scale) apply(x, 2L, stats::sd) else FALSE
      base::scale(x, center = center, scale = scale_values)
    }
  )

  # ---------------------------------------------------------
  # Build design blocks
  # ---------------------------------------------------------
  if (verbose) cat("Building design blocks ...\n")

  Y_time <- NULL
  Y_group <- NULL
  Y_design <- NULL

  if (multiple_groups) {
    if (no_interaction) {
      # Common rhythm across groups, plus separate group block
      Y_time <- stats::model.matrix(~ 0 + inphase + outphase, data = exp_design)
      Y_group <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", group)), data = exp_design)
    } else {
      # Group-specific rhythm
      Y_time <- stats::model.matrix(
        stats::as.formula(paste0("~ 0 + ", group, ":inphase + ", group, ":outphase")),
        data = exp_design
      )
      Y_group <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", group)), data = exp_design)
    }

    if (merge) {
      if (is.null(Y_group)) {
        Y_design <- Y_time
      } else {
        Y_design <- cbind(Y_group, Y_time)
      }
      Y_time <- NULL
      Y_group <- NULL
    }
  } else {
    Y_time <- stats::model.matrix(~ 0 + inphase + outphase, data = exp_design)
  }

  # ---------------------------------------------------------
  # Prepare SGCCA input blocks
  # ---------------------------------------------------------
  blocks <- data_input

  if (!is.null(Y_design)) {
    blocks$Y_design <- Y_design
  } else {
    blocks$Y_time <- Y_time
    if (!is.null(Y_group)) {
      blocks$Y_group <- Y_group
    }
  }

  if (verbose) cat("Building SGCCA design matrix ...\n")

  block_names <- names(blocks)
  design_matrix <- matrix(0, nrow = length(block_names), ncol = length(block_names),
                          dimnames = list(block_names, block_names))
  diag(design_matrix) <- 0

  data_block_names <- names(data_input)

  if (!is.null(Y_design)) {
    design_matrix[data_block_names, "Y_design"] <- 1
    design_matrix["Y_design", data_block_names] <- 1
  } else {
    design_matrix[data_block_names, "Y_time"] <- 1
    design_matrix["Y_time", data_block_names] <- 1

    if (!is.null(Y_group)) {
      design_matrix[data_block_names, "Y_group"] <- 1
      design_matrix["Y_group", data_block_names] <- 1

      if (!no_correlation) {
        design_matrix["Y_time", "Y_group"] <- 1
        design_matrix["Y_group", "Y_time"] <- 1
      }
    }
  }

  # ---------------------------------------------------------
  # Build keepX
  # ---------------------------------------------------------
  if (is.null(variables)) {
    keepX <- lapply(blocks, function(x) rep(ncol(x), ncomp))
  } else {
    keepX <- vector("list", length(blocks))
    names(keepX) <- names(blocks)

    for (nm in names(data_input)) {
      if (nm %in% names(variables)) {
        keepX[[nm]] <- as.integer(variables[[nm]][seq_len(ncomp)])
      } else {
        keepX[[nm]] <- rep(ncol(blocks[[nm]]), ncomp)
      }
    }

    for (nm in setdiff(names(blocks), names(data_input))) {
      keepX[[nm]] <- rep(ncol(blocks[[nm]]), ncomp)
    }
  }

  # ---------------------------------------------------------
  # Fit SGCCA
  # ---------------------------------------------------------
  if (verbose) cat("Running SGCCA ...\n")

  sgcca_model <- tryCatch(
    mixOmics::wrapper.sgcca(
      X = blocks,
      ncomp = ncomp,
      scale = FALSE,
      design = design_matrix,
      mode = "canonical",
      keepX = keepX,
      max.iter = max.iter
    ),
    error = function(e) {
      stop(
        "SGCCA fitting failed: ", conditionMessage(e),
        ". Check block variation, component count, and design settings.",
        call. = FALSE
      )
    }
  )

  if (verbose) cat("Extracting scores and loadings ...\n")

  design_block_names <- intersect(c("Y_time", "Y_group", "Y_design"), names(sgcca_model$variates))
  data_block_names <- setdiff(names(sgcca_model$variates), design_block_names)

  partial_scores <- sgcca_model$variates[data_block_names]
  partial_loadings <- sgcca_model$loadings[data_block_names]

  # ---------------------------------------------------------
  # Average block scores across biological datasets
  # ---------------------------------------------------------
  if (verbose) cat("Calculating average scores ...\n")

  get_average_variates <- function(object, X_blocks) {
    arrays <- object$variates[X_blocks]
    summed <- Reduce(`+`, arrays)
    summed / length(arrays)
  }

  average_score <- get_average_variates(sgcca_model, data_block_names)

  output <- new(
    "mixedcirc_integration",
    partial = partial_scores,
    average = average_score,
    loadings = partial_loadings,
    model = sgcca_model
  )

  if (verbose) cat("Finished!\n")
  output
}

#' Extract scores, loadings, or model from a mixedcirc integration object
#'
#' Retrieves partial scores, average scores, loadings, or the full fitted SGCCA
#' model from an object of class `mixedcirc_integration`.
#'
#' @param input A `mixedcirc_integration` object.
#' @param type One of `"partial"`, `"average"`, `"loadings"`, or `"model"`.
#' @param dataset Optional dataset name to extract when `merge = FALSE` and
#'   `type` is `"partial"` or `"loadings"`.
#' @param merge Logical. If `TRUE`, partial scores are combined column-wise
#'   (samples by block-components), while loadings are combined row-wise
#'   (block-features by components). Block names are added to the merged
#'   column or row names, respectively. Default `FALSE`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default `FALSE`.
#'
#' @return
#' A matrix for `"partial"`, `"average"`, or `"loadings"`, or the full fitted
#' model object for `"model"`.
#'
#' @examples
#' library(mixedcirc)
#' data("circa_data")
#'
#' data_input <- list(
#'   a = circa_data$data_matrix[, 1:3],
#'   b = circa_data$data_matrix[, 4:7]
#' )
#'
#' results <- mixedcirc_integrate(
#'   data_input = data_input,
#'   meta_input = data.frame(
#'     time = circa_data$time,
#'     group = circa_data$group,
#'     id = circa_data$id
#'   ),
#'   time = "time",
#'   group = "group",
#'   id = "id"
#' )
#'
#' score_matrix <- mixedcirc_getscore(results, type = "average")
#'
#' @export
mixedcirc_getscore <- function(input = NULL,
                               type = NULL,
                               dataset = NULL,
                               merge = FALSE,
                               verbose = FALSE) {

  if (!is.logical(merge) || length(merge) != 1L || is.na(merge)) {
    stop("merge must be TRUE or FALSE")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  if (verbose) cat("Checking parameters ...\n")

  if (!is(input, "mixedcirc_integration")) {
    stop("input must be a mixedcirc_integration object")
  }

  if (is.null(type)) {
    stop("type must be provided")
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be one of: partial, average, loadings, model")
  }

  type <- match.arg(type, c("partial", "average", "loadings", "model"))

  if (type == "model") {
    return(input@model)
  }

  if (type == "average") {
    return(as.matrix(input@average))
  }

  slot_data <- switch(
    type,
    partial = input@partial,
    loadings = input@loadings
  )

  if (!is.list(slot_data) || length(slot_data) == 0) {
    stop("No data available for type = '", type, "'")
  }

  if (merge) {
    if (type == "partial") {
      row_counts <- vapply(slot_data, NROW, integer(1L))
      if (length(unique(row_counts)) != 1L) {
        stop("Partial-score matrices must have the same number of samples to merge")
      }
      merged <- Map(function(x, nm) {
        x <- as.matrix(x)
        component_names <- colnames(x)
        if (is.null(component_names)) {
          component_names <- paste0("comp", seq_len(ncol(x)))
        }
        colnames(x) <- paste0(nm, "_", component_names)
        x
      }, slot_data, names(slot_data))
      return(as.matrix(do.call(cbind, merged)))
    }

    component_counts <- vapply(slot_data, NCOL, integer(1L))
    if (length(unique(component_counts)) != 1L) {
      stop("Loading matrices must have the same number of components to merge")
    }
    merged <- Map(function(x, nm) {
      x <- as.matrix(x)
      feature_names <- rownames(x)
      if (is.null(feature_names)) {
        feature_names <- paste0("feature", seq_len(nrow(x)))
      }
      rownames(x) <- paste0(nm, "_", feature_names)
      x
    }, slot_data, names(slot_data))
    return(as.matrix(do.call(rbind, merged)))
  }

  if (is.null(dataset)) {
    stop("dataset must be provided when merge = FALSE and type is 'partial' or 'loadings'")
  }
  if (!is.character(dataset) || length(dataset) != 1L || is.na(dataset) || dataset == "") {
    stop("dataset must be one block name")
  }

  if (!dataset %in% names(slot_data)) {
    stop("dataset must be one of: ", paste(names(slot_data), collapse = ", "))
  }

  if (verbose) cat("Extracting data ...\n")
  as.matrix(slot_data[[dataset]])
}
