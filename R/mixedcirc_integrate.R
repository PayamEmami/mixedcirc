#' Integrate circadian datasets using SGCCA
#'
#' Performs multi-block data integration using sparse generalized canonical
#' correlation analysis (SGCCA). The biological datasets are
#' integrated together with one or two design blocks representing circadian
#' rhythm and, optionally, group structure.
#'
#' @param data_input A named list of numeric matrices or data.frames. Each element
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
#'   per component in each dataset. Each element should be an integer vector of
#'   length `ncomp`. Names must match `data_input`.
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
#'   separate blocks, the covariance between the group-design block and the
#'   rhythm-design block is fixed to zero. Default `FALSE`.
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
#' The resulting partial or average scores can be used as integrated features in
#' downstream analyses such as `mixedcirc_detect()`.
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
#' @import stats
#' @import methods
#' @import mixOmics
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

  if (verbose) cat("Checking inputs ...\n")

  if (is.null(data_input) || !is.list(data_input) || length(data_input) == 0) {
    stop("data_input must be a non-empty named list of matrices or data.frames")
  }

  if (is.null(names(data_input)) || any(names(data_input) == "")) {
    stop("data_input must be a named list")
  }

  ok_blocks <- vapply(data_input, function(x) is.matrix(x) || is.data.frame(x), logical(1))
  if (any(!ok_blocks)) {
    stop(
      "All elements of data_input must be matrices or data.frames. Problematic blocks: ",
      paste(names(data_input)[!ok_blocks], collapse = ", ")
    )
  }

  data_input <- lapply(data_input, as.matrix)

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

  if (is.null(time) || !time %in% colnames(meta_input)) {
    stop("time must be the name of a column in meta_input")
  }

  if (!is.numeric(period) || length(period) != 1 || is.na(period) || is.infinite(period)) {
    stop("period must be a single finite numeric value")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1 || is.na(ncomp) || ncomp <= 0 || ncomp %% 1 != 0) {
    stop("ncomp must be a positive integer")
  }

  multiple_groups <- FALSE
  if (!is.null(group)) {
    if (!group %in% colnames(meta_input)) {
      stop("group must be the name of a column in meta_input")
    }
    multiple_groups <- TRUE
  }

  if (decompose) {
    if (is.null(id) || !id %in% colnames(meta_input)) {
      stop("If decompose = TRUE, id must be provided as a column name in meta_input")
    }
  }

  if (!is.null(variables)) {
    if (!is.list(variables) || is.null(names(variables))) {
      stop("variables must be a named list")
    }

    if (!all(names(variables) %in% names(data_input))) {
      stop("variables must have names matching data_input")
    }

    bad_vars <- !vapply(variables, is.numeric, logical(1))
    if (any(bad_vars)) {
      stop("Each element of variables must be a numeric vector")
    }

    short_vars <- vapply(variables, length, integer(1)) < ncomp
    if (any(short_vars)) {
      warning(
        "Each element of variables should have length ncomp = ", ncomp,
        ". Using all available variables for: ",
        paste(names(variables)[short_vars], collapse = ", ")
      )
    }
  }

  if (verbose) cat("Building experiment design ...\n")

  exp_design <- meta_input
  exp_design$inphase <- cos(2 * pi * exp_design[[time]] / period)
  exp_design$outphase <- sin(2 * pi * exp_design[[time]] / period)

  if (multiple_groups) {
    exp_design[[group]] <- factor(exp_design[[group]])
  }

  if (decompose) {
    if (verbose) cat("Performing within-subject decomposition ...\n")
    id_df <- as.data.frame(meta_input[[id]])
    names(id_df) <- id
    data_input <- lapply(data_input, function(x) mixOmics::withinVariation(x, id_df))
  }

  if (verbose) cat("Centering / scaling datasets ...\n")
  data_input <- lapply(data_input, function(x) scale(x, center = center, scale = scale))

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
      if (nm %in% names(variables) && length(variables[[nm]]) >= ncomp) {
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

  sgcca_model <- mixOmics::wrapper.sgcca(
    X = blocks,
    ncomp = ncomp,
    scale = FALSE,
    design = design_matrix,
    mode = "canonical",
    keepX = keepX,
    max.iter = max.iter
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
#' @param merge Logical. If `TRUE`, all datasets are combined column-wise into one
#'   matrix when `type` is `"partial"` or `"loadings"`. Default `FALSE`.
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
#' @import methods
#' @export
mixedcirc_getscore <- function(input = NULL,
                               type = NULL,
                               dataset = NULL,
                               merge = FALSE,
                               verbose = FALSE) {

  if (verbose) cat("Checking parameters ...\n")

  if (!is(input, "mixedcirc_integration")) {
    stop("input must be a mixedcirc_integration object")
  }

  if (is.null(type)) {
    stop("type must be provided")
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
    merged <- Map(function(x, nm) {
      x <- as.matrix(x)
      colnames(x) <- paste0(nm, "_", colnames(x))
      x
    }, slot_data, names(slot_data))

    return(as.matrix(do.call(cbind, merged)))
  }

  if (is.null(dataset)) {
    stop("dataset must be provided when merge = FALSE and type is 'partial' or 'loadings'")
  }

  if (!dataset %in% names(slot_data)) {
    stop("dataset must be one of: ", paste(names(slot_data), collapse = ", "))
  }

  if (verbose) cat("Extracting data ...\n")
  as.matrix(slot_data[[dataset]])
}
