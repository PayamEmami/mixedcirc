#' Adjust p-values across fitted features
#'
#' Applies a multiple-testing correction to matching result columns in a
#' `mixedcirc_fit_list`. Each selected column is adjusted separately across
#' features, and the adjusted values are appended to the corresponding feature
#' results.
#'
#' @param input A `mixedcirc_fit_list` containing at least two fitted features.
#' @param method A multiple-testing method accepted by [stats::p.adjust()].
#' @param pattern A regular expression used to select result columns. The
#'   default, `"p.value"`, selects the p-value columns produced by
#'   [mixedcirc_detect()].
#' @param verbose Logical; print progress messages.
#' @param ... Reserved for future use.
#' @examples
#' simulated <- mixedcirc_simulate_data(
#'   n_subjects_per_group = 3,
#'   n_timepoints = 6,
#'   n_features = 3,
#'   groups = c("A", "B"),
#'   seed = 10
#' )
#' results <- mixedcirc_detect(
#'   data_input = simulated$data_input,
#'   meta_input = simulated$meta_input,
#'   time = "time",
#'   group = "group",
#'   id = "id",
#'   lm_method = "lm",
#'   random_effect_type = "none"
#' )
#' results <- mixedcirc_adjust(results, method = "bonferroni")
#'
#' @return A `mixedcirc_fit_list` with adjusted p-value columns appended to each
#'   feature's result table.
#'
#' @export
mixedcirc_adjust <- function(input, method = "BH", pattern = "p.value",
                             verbose = FALSE, ...) {

  if (isTRUE(verbose)) {
    cat("Checking inputs ...\n")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }

  if (!is(input, "mixedcirc_fit_list")) {
    stop("input must be a mixedcirc_fit_list object", call. = FALSE)
  }

  if (length(input) < 2) {
    stop("input must contain at least two fitted features", call. = FALSE)
  }

  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
    stop(
      "method must be one of: 'holm', 'hochberg', 'hommel', 'bonferroni', ",
      "'BH', 'BY', 'fdr', or 'none'",
      call. = FALSE
    )
  }

  if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern) || pattern == "") {
    stop("pattern must be one non-empty character string", call. = FALSE)
  }

  if (verbose) {
    cat("Gathering the data ...\n")
  }

  data_gathered <- input@results
  data_gathered <- do.call(rbind, lapply(data_gathered, function(x) { x@results }))

  selected_pvalues <- grepl(pattern = pattern, x = colnames(data_gathered))

  if (sum(selected_pvalues) == 0) {
    stop("no result column matches pattern", call. = FALSE)
  }

  selected_data <- data_gathered[, selected_pvalues, drop = FALSE]

  if (verbose) {
    cat("Performing adjustment ...\n")
  }

  adjusted_pvalues <- vapply(
    selected_data,
    function(x) stats::p.adjust(p = x, method = method),
    numeric(nrow(selected_data))
  )
  adjusted_pvalues <- matrix(
    adjusted_pvalues,
    nrow = nrow(selected_data),
    dimnames = list(rownames(selected_data), colnames(selected_data))
  )

  colnames(adjusted_pvalues) <- paste0(colnames(selected_data), "_adjusted_", method)

  if (verbose) {
    cat("Preparing the results ...\n")
  }

  for (i in seq_along(input)) {
    input@results[[i]]@results <- cbind(
      input@results[[i]]@results[, , drop = FALSE],
      adjusted_pvalues[i, , drop = FALSE]
    )
  }

  if (verbose) {
    cat("Finished\n")
  }

  input
}
