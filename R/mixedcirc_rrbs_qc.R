# RRBS quality control -------------------------------------------------------

#' Summarize paired RRBS count quality without changing the analysis
#'
#' Checks the paired methylated/unmethylated representation expected by
#' [mixedcirc_detect()] and summarizes coverage and methylation fractions. The
#' data are not transformed, filtered, reordered or otherwise modified.
#'
#' @param data_input Numeric matrix or data frame whose consecutive rows are
#'   methylated and unmethylated counts for each biological observation.
#' @param meta_input Sample metadata with one row per biological observation.
#' @param replicate_id Name of the unique sample/replicate identifier column.
#' @param min_coverage Coverage threshold used only for QC summaries. The edgeR
#'   RRBS workflow commonly uses eight reads for CpG-level analyses.
#'
#' @return A list of class `mixedcirc_rrbs_qc` with `feature`, `sample`, and
#'   `pairing` summaries. Input data are never modified.
#' @export
mixedcirc_rrbs_qc <- function(data_input,
                              meta_input,
                              replicate_id,
                              min_coverage = 8) {
  if (!is.matrix(data_input) && !is.data.frame(data_input)) {
    stop("data_input must be a numeric matrix or data frame", call. = FALSE)
  }
  counts <- as.matrix(data_input)
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("data_input must contain finite non-negative counts", call. = FALSE)
  }
  if (nrow(counts) %% 2L != 0L) {
    stop("data_input must have an even number of paired rows", call. = FALSE)
  }
  meta_input <- as.data.frame(meta_input)
  if (!is.character(replicate_id) || length(replicate_id) != 1L ||
      is.na(replicate_id) || !replicate_id %in% names(meta_input)) {
    stop("replicate_id must name a column in meta_input", call. = FALSE)
  }
  if (nrow(meta_input) * 2L != nrow(counts)) {
    stop("meta_input must contain one row per methylated/unmethylated pair", call. = FALSE)
  }
  if (anyNA(meta_input[[replicate_id]]) || anyDuplicated(meta_input[[replicate_id]])) {
    stop("replicate_id must be complete and unique", call. = FALSE)
  }
  if (!is.numeric(min_coverage) || length(min_coverage) != 1L ||
      !is.finite(min_coverage) || min_coverage < 0) {
    stop("min_coverage must be one non-negative number", call. = FALSE)
  }

  methylated <- counts[seq.int(1L, nrow(counts), by = 2L), , drop = FALSE]
  unmethylated <- counts[seq.int(2L, nrow(counts), by = 2L), , drop = FALSE]
  coverage <- methylated + unmethylated
  fraction <- methylated / coverage
  fraction[coverage == 0] <- NA_real_
  feature_names <- colnames(counts)
  if (is.null(feature_names)) feature_names <- paste0("feature_", seq_len(ncol(counts)))
  feature_names <- make.unique(as.character(feature_names))
  sample_names <- as.character(meta_input[[replicate_id]])

  feature_summary <- data.frame(
    feature = feature_names,
    n_samples = nrow(coverage),
    minimum_coverage = apply(coverage, 2L, min),
    median_coverage = apply(coverage, 2L, stats::median),
    mean_coverage = colMeans(coverage),
    maximum_coverage = apply(coverage, 2L, max),
    proportion_covered = colMeans(coverage >= min_coverage),
    mean_methylation = colMeans(fraction, na.rm = TRUE),
    zero_coverage_pairs = colSums(coverage == 0),
    usable_all_samples = colSums(coverage >= min_coverage) == nrow(coverage),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  sample_summary <- data.frame(
    replicate_id = sample_names,
    minimum_coverage = apply(coverage, 1L, min),
    median_coverage = apply(coverage, 1L, stats::median),
    mean_coverage = rowMeans(coverage),
    maximum_coverage = apply(coverage, 1L, max),
    proportion_covered = rowMeans(coverage >= min_coverage),
    mean_methylation = rowMeans(fraction, na.rm = TRUE),
    zero_coverage_features = rowSums(coverage == 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  input_row_names <- rownames(counts)
  row_name_verification <- NA
  if (!is.null(input_row_names)) {
    methylated_names <- input_row_names[seq.int(1L, nrow(counts), by = 2L)]
    unmethylated_names <- input_row_names[seq.int(2L, nrow(counts), by = 2L)]
    methylated_base <- sub("([_.-])(methylated|meth)$", "", methylated_names,
                           ignore.case = TRUE)
    unmethylated_base <- sub("([_.-])(unmethylated|unmeth)$", "", unmethylated_names,
                             ignore.case = TRUE)
    recognized <- methylated_base != methylated_names & unmethylated_base != unmethylated_names
    if (all(recognized)) row_name_verification <- identical(methylated_base, unmethylated_base)
  }

  structure(
    list(
      feature = feature_summary,
      sample = sample_summary,
      pairing = list(
        n_pairs = nrow(coverage),
        n_features = ncol(coverage),
        alternating_rows_assumed = TRUE,
        row_names_verify_pairing = row_name_verification,
        counts_are_integer = all(counts == round(counts)),
        min_coverage = min_coverage
      )
    ),
    class = "mixedcirc_rrbs_qc"
  )
}

#' Print an RRBS quality-control summary
#'
#' @param x An object returned by [mixedcirc_rrbs_qc()].
#' @param ... Additional arguments; currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.mixedcirc_rrbs_qc <- function(x, ...) {
  cat("mixedcirc RRBS quality control\n")
  cat("------------------------------\n")
  cat("Paired samples:", x$pairing$n_pairs, "\n")
  cat("Features:", x$pairing$n_features, "\n")
  cat("Coverage threshold:", x$pairing$min_coverage, "\n")
  cat(
    "Features covered in all samples:",
    sum(x$feature$usable_all_samples), "of", nrow(x$feature), "\n"
  )
  if (!is.na(x$pairing$row_names_verify_pairing)) {
    cat("Pair order verified from row names:", x$pairing$row_names_verify_pairing, "\n")
  } else {
    cat("Pair order: assumed (row names did not encode methylated/unmethylated labels)\n")
  }
  invisible(x)
}

#' Plot RRBS coverage and methylation quality summaries
#'
#' @param x An object returned by [mixedcirc_rrbs_qc()].
#' @param type Summary to plot.
#'
#' @return A ggplot object.
#' @export
mixedcirc_plot_rrbs_qc <- function(
    x,
    type = c("feature_coverage", "sample_coverage", "methylation")) {
  if (!inherits(x, "mixedcirc_rrbs_qc")) {
    stop("x must be returned by mixedcirc_rrbs_qc", call. = FALSE)
  }
  type <- match.arg(type)
  threshold <- x$pairing$min_coverage
  if (type == "feature_coverage") {
    return(
      ggplot2::ggplot(
        x$feature,
        ggplot2::aes(x = .data$median_coverage, fill = .data$usable_all_samples)
      ) +
        ggplot2::geom_histogram(bins = 30, color = "white") +
        ggplot2::geom_vline(xintercept = threshold, linetype = 2) +
        ggplot2::labs(
          x = "Median feature coverage", y = "Features",
          fill = "Covered in all samples", title = "RRBS feature coverage"
        ) +
        ggplot2::theme_classic(base_size = 11)
    )
  }
  if (type == "sample_coverage") {
    return(
      ggplot2::ggplot(
        x$sample,
        ggplot2::aes(x = .data$median_coverage, y = .data$proportion_covered)
      ) +
        ggplot2::geom_point(color = "#2166AC", alpha = 0.8) +
        ggplot2::labs(
          x = "Median sample coverage",
          y = paste0("Proportion features with coverage >= ", threshold),
          title = "RRBS sample coverage"
        ) +
        ggplot2::theme_classic(base_size = 11)
    )
  }
  ggplot2::ggplot(x$feature, ggplot2::aes(x = .data$mean_methylation)) +
    ggplot2::geom_histogram(bins = 30, color = "white", fill = "#4D9221") +
    ggplot2::labs(
      x = "Mean methylation fraction", y = "Features",
      title = "RRBS methylation distribution"
    ) +
    ggplot2::theme_classic(base_size = 11)
}
