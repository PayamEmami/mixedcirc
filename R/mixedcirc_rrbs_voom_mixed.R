# Chunked methylBaseDB weights ----------------------------------------------

.mc_methylkit_internals <- function(functions = c(
    "applyTbxByChunk", "makeTabixHeader", ".formatTabixHeader",
    "readMethylBaseDB"
  )) {
  if (!requireNamespace("methylKit", quietly = TRUE)) {
    stop(
      "mixedcirc_rrbs_voom_mixed() requires the Bioconductor package methylKit.",
      call. = FALSE
    )
  }

  namespace <- asNamespace("methylKit")
  missing <- functions[!vapply(
    functions,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )]
  if (length(missing)) {
    stop(
      "The installed methylKit version is incompatible with the chunked RRBS ",
      "workflow; missing internal function(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  stats::setNames(
    lapply(functions, get, envir = namespace, inherits = FALSE),
    functions
  )
}

.mc_methylbase_slots <- function(object,
                                 sample.ids = object@sample.ids,
                                 treatment = object@treatment,
                                 coverage.index = object@coverage.index,
                                 numCs.index = object@numCs.index,
                                 numTs.index = object@numTs.index) {
  list(
    dbtype = object@dbtype,
    sample.ids = sample.ids,
    assembly = object@assembly,
    context = object@context,
    resolution = object@resolution,
    treatment = treatment,
    destranded = object@destranded,
    coverage.index = coverage.index,
    numCs.index = numCs.index,
    numTs.index = numTs.index,
    num.records = object@num.records
  )
}

.mc_read_methylbase <- function(path, object, internals,
                                sample.ids = object@sample.ids,
                                treatment = object@treatment) {
  internals$readMethylBaseDB(
    dbpath = path,
    sample.ids = sample.ids,
    assembly = object@assembly,
    dbtype = "tabix",
    resolution = object@resolution,
    context = object@context,
    treatment = treatment,
    destranded = object@destranded
  )
}

#' Transform an on-disk RRBS data set for mixed modelling (experimental)
#'
#' Estimates observation-level precision weights for a `methylKit`
#' `methylBaseDB` object without loading the full locus-by-sample matrix into
#' memory. Methylated and unmethylated counts are modelled as paired
#' observations, and the mean-variance trend is estimated in chunks.
#'
#' @param counts A `methylKit::methylBaseDB` object.
#' @param formula A one-sided fixed- or mixed-effects formula describing the
#'   paired methylated/unmethylated observations.
#' @param data A data frame with one row for every methylated or unmethylated
#'   observation. Its rows must follow the order sample 1 methylated, sample 1
#'   unmethylated, sample 2 methylated, sample 2 unmethylated, and so on.
#' @param lib.size Optional library sizes. Supply either one value per
#'   biological sample or one value per paired model observation. When omitted,
#'   library sizes are calculated from the coverage columns.
#' @param chunk.size Positive integer number of loci processed per chunk.
#' @param span Width of the lowess smoothing window as a proportion.
#' @param plot Logical; draw the mean-variance trend.
#' @param save.plot Logical; include `voom.xy` and `voom.line` in the result.
#' @param quiet Logical; suppress model-fitting progress messages.
#' @param BPPARAM A `BiocParallelParam` object used while fitting locus models.
#' @param ignore_missing Logical. If `FALSE`, missing counts are replaced by
#'   zero for compatibility with the original workflow. If `TRUE`, models use
#'   the available observations and retain missing fitted values.
#' @param ... Additional arguments passed to `lm()` or `lmer()` through
#'   `variance_fit()` (internal).
#'
#' @return A list with `E` and `weights`, both on-disk `methylBaseDB` objects.
#'   `E` stores log2 counts per million in the methylated/unmethylated count
#'   columns; `weights` stores the corresponding precision weights. If
#'   `save.plot = TRUE`, the list also contains `voom.xy` and `voom.line`.
#'
#' @details
#' This is a compatibility workflow for very large RRBS data sets already
#' stored by `methylKit`. It relies on methylKit's internal tabix chunking API,
#' which mixedcirc checks explicitly at run time. The method is experimental;
#' generated files are written beside the input database.
#'
#' A typical circadian design contains `scaler` (one for methylated and zero for
#' unmethylated observations), harmonic `inphase` and `outphase` variables, and a sample
#' or subject identifier. For example,
#' `~ scaler + scaler:inphase + scaler:outphase + (1 | replicate_id)`.
#'
#' @export
mixedcirc_rrbs_voom_mixed <- function(
    counts,
    formula,
    data,
    lib.size = NULL,
    chunk.size = 100,
    span = 0.5,
    plot = FALSE,
    save.plot = FALSE,
    quiet = FALSE,
    BPPARAM = BiocParallel::SerialParam(),
    ignore_missing = FALSE,
    ...) {
  internals <- .mc_methylkit_internals()
  if (!methods::is(counts, "methylBaseDB")) {
    stop("counts must be a methylKit methylBaseDB object", call. = FALSE)
  }

  formula <- stats::as.formula(formula)
  if (length(formula) != 2L) {
    stop("formula must be one-sided", call. = FALSE)
  }
  data <- as.data.frame(data)
  sample_count <- length(counts@sample.ids)
  observation_count <- 2L * sample_count
  if (nrow(data) != observation_count) {
    stop(
      "data must have exactly two rows per methylBaseDB sample (methylated then unmethylated)",
      call. = FALSE
    )
  }
  missing_variables <- setdiff(all.vars(formula), names(data))
  if (length(missing_variables)) {
    stop(
      "Variable in formula not found in data: ",
      paste(missing_variables, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(chunk.size) || length(chunk.size) != 1L ||
      !is.finite(chunk.size) || chunk.size < 1 || chunk.size %% 1 != 0) {
    stop("chunk.size must be a positive integer", call. = FALSE)
  }
  chunk.size <- as.integer(chunk.size)
  if (!is.numeric(span) || length(span) != 1L || !is.finite(span) ||
      span <= 0 || span > 1) {
    stop("span must be a number greater than zero and no greater than one", call. = FALSE)
  }
  logical_arguments <- list(
    plot = plot,
    save.plot = save.plot,
    quiet = quiet,
    ignore_missing = ignore_missing
  )
  invalid_logical <- !vapply(
    logical_arguments,
    function(value) is.logical(value) && length(value) == 1L && !is.na(value),
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(names(logical_arguments)[which(invalid_logical)[1L]],
         " must be TRUE or FALSE", call. = FALSE)
  }
  if (counts@num.records < 2L) {
    stop("Need at least two loci to fit a mean-variance trend", call. = FALSE)
  }

  source_path <- methylKit::getDBPath(counts)
  output_directory <- dirname(source_path)
  source_stem <- tools::file_path_sans_ext(basename(source_path))
  value_indices <- c(rbind(counts@numCs.index, counts@numTs.index))
  original_column_count <- max(c(
    4L, counts@coverage.index, counts@numCs.index, counts@numTs.index
  ))

  raw_totals <- internals$applyTbxByChunk(
    tbxFile = source_path,
    chunk.size = chunk.size,
    return.type = "data.frame",
    FUN = function(chunk) {
      values <- as.matrix(chunk[, value_indices, drop = FALSE])
      data.frame(total = rowSums(values, na.rm = TRUE))
    }
  )$total

  if (is.null(lib.size)) {
    library_parts <- internals$applyTbxByChunk(
      tbxFile = source_path,
      chunk.size = chunk.size,
      return.type = "data.frame",
      FUN = function(chunk) {
        coverage <- as.matrix(chunk[, counts@coverage.index, drop = FALSE])
        if (!ignore_missing) coverage[is.na(coverage)] <- 0
        as.data.frame(t(colSums(coverage, na.rm = TRUE)), check.names = FALSE)
      }
    )
    sample_library_sizes <- colSums(library_parts, na.rm = TRUE)
    lib.size <- rep(sample_library_sizes, each = 2L)
  } else {
    if (!is.numeric(lib.size) || any(!is.finite(lib.size)) || any(lib.size < 0)) {
      stop("lib.size must contain finite non-negative values", call. = FALSE)
    }
    if (length(lib.size) == sample_count) {
      lib.size <- rep(lib.size, each = 2L)
    } else if (length(lib.size) != observation_count) {
      stop(
        "lib.size must have one value per sample or one value per paired model observation",
        call. = FALSE
      )
    }
  }

  base_slots <- .mc_methylbase_slots(counts)
  base_header <- internals$makeTabixHeader(base_slots)
  base_header <- internals$.formatTabixHeader(
    class = "methylBaseDB",
    tabixHead = base_header
  )

  normalized_name <- paste0(source_stem, "_mixedcirc_libnorm")
  normalized_path <- internals$applyTbxByChunk(
    tbxFile = source_path,
    tabixHead = base_header,
    chunk.size = chunk.size,
    return.type = "tabix",
    dir = output_directory,
    filename = normalized_name,
    FUN = function(chunk) {
      values <- as.matrix(chunk[, value_indices, drop = FALSE])
      if (!ignore_missing) values[is.na(values)] <- 0
      chunk[, value_indices] <- t(
        log2(t(values + 0.5) / (lib.size + 1) * 1e6)
      )
      chunk
    }
  )
  normalized <- .mc_read_methylbase(
    normalized_path, counts, internals
  )

  auxiliary_names <- c("sigma", "resid", "Ameans")
  maximum_index <- original_column_count
  statistics_slots <- .mc_methylbase_slots(
    counts,
    sample.ids = c(counts@sample.ids, auxiliary_names),
    treatment = c(counts@treatment, rep(NA_real_, length(auxiliary_names))),
    coverage.index = c(
      counts@coverage.index,
      maximum_index + c(1L, 4L, 7L)
    ),
    numCs.index = c(
      counts@numCs.index,
      maximum_index + c(2L, 5L, 8L)
    ),
    numTs.index = c(
      counts@numTs.index,
      maximum_index + c(3L, 6L, 9L)
    )
  )
  statistics_header <- internals$.formatTabixHeader(
    class = "methylBaseDB",
    tabixHead = internals$makeTabixHeader(statistics_slots)
  )
  fit_dots <- list(...)
  if (ignore_missing && is.null(fit_dots$na.action)) {
    fit_dots$na.action <- stats::na.exclude
  }

  statistics_name <- paste0(source_stem, "_mixedcirc_stats")
  statistics_path <- internals$applyTbxByChunk(
    tbxFile = methylKit::getDBPath(normalized),
    tabixHead = statistics_header,
    chunk.size = chunk.size,
    return.type = "tabix",
    dir = output_directory,
    filename = statistics_name,
    FUN = function(chunk) {
      expression <- as.matrix(chunk[, value_indices, drop = FALSE])
      average_log_cpm <- rowMeans(expression, na.rm = TRUE)
      model_arguments <- c(
        list(
          exprObj = expression,
          formula = formula,
          data = data,
          showWarnings = FALSE,
          fxn = function(fit) {
            list(
              sd = stats::sigma(fit),
              df.residual = stats::df.residual(fit),
              fitted.values = stats::fitted(fit)
            )
          },
          BPPARAM = BPPARAM,
          ignore_na = ignore_missing,
          quiet = quiet
        ),
        fit_dots
      )
      fitted_models <- do.call(variance_fit, model_arguments)
      fitted_values <- do.call(
        rbind,
        lapply(fitted_models, function(model) model$fitted.values)
      )
      if (!identical(dim(fitted_values), dim(expression))) {
        stop(
          "Fitted RRBS values did not retain the input observation layout",
          call. = FALSE
        )
      }
      sigma_values <- vapply(fitted_models, `[[`, numeric(1), "sd")
      residual_df <- vapply(
        fitted_models, `[[`, numeric(1), "df.residual"
      )

      chunk[, value_indices] <- fitted_values
      cbind(
        chunk,
        sigma_values, sigma_values, sigma_values,
        residual_df, residual_df, residual_df,
        average_log_cpm, average_log_cpm, average_log_cpm
      )
    }
  )

  statistics <- .mc_read_methylbase(
    statistics_path,
    counts,
    internals,
    sample.ids = statistics_slots$sample.ids,
    treatment = statistics_slots$treatment
  )
  real_samples <- !statistics@sample.ids %in% auxiliary_names
  sigma_column <- statistics@numCs.index[statistics@sample.ids == "sigma"]
  average_column <- statistics@numCs.index[statistics@sample.ids == "Ameans"]
  statistic_summary <- internals$applyTbxByChunk(
    tbxFile = methylKit::getDBPath(statistics),
    chunk.size = chunk.size,
    return.type = "data.frame",
    FUN = function(chunk) {
      data.frame(
        sigma = chunk[, sigma_column],
        average = chunk[, average_column]
      )
    }
  )

  trend_x <- statistic_summary$average +
    mean(log2(lib.size + 1)) - log2(1e6)
  trend_y <- sqrt(statistic_summary$sigma)
  keep <- raw_totals > 0 & is.finite(trend_x) & is.finite(trend_y)
  if (sum(keep) < 2L || length(unique(trend_x[keep])) < 2L) {
    stop(
      "Too few informative loci to estimate the RRBS mean-variance trend",
      call. = FALSE
    )
  }
  lowess_line <- stats::lowess(trend_x[keep], trend_y[keep], f = span)
  trend_function <- stats::approxfun(lowess_line, rule = 2)

  if (plot) {
    graphics::plot(
      trend_x[keep], trend_y[keep],
      xlab = "log2(count size + 0.5)",
      ylab = "Sqrt(standard deviation)",
      pch = 16,
      cex = 0.25
    )
    graphics::title("voom: Mean-variance trend")
    graphics::lines(lowess_line, col = "red")
  }

  fitted_indices <- c(rbind(
    statistics@numCs.index[real_samples],
    statistics@numTs.index[real_samples]
  ))
  weights_name <- paste0(source_stem, "_mixedcirc_weights")
  weights_path <- internals$applyTbxByChunk(
    tbxFile = methylKit::getDBPath(statistics),
    tabixHead = base_header,
    chunk.size = chunk.size,
    return.type = "tabix",
    dir = output_directory,
    filename = weights_name,
    FUN = function(chunk) {
      fitted_cpm <- 2^as.matrix(chunk[, fitted_indices, drop = FALSE])
      fitted_count <- 1e-6 * t(t(fitted_cpm) * (lib.size + 1))
      fitted_log_count <- log2(pmax(fitted_count, .Machine$double.xmin))
      weights <- 1 / trend_function(fitted_log_count)^4
      if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop("The fitted RRBS trend produced invalid precision weights", call. = FALSE)
      }

      result <- chunk[, seq_len(original_column_count), drop = FALSE]
      result[, value_indices] <- weights
      result
    }
  )
  weights <- .mc_read_methylbase(weights_path, counts, internals)

  result <- list(E = normalized, weights = weights)
  if (save.plot) {
    result$voom.xy <- list(x = trend_x[keep], y = trend_y[keep])
    result$voom.line <- lowess_line
  }
  result
}
