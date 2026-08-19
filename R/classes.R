#' Bootstrap results from a mixedcirc model
#'
#' @slot results A data frame containing bootstrap estimates, confidence
#'   intervals, standard errors, and p-values.
#' @slot metadata A list describing the bootstrap procedure.
#' @slot original_model The `mixedcirc_fit` object that was bootstrapped.
#'
#' @exportClass mixedcirc_boot
mixedcirc_boot <- methods::setClass(
  "mixedcirc_boot",
  slots = list(
    results = "data.frame",
    metadata = "list",
    original_model = "ANY"
  )
)

#' @param object A `mixedcirc_boot` object.
#' @rdname mixedcirc_boot-class
#' @export
methods::setMethod("show", "mixedcirc_boot", function(object) {
  cat("mixedcirc bootstrap result\n")
  cat("--------------------------\n")
  print(object@results)
  if (length(object@metadata)) {
    cat("\nMetadata:\n")
    print(object@metadata)
  }
  invisible(object)
})

#' @rdname mixedcirc_boot-class
#' @export
methods::setMethod("summary", "mixedcirc_boot", function(object, ...) {
  cat("Summary of mixedcirc bootstrap result\n")
  if (!is.null(object@metadata$nsim)) {
    cat("Simulations:", object@metadata$nsim, "\n")
  }
  if (!is.null(object@metadata$conf_level)) {
    cat("Confidence level:", object@metadata$conf_level, "\n")
  }
  print(object@results, ...)
  invisible(object)
})

#' @param x A `mixedcirc_boot` object.
#' @param row.names,optional Retained for compatibility with
#'   [as.data.frame()]; currently ignored.
#' @param ... Additional arguments; currently ignored by `as.data.frame()`.
#' @rdname mixedcirc_boot-class
#' @export
methods::setMethod("as.data.frame", "mixedcirc_boot", function(x,
                                                                row.names = NULL,
                                                                optional = FALSE,
                                                                ...) {
  x@results
})

#' Result from fitting one circadian feature
#'
#' @slot results A one-row data frame containing estimates and test results.
#' @slot fit The fitted `lm` or `lmerMod` model.
#' @slot exp_design The experimental design created by `mixedcirc_detect()`.
#' @slot type Analysis type, currently `"expression"` or `"RRBS"`.
#' @slot call The call used to create the object.
#' @slot params Additional parameters used by downstream methods.
#'
#' @exportClass mixedcirc_fit
mixedcirc_fit <- methods::setClass(
  "mixedcirc_fit",
  slots = list(
    results = "data.frame",
    fit = "ANY",
    exp_design = "data.frame",
    type = "character",
    call = "call",
    params = "list"
  )
)

#' Results from fitting multiple circadian features
#'
#' @slot results A list of `mixedcirc_fit` objects.
#'
#' @exportClass mixedcirc_fit_list
mixedcirc_fit_list <- methods::setClass(
  "mixedcirc_fit_list",
  slots = list(results = "list")
)

#' Print a single-feature result
#'
#' @param x A `mixedcirc_fit` object.
#' @param ... Additional arguments passed to [print()].
#'
#' @return `x`, invisibly.
#' @export
print.mixedcirc_fit <- function(x, ...) {
  print(x@results, ...)
  invisible(x)
}

#' @param object A `mixedcirc_fit` object.
#' @rdname mixedcirc_fit-class
#' @export
methods::setMethod("show", "mixedcirc_fit", function(object) {
  print(object@results)
  invisible(object)
})

#' Print a multi-feature result
#'
#' @param x A `mixedcirc_fit_list` object.
#' @param ... Additional arguments; currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.mixedcirc_fit_list <- function(x, ...) {
  cat("Total number of variables:", length(x@results), "\n")
  invisible(x)
}

#' @param object A `mixedcirc_fit_list` object.
#' @rdname mixedcirc_fit_list-class
#' @export
methods::setMethod("show", "mixedcirc_fit_list", function(object) {
  print.mixedcirc_fit_list(object)
})

#' Experimental integration result for multiple data sets
#'
#' @slot partial A named list of sample score matrices, one per data block.
#' @slot average A matrix of average scores across blocks.
#' @slot loadings A named list of feature loading matrices, one per data block.
#' @slot model The fitted SGCCA model.
#'
#' @details
#' This class and the integration workflow that creates it are experimental.
#'
#' @exportClass mixedcirc_integration
mixedcirc_integration <- methods::setClass(
  "mixedcirc_integration",
  slots = list(
    partial = "list",
    average = "matrix",
    loadings = "list",
    model = "ANY"
  )
)

#' @param object A `mixedcirc_integration` object.
#' @rdname mixedcirc_integration-class
#' @export
methods::setMethod("show", "mixedcirc_integration", function(object) {
  cat("Experimental mixedcirc integration result\n")
  cat("Omics blocks:", length(object@partial), "\n")
  cat("Use for exploratory analysis; validate learned scores on independent data.\n")
  invisible(object)
})

#' Print an experimental integration result
#'
#' @param x A `mixedcirc_integration` object.
#' @param ... Additional arguments; currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.mixedcirc_integration <- function(x, ...) {
  methods::show(x)
  invisible(x)
}

#' Subset mixedcirc result objects
#'
#' `[` on a `mixedcirc_fit_list` returns a single `mixedcirc_fit` when exactly
#' one feature is selected, preserving the package's original API. Selecting
#' multiple features returns another `mixedcirc_fit_list`. `[[` always extracts
#' one `mixedcirc_fit`. `[` on a `mixedcirc_fit` subsets its results data frame.
#'
#' @param x A `mixedcirc_fit_list` or `mixedcirc_fit` object.
#' @param i Row or feature indices.
#' @param j Column indices for a single-feature results data frame.
#' @param ... Additional subscripting arguments.
#' @param drop Logical; passed to data-frame subsetting for `mixedcirc_fit`.
#'
#' @return A `mixedcirc_fit_list`, `mixedcirc_fit`, or data frame, depending on
#'   the method used.
#' @name mixedcirc-subset
NULL

#' @rdname mixedcirc-subset
#' @export
methods::setMethod(
  "[",
  signature(x = "mixedcirc_fit_list"),
  function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) {
      i <- seq_along(x@results)
    }
    selected <- x@results[i]
    if (length(selected) == 1L) {
      return(selected[[1L]])
    }
    methods::new("mixedcirc_fit_list", results = selected)
  }
)

#' @rdname mixedcirc-subset
#' @export
methods::setMethod(
  "[[",
  signature(x = "mixedcirc_fit_list"),
  function(x, i, j, ...) {
    if (missing(i)) {
      stop("i must identify one fitted feature", call. = FALSE)
    }
    x@results[[i]]
  }
)

#' @rdname mixedcirc-subset
#' @export
methods::setMethod(
  "[",
  signature(x = "mixedcirc_fit"),
  function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) i <- seq_len(nrow(x@results))
    if (missing(j)) j <- seq_len(ncol(x@results))
    x@results[i, j, ..., drop = drop]
  }
)

#' Number of fitted features
#'
#' @param x A `mixedcirc_fit_list` object.
#'
#' @return An integer count.
#' @export
methods::setMethod("length", "mixedcirc_fit_list", function(x) {
  length(x@results)
})

#' Plot a fitted circadian feature
#'
#' @param x A `mixedcirc_fit` object.
#' @param y Missing; present for compatibility with the `plot` generic.
#' @param ... Arguments passed to [mixedcirc_fit_plot()].
#'
#' @return The ggplot object returned by [mixedcirc_fit_plot()].
#' @export
methods::setMethod(
  "plot",
  signature(x = "mixedcirc_fit", y = "missing"),
  function(x, y, ...) mixedcirc_fit_plot(x, ...)
)
