#' Represents the result of model fitting on a single variable.
#'
#' @slot results A one row data frame that contains statistics of fitted model
#' @slot fit is the actual fitted model
#' @slot exp_design is the experimental design file created by mixedcirc_detect
#'
mixedcirc_fit <- setClass("mixedcirc_fit",
                          slots = list(results = "data.frame",
                                       fit = "ANY",
                                       exp_design = "data.frame")
)

#' Represents the result of model fitting on several variables.
#'
#' @slot results A list of mixedcirc_fit
#'
mixedcirc_fit_list <- setClass("mixedcirc_fit_list",
                          slots = list(results = "list")
)

#' Plots circadian rhythm
#'
#' A simple function for printing the mixedcirc_fit object
#'
#'
#' @import stats
#' @import multcomp
#' @import doFuture
#' @import future
#' @import nlme
#' @import future.apply
#' @import lme4
#' @import limma
#' @import lmerTest
#' @import foreach
#' @import ggplot2
#' @import ggsci

print.mixedcirc_fit <- function(x) { print(x@results) }
setMethod(f = "show",
          signature = "mixedcirc_fit",
          definition = function(object) { print(x@results) } )



#' Prints circadian rhythm
#'
#' A simple function for printing the mixedcirc_fit object
#'
#'
#' @import stats
#' @import multcomp
#' @import doFuture
#' @import future
#' @import nlme
#' @import future.apply
#' @import lme4
#' @import limma
#' @import lmerTest
#' @import foreach
#' @import ggplot2
#' @import ggsci

setMethod(f = "show",
          signature = "mixedcirc_fit",
          definition = function(object) { print(object@results) }
)

print.mixedcirc_fit_list <- function(x, ...) {
	cat('Total number of variables: ', length(x@results))
}

print.mixedcirc_fit <- function(x, ...) {
	print(x@results, ...)
}

setMethod(show, "mixedcirc_integration", function(object) {
  cat('Total number of Omics: ', length(object@partial), '! Use partial or average scores in mixedcirc_detect!')
})

print.mixedcirc_integration <- function(x, ...) {
	cat('Total number of Omics: ', length(object@partial),'! Use partial or average scores in mixedcirc_detect!')
}

setMethod(f = "show",
          signature = "mixedcirc_fit_list",
          definition = function(object) { cat("Total number of variables: ", length(object@results)) }
)

setMethod(f = `[`,
          signature = "mixedcirc_fit_list",
          definition = function(x, i = NULL, j = NULL) { x@results[[i]] }
)

setMethod(f = `[`,
          signature = "mixedcirc_fit",
          definition = function(x, i = NULL, j = NULL) { x@results[j] }
)

setMethod(f = "length",
          signature = "mixedcirc_fit_list",
          definition = function(x) { length(x@results) }
)

setGeneric(name = "plot",
           def = function(x, y, ...) standardGeneric("plot")
)

#' @export
setMethod(f = "plot",
          signature = c(x = "mixedcirc_fit", y = "missing"),
          definition = function(x, y, ...){ mixedcirc_fit_plot(x, ...) }
         )

#' Represents the result of integration on multiple data sets
#'
#' @slot partial A list of partial scores for individual data sets
#' @slot average A matrix of average scores
#' @slot loadings A List of loadings for individual data sets
#' @slot model The original model fitted
#'
mixedcirc_integration <- setClass('mixedcirc_integration',
	slots = list(partial = 'list',
		average = 'matrix',
		loadings = 'list', model = 'ANY')
)
