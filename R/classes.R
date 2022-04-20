#' Represents the result of model fitting on a single variable.
#'
#' @slot results A one row data frame that contains statistics of fitted model
#' @slot fit is the actual fitted model
#' @slot exp_design is the experimental design file created by mixedcirc_detect
#'
mixedcirc_fit <- setClass("mixedcirc_fit",
                          slots = list(results = "data.frame",fit="ANY",exp_design="data.frame")
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

print.mixedcirc_fit<- function(x) {print(x@results)}
setMethod(show, "mixedcirc_fit", function(object) {print(x@results)})



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

print.mixedcirc_fit_list<-function(x,...) {cat("Total number of variables: ",length(x@results))}
print.mixedcirc_fit<-function(x,...) {print(x@results,...)}

setMethod(show, "mixedcirc_fit", function(object) {print(object@results)})
setMethod(show, "mixedcirc_fit_list", function(object) {cat("Total number of variables: ",length(object@results))})


setMethod(`[`, "mixedcirc_fit_list", function(x,i=NULL,j=NULL) {
  x@results[[i]]
}
)

setMethod(`[`, "mixedcirc_fit", function(x,i=NULL,j=NULL) {
  x@results[j]
}
)
setMethod(length, "mixedcirc_fit_list", function(x) {length(x@results)})

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
#' @export
setMethod("plot", signature = c(x="mixedcirc_fit",y="missing"), function(x,y,...) {
  mixedcirc_fit_plot(x,...)
})
