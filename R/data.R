#' Example longitudinal circadian data
#'
#' A simulated data set for demonstrating grouped, repeated-measures circadian
#' analyses. It contains 20 subjects observed hourly over 72 hours in each of
#' two groups, with eight simulated response variables.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{data_matrix}{A numeric matrix with 2,880 observations in rows and
#'     eight response variables in columns.}
#'   \item{time}{An integer vector giving time from 1 to 72 hours.}
#'   \item{group}{A character vector with levels `"g1"` and `"g2"`.}
#'   \item{id}{A character vector identifying the 20 repeatedly measured
#'     subjects.}
#' }
#'
#' @source Simulated with mixedcirc for package examples.
#' @usage data(circa_data)
"circa_data"
