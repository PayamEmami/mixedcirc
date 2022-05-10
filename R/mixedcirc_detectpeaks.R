#' Finds peak of a rhythm
#'
#' This functions finds peak of a rhythm
#'
#' @param phi phase of a peak
#' @param period Period of circadian rhythm. Default: 24
#' @param min_time Minimum time considered
#' @param correct_period If TRUE, period will be forced to be between the min and max time
#' @param correct_min minimum time will be adjusted
#' @export
#' @examples
#' data("circa_data")
#'
#'results <- mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,
#'group = circa_data$group,
#'id = circa_data$id,
#'period = 24,
#'verbose = TRUE)
#'mixedcirc_detectpeaks(results[1][,5], min_time=0)
#'
#' @return
#' the peak time
#'
#' @details
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
#' @import variancePartition
#' @import multtest


mixedcirc_detectpeaks <- function(phi,
                                  period = 24,
                                  min_time = 0,
                                  correct_period = T,
                                  correct_min = F,) {
  if(!is.numeric(phi))
    stop("phi must be numeric")

  if(!is.numeric(period))
    stop("period must be numeric")

  if(!is.numeric(min_time))
    stop("min_time must be numeric")

  if(!is.logical(correct_period))
    stop("correct_period must be TRUE or FALSE")

  if(!is.logical(correct_min))
    stop("correct_min must be TRUE or FALSE")

  tau <- period
  peak_time <- phi * tau / (2 * pi)
  #peak_time<-round(peak_time)
  if(correct_period) {
    while (peak_time > tau | peak_time < 0){
      if (peak_time > tau) {
        peak_time <- peak_time - tau
      }
      if (peak_time < 0) {
        peak_time <- peak_time + tau
      }
    }
    if (correct_min & peak_time < min_time) {
      peak_time <- peak_time + tau
    }
  }
  return(peak_time)
}
