#' Extract the fitted circadian rhythm
#'
#' This function extract the fitted values of circadian rhythm.
#'
#' @param fit an object of class mixedcirc_fit or mixedcirc_fit_list
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @export
#' @examples
#' library(mixedcirc)
#' data("circa_data")
#' results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'                          time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#' detrended<-mixedcirc_fitted(results)
#'
#'
#' @return
#' A matrix of detrended data
#'
#' @details
#'
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
#' @import mixOmics
#' @import dplyr

mixedcirc_fitted<-function (fit = NULL, per_group = FALSE, verbose = FALSE, ...)
{
  if (verbose)
    cat("Checking inputs ...\n")
  if (!(is(fit, "mixedcirc_fit") | is(fit, "mixedcirc_fit_list")))
    stop("input must be mixedcirc_fit_list or mixedcirc_fit!")
  if (verbose)
    cat("Performing detrending ...\n")
  if (is(fit, "mixedcirc_fit_list")) {
    fitted_values <- matrix(NA, nrow = nrow(fit[1]@exp_design),
                            ncol = length(fit))
    gr_pos <- rep("gr", nrow(fit[1]@exp_design))
    ft_names <- c()
    for (i in 1:length(fit)) {
      ft_names <- c(ft_names, rownames(fit[i]@results))
      if (verbose)
        cat("Performing detrending on", rownames(fit[i]@results),
            "...\n")
      for (grs in unique(gr_pos)) {
        object <- fit[i]
        fit_obj <- object@fit
        gr_pos_l<-gr_pos[rownames(object@exp_design)%in%rownames(object@fit@frame)]
        indexe <- which(gr_pos_l == grs)

        fitted_v <- fitted(fit_obj)[indexe]
        fitted_values[which(gr_pos == grs)[rownames(object@exp_design[gr_pos == grs,])%in%rownames(object@fit@frame)], i] <- fitted_v
      }
    }
    colnames(fitted_values) <- ft_names
  }
  else {
    fitted_values <- matrix(NA, nrow = nrow(fit@exp_design),
                            ncol = 1)
    gr_pos <- rep("gr", nrow(fit@exp_design))
    ft_names <- c()
    if (verbose)
      cat("Performing detrending on", rownames(fit@results),
          "...\n")
    for (grs in unique(gr_pos)) {

      fit_obj <- fit@fit
      gr_pos_l<-gr_pos[rownames(fit@exp_design)%in%rownames(fit@fit@frame)]
      indexe <- which(gr_pos_l == grs)

      fitted_v <- fitted(fit_obj)[indexe]

      fitted_values[which(gr_pos == grs)[rownames(fit@exp_design[gr_pos == grs,])%in%rownames(fit@fit@frame)]] <- as.matrix(fitted_v)
    }
    colnames(fitted_values) <- rownames(fit@results)
  }
  return(fitted_values)
}
