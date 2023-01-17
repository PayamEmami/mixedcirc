#' Removes linear trend from the fitted circadian rhythm
#'
#' This function removes linear trend from the fitted values of circadian rhythm.
#'
#' @param fit an object of class mixedcirc_fit or mixedcirc_fit_list
#' @param per_group If TRUE, detrending is performed per group (default: FALSE)
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @param ... additionl arguments to the regression function
#' @export
#' @examples
#' library(mixedcirc)
#' data("circa_data")
#' results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'                          time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#' detrended<-mixedcirc_detrend(results)
#'
#'
#' @return
#' A matrix of detrended data
#'
#' @details
#' This method will first calculate the fitted values of the provided model followed by fitting a regression of form x~time. The residues of this regression model will be returned as the detrended data.
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


############ mixedcirc_detrend must be changed so it will only output the fitted value. Nothing else!!
############ mixedcirc_detrend must be changed so it will only output the fitted value. Nothing else!!

############ mixedcirc_detrend must be changed so it will only output the fitted value. Nothing else!!
mixedcirc_fitted<-function (fit = NULL, per_group = FALSE, verbose = FALSE, ...)
{
  if (verbose)
    cat("Checking inputs ...\n")
  if (!(is(fit, "mixedcirc_fit") | is(fit, "mixedcirc_fit_list")))
    stop("input must be mixedcirc_fit_list or mixedcirc_fit!")
  if (verbose)
    cat("Performing detrending ...\n")
  if (verbose)
    if (per_group)
      cat("Performing regression on each group sepearately ...\n")
  if (is(fit, "mixedcirc_fit_list")) {
    fitted_values <- matrix(NA, nrow = nrow(fit[1]@exp_design),
                            ncol = length(fit))
    gr_pos <- rep("gr", nrow(fit[1]@exp_design))
    if ("group" %in% colnames(fit[1]@exp_design) & per_group) {
      gr_pos <- fit[1]@exp_design[, "group"]
    }
    else if (per_group) {
      warning("There is no group information in the experimental design! Performing detrending without grouping!\n")
    }
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
        de_tr <- NA
        if (class(fit_obj) == "lm") {
          de_tr <- residuals(lm(fitted_v ~ object@exp_design$time[indexe],
                                weights = fit_obj@weights[indexe]))
        }
        else {
          weights <- NULL
          if ("(weights)" %in% colnames(fit_obj@frame))
            weights <- fit_obj@frame[indexe, "(weights)"]
          de_tr <- residuals(lme4::lmer(measure ~ time +
                                          (1 | rep), data = cbind.data.frame(measure = fitted_v,
                                                                             time = object@exp_design$time[indexe], rep = as.character(object@exp_design$rep)[indexe]),
                                        weights = weights))
        }
        fitted_values[which(gr_pos == grs)[rownames(object@exp_design[gr_pos == grs,])%in%rownames(object@fit@frame)], i] <- de_tr
      }
    }
    colnames(fitted_values) <- ft_names
  }
  else {
    fitted_values <- matrix(NA, nrow = nrow(fit@exp_design),
                            ncol = 1)
    gr_pos <- rep("gr", nrow(fit@exp_design))
    if ("group" %in% colnames(fit@exp_design) & per_group) {
      gr_pos <- fit@exp_design[, "group"]
    }
    else if (per_group) {
      warning("There is no group information in the experimental design! Performing detrending without grouping!\n")
    }
    ft_names <- c()
    if (verbose)
      cat("Performing detrending on", rownames(fit@results),
          "...\n")
    for (grs in unique(gr_pos)) {

      fit_obj <- fit@fit
      gr_pos_l<-gr_pos[rownames(fit@exp_design)%in%rownames(fit@fit@frame)]
      indexe <- which(gr_pos_l == grs)

      fitted_v <- fitted(fit_obj)[indexe]
      de_tr <- NA
      if (class(fit_obj) == "lm") {
        de_tr <- residuals(lm(fitted_v ~ fit@exp_design$time[indexe],
                              weights = fit_obj@weights[indexe], ...))
      }
      else {
        weights <- NULL
        if ("(weights)" %in% colnames(fit_obj@frame))
          weights <- fit_obj@frame[indexe, "(weights)"]
        de_tr <- residuals(lme4::lmer(measure ~ time +
                                        (1 | rep), data = cbind.data.frame(measure = fitted_v,
                                                                           time = fit@exp_design$time[indexe], rep = as.character(fit@exp_design$rep)[indexe]),
                                      weights = weights, ...))
      }
      fitted_values[which(gr_pos == grs)[rownames(fit@exp_design[gr_pos == grs,])%in%rownames(fit@fit@frame)]] <- as.matrix(de_tr)
    }
    colnames(fitted_values) <- rownames(fit@results)
  }
  return(fitted_values)
}
