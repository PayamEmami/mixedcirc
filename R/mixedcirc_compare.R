#' Compares circadian rhythm
#'
#'This functions compare two circadian rhythm fitted using mixedcirc_detect function using standard or cross-correlation.
#'
#' @param x A class of mixedcirc_fit
#' @param y A class of mixedcirc_fit
#' @param cor_type It can be one of "standard" or "crosscorrelation". Default:crosscorrelation  See details.
#' @param type It can be one of "original", "fitted", "simulate". original: original data will used. fitted values will be used. simulate: data will be simulated based on the models and will be used to do cross correlation. Default: original
#' @param merge_groups If TRUE, for multi-group experiments, groups are NOT correlated separately. Default: FALSE
#' @param lag.max maximum lag at which to calculate the cross-correlation. Will be automatically limited as in ccf. Default: NULL
#' @param level confidence level, from 0 to 1. Default is 0.95, that is, 95% confidence.
#' @param number of bootstrap simulations to obtain empirical critical values. Default is 1000.
#' @param period The rhythm period. Default: 24
#' @param min_time Minimum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @param max_time Maximum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#' mixedcirc_compare(results[1],results[2])
#'
#' @return
#' A list containing:
#' A ggplot object.
#' A statistical components of standard and cross correlation.
#'
#' @details
#'
#' cor_type
#'
#' One can doe standard and cross correlation between two series.
#' In standard mode, we fit a model of y~x and calculate r-squared (marginal r-squared  in case of repeated measures).
#' The correlation is then simply assumed to be sign of coeffect of X * sqrt(r-squared)
#'
#' This is a typical cross-correlation. In case of repeated measures specially the p-values cannot be trusted.
#'
#'
#' In case of RRBS, the data is assumed to be log2. The M-values are calculated as Mathylated-Unmethylated.
#' In any case, the data in x and y will be matched using all the columns in the exp_design (except the measurement)
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
#' @import ggplot2
#' @import ggsci
#' @import ggpubr
#' @import funtimes
#' @importFrom graphics plot
#' @export
mixedcirc_compare<-function (x, y,cor_type = c("standard","crosscorrelation")[2], type = c("original", "fitted", "simulate")[1],
          merge_groups = FALSE, period = 24, lag.max = NULL, level = 0.95,
          B = 1000, xlab = "x", ylab = "y", min_time = NULL, max_time = NULL)
{


  if (missing(x))
    stop("x has to be provided")
  if (missing(y))
    stop("y has to be provided")
  if (!is(x, "mixedcirc_fit"))
    stop("x has to be a mixedcirc_fit object")
  if (!is(y, "mixedcirc_fit"))
    stop("y has to be a mixedcirc_fit object")

  object <- x
  object2 <- y

  if(class(object@fit)!=class(object2@fit))
  {
    stop("x and y must have been fitted using identical methods eg. lm or lmer!")
  }

  data_for_cff <- NULL
  type_cff <- match.arg(type, c("original", "fitted", "simulate"))
  type_cor <- match.arg(cor_type, c("standard","crosscorrelation"))


  if (type_cff %in% c("original", "fitted")) {
    type_of_analysis <- object@type
    fit <- object@fit
    pr_rows <- c()
    if (class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    }
    else {
      pr_rows <- rownames(fit@frame)
    }
    exp_design <- object@exp_design[rownames(object@exp_design) %in%
                                      pr_rows, ]
    if (class(fit) == "lm" & type_cff == "original") {
      raw_data <- data.frame(measure = fit$model[, "measure"],
                             exp_design)
    }else if(type_cff == "original"){
      raw_data <- data.frame(measure = fit@frame[, "measure"],
                             exp_design)
    }
    if (class(fit) == "lm" & type_cff == "fitted") {
      raw_data <- data.frame(measure = fitted(fit), exp_design)
    }else if(class(fit) == "lmerMod" & type_cff == "fitted"){
      raw_data <- data.frame(measure = fitted(fit), exp_design)
    }
    if (type_of_analysis == "RRBS") {
      raw_data_2 <- raw_data
      raw_data <- raw_data[raw_data$scaler == 1, ]
      raw_data$measure <- raw_data_2[raw_data_2$scaler ==
                                       1, ]$measure - raw_data_2[raw_data_2$scaler ==
                                                                   0, ]$measure
    }
    raw_data_x <- raw_data
    type_of_analysis <- object2@type
    fit <- object2@fit
    pr_rows <- c()
    if (class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    }
    else {
      pr_rows <- rownames(fit@frame)
    }
    exp_design <- object2@exp_design[rownames(object2@exp_design) %in%
                                       pr_rows, ]
    if (class(fit) == "lm" & type_cff == "original") {
      raw_data <- data.frame(measure = fit$model[, "measure"],
                             exp_design)
    }else if(type_cff == "original"){
      raw_data <- data.frame(measure = fit@frame[, "measure"],
                             exp_design)
    }
    if (class(fit) == "lm" & type_cff == "fitted") {
      raw_data <- data.frame(measure = fitted(fit), exp_design)
    }else if(class(fit) == "lmerMod" & type_cff == "fitted"){
      raw_data <- data.frame(measure = fitted(fit), exp_design)
    }
    if (type_of_analysis == "RRBS") {
      raw_data_2 <- raw_data
      raw_data <- raw_data[raw_data$scaler == 1, ]
      raw_data$measure <- raw_data_2[raw_data_2$scaler ==
                                       1, ]$measure - raw_data_2[raw_data_2$scaler ==
                                                                   0, ]$measure
    }
    raw_data_y <- raw_data
    data_for_cff <- merge.data.frame(raw_data_x[, c("time",
                                                    "group", "rep", "measure")], raw_data_y[, c("time",
                                                                                                "group", "rep", "measure")], by = c("time", "group",
                                                                                                                                    "rep"))
  }
  if (type_cff %in% c("simulate")) {
    suppressMessages({
      to_be_predited<-mixedcirc_interpolate(object,period = period,min_time = min_time,max_time = max_time)
    })
    to_be_predited$measure <- to_be_predited$Y.hat
    raw_data_x <- to_be_predited
    suppressMessages({
      to_be_predited<-mixedcirc_interpolate(object2,period = period,min_time = min_time,max_time = max_time)
    })
    to_be_predited$measure <- to_be_predited$Y.hat

    to_be_predited$measure <- to_be_predited$Y.hat
    raw_data_y <- to_be_predited
    data_for_cff <- merge.data.frame(raw_data_x[, c("time",
                                                    "group", "measure")], raw_data_y[, c("time", "group",
                                                                                         "measure")], by = c("time", "group"))
  }

  list_output_values<-list()
  if(type_cor=="standard")
  {
    if (length(unique(data_for_cff$group)) == 1)
      merge_groups <- TRUE


    if(class(fit)=="lm" | type_cff=="simulate")
    {
      if(merge_groups)
      {
        mdl<-summary(lm(measure.y~measure.x,data = data_for_cff))
        r_est<-sign(mdl$coefficients[2,"Estimate"])*sqrt(mdl$r.squared)
        output<-data.frame(Estimate=mdl$coefficients[2,"Estimate"],p.value=mdl$coefficients[2,"Pr(>|t|)"],r.squared=(mdl$r.squared),r=r_est)
        list_output_values[["merged_information"]] <- output
        return(list_output_values)

      }else{
        for (gr in unique(data_for_cff$group)) {



          mdl<-summary(lm(measure.y~measure.x,data = data_for_cff[data_for_cff$group==gr,]))
          r_est<-sign(mdl$coefficients[2,"Estimate"])*sqrt(mdl$r.squared)
          output<-data.frame(Estimate=mdl$coefficients[2,"Estimate"],p.value=mdl$coefficients[2,"Pr(>|t|)"],r.squared=(mdl$r.squared),r=r_est)


          list_output_values[[paste0(gr, "_information")]] <- output
        }

        return(list_output_values)
      }

    }else{

      if(merge_groups)
      {
        suppressMessages({
          lme_m<-lmerTest::lmer(measure.y~measure.x+(1|rep),data = data_for_cff, REML=F)
        })
        mdl<-summary(lme_m)

        fixed_eff <- lme4::fixef(lme_m)
        estm <- !is.na(fixed_eff)
        fittedvar <- var((model.matrix(lme_m)[, estm, drop = FALSE] %*% fixed_eff[estm])[,
                                                                                         1L])


        rsq<-fittedvar/(as.data.frame(VarCorr(lme_m))[1,"vcov"]+
                          as.data.frame(VarCorr(lme_m))[2,"vcov"]+fittedvar)
        r_es<-sign(mdl$coefficients[2,"Estimate"])*sqrt(rsq)
        output<-data.frame(Estimate=mdl$coefficients[2,"Estimate"],p.value=mdl$coefficients[2,"Pr(>|t|)"],r.squared=rsq,r=r_es)
        list_output_values[["merged_information"]] <- output
        return(list_output_values)

      }else{
        for (gr in unique(data_for_cff$group)) {



          suppressMessages({
            lme_m<-lmerTest::lmer(measure.y~measure.x+(1|rep),data = data_for_cff[data_for_cff$group==gr,], REML=F)
          })
          mdl<-summary(lme_m)
          fixed_eff <- lme4::fixef(lme_m)
          estm <- !is.na(fixed_eff)
          fittedvar <- var((model.matrix(lme_m)[, estm, drop = FALSE] %*% fixed_eff[estm])[,
                                                                                           1L])


          rsq<-fittedvar/(as.data.frame(VarCorr(lme_m))[1,"vcov"]+
                            as.data.frame(VarCorr(lme_m))[2,"vcov"]+fittedvar)
          r_es<-sign(mdl$coefficients[2,"Estimate"])*sqrt(rsq)
          output<-data.frame(Estimate=mdl$coefficients[2,"Estimate"],p.value=mdl$coefficients[2,"Pr(>|t|)"],r.squared=rsq,r=r_es)



          list_output_values[[paste0(gr, "_information")]] <- output
        }

        return(list_output_values)
      }
    }

  }else if(type_cor=="crosscorrelation")
  {
    if (length(unique(data_for_cff$group)) == 1)
      merge_groups <- TRUE
    list_output_values <- list()
    if (merge_groups == TRUE) {
      output <- funtimes::ccf_boot(data_for_cff[, "measure.x"],
                                   data_for_cff[, "measure.y"], lag.max = lag.max, level = level,
                                   B = B, plot = "none")
      list_output_values[["merged_information"]] <- output
      sig_groups <- rep("N.S", nrow(output))
      sig_groups[output$pS < (1 - level)] <- paste0("<", 1 -
                                                      level)
      plt <- ggplot(output) + geom_point(aes(Lag, rS, color = sig_groups,
                                             fill = sig_groups)) + geom_segment(aes(x = Lag, y = 0,
                                                                                    xend = Lag, yend = rS)) + geom_ribbon(aes(ymin = lowerS,
                                                                                                                              ymax = upperS, x = Lag, y = rS), linetype = 2, alpha = 0.1) +
        theme_classic2() + scale_color_aaas() + xlab("Lag") +
        ylab(paste0("Spearman", " correlation of ", xlab,
                    "(t + Lag)", " and ", ylab, "(t)\n", "with ",
                    level * 100, "% bootstrapped confidence region")) +
        guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
      list_output_values[["plot"]] <- plt
      return(list_output_values)
    }
    else {
      for (gr in unique(data_for_cff$group)) {
        output <- funtimes::ccf_boot(data_for_cff[data_for_cff$group ==
                                                    gr, "measure.x"], data_for_cff[data_for_cff$group ==
                                                                                     gr, "measure.y"], lag.max = lag.max, level = level,
                                     B = B, plot = "none")
        list_output_values[[paste0(gr, "_information")]] <- output
        sig_groups <- rep("N.S", nrow(output))
        sig_groups[output$pS < (1 - level)] <- paste0("<",
                                                      1 - level)
        assign(x = paste0("var_", gr, "_output"), value = cbind(output,
                                                                sig_groups = sig_groups))
        list_output_values[[paste0(gr, "_plot")]] <- ggplot(get(paste0("var_",
                                                                       gr, "_output"))) + geom_point(aes(Lag, rS, color = sig_groups,
                                                                                                         fill = sig_groups)) + geom_segment(aes(x = Lag,
                                                                                                                                                y = 0, xend = Lag, yend = rS)) + geom_ribbon(aes(ymin = lowerS,
                                                                                                                                                                                                 ymax = upperS, x = Lag, y = rS), linetype = 2,
                                                                                                                                                                                             alpha = 0.1) + theme_classic2() + scale_color_aaas() +
          xlab("Lag") + ylab(paste0("Spearman", " correlation of ",
                                    xlab, "(t + Lag)", " and ", ylab, "(t)\n", "with ",
                                    level * 100, "% bootstrapped confidence region")) +
          guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
      }
      return(list_output_values)
    }
  }


  return(list_output_values)
}
