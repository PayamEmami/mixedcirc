#' Generate circadian rhythm
#'
#' This functions gets a mixedcirc_fit circadian rhythm fitted using mixedcirc_detect function and interpolate data points between two times.
#'
#' @param x A class of mixedcirc_fit
#' @param period The rhythm period. Default: 24
#' @param min_time Minimum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @param max_time Maximum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @param npoints Maximun number of data points to generate Default: 200
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#'int_data<-mixedcirc_interpolate(results)
#' @return
#' A data.frame with all the elements of the model including "Y.hat" that is the interpolated data
#'
#' @details
#' In case of RRBS, the data is assumed to be log2. The M-values are calculated as Mathylated-Unmethylated
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
#' @importFrom graphics plot
#' @export
mixedcirc_interpolate<- function(x, period=24,
                              min_time=NULL,
                              max_time=NULL,
                              npoints=200) {

  if(npoints<=0)
  {
    stop("npoints must be higher than 0!")
  }

  object <- x
  if(is.null(plot_title)) {
    plot_title<-rownames(object@results)
    if(is.null(plot_title)) { plot_title <- "Variable" }
  }
  type_of_analysis<-object@type
  fit <- object@fit
  pr_rows <- c()
  if(class(fit) == "lm") {
    pr_rows <- rownames(fit[[1]]@fit$model)
  } else {
    pr_rows <- rownames(fit@frame)
  }

  exp_design <- object@exp_design[rownames(object@exp_design)%in%pr_rows,]
  all_combs <- expand.grid(group=unique(object@exp_design[,"group"]))
  to_be_predited <- c()
  if(is.null(min_time)) { min_time <- min(exp_design$time) }
  if(is.null(max_time)) { max_time <- max(exp_design$time) }
  replicate_id<-NULL
  if(type_of_analysis=="RRBS")
  {
    set.seed(1)
    replicate_id<-sample(fit@frame$replicate_id,size = 1)
  }
  for(i in 1:nrow(all_combs)) {
    gr<-as.character(all_combs[i,"group"])

    timeax <- seq(min_time, max(c(period,max_time)), length.out = npoints)
    newdata <- data.frame(time = timeax,
                          inphase = cos(2 * pi * timeax/period),
                          outphase = sin(2 * pi * timeax/period),
                          group=gr)

    if(type_of_analysis=="RRBS")
    {
      newdata_1<-cbind(newdata,scaler=1,replicate_id=as.character(replicate_id))
      newdata_2<-cbind(newdata,scaler=0,replicate_id=as.character(replicate_id))
      to_be_predited <- rbind(to_be_predited,rbind(newdata_1,newdata_2))
    }else{
      to_be_predited <- rbind(to_be_predited,newdata)
    }

  }

  to_be_predited$Y.hat <- predict(object = fit,
                                  newdata = to_be_predited,
                                  re.form=NA,
                                  level = 0)

  if(type_of_analysis=="RRBS")
  {
    to_be_predited_2<-to_be_predited
    to_be_predited<-to_be_predited[to_be_predited$scaler==1,]
    to_be_predited$Y.hat<-to_be_predited_2$Y.hat[to_be_predited_2$scaler==1]-to_be_predited_2$Y.hat[to_be_predited_2$scaler==0]
  }

  return(to_be_predited)
}
