#' Plots circadian rhythm
#'
#' This functions plots circadian rhythm fitted using mixedcirc_detect function.
#'
#' @param x A class of mixedcirc_fit
#' @param y Not used
#' @param period The rhythm period. Default: 24
#' @param min_time Minimum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @param max_time Maximum time span to do the prediction. If NULL, it will be taken from the fit. Default: NULL
#' @param plot_title Title of the plot.
#' @param plot_points whether to plot the individual points or not. Default TRUE
#' @param plot_smooth Whether to plot smoothing line or not. Default: FALSE
#' @param plot_trend Whether to plot trend lines or not. Default: TRUE
#' @param xlab Label on the x-axis
#' @param ylab Label on the y-axis
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#' plot(results[[1]])
#' @return
#' A ggplot object.
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
#' @importFrom graphics plot
#' @export
mixedcirc_fit_plot<- function(x,period=24,
                              min_time=NULL,max_time=NULL,plot_title=NULL,
                              plot_points=TRUE,plot_smooth=FALSE,plot_trend=TRUE,xlab="time",ylab="Expression") {
  object<-x
  if(is.null(plot_title))
  {
    plot_title<-rownames(object@results)
    if(is.null(plot_title))plot_title<-"Variable"
  }

  fit<-object@fit
  pr_rows<-c()
  if(class(fit)=="lm")
    pr_rows<-rownames(fit[[1]]@fit$model)
  else
    pr_rows<-rownames(fit@frame)

  exp_design<-object@exp_design[rownames(object@exp_design)%in%pr_rows,]
  all_combs<-expand.grid(group=unique(object@exp_design[,"group"]))
  to_be_predited<-c()
  if(is.null(min_time))
    min_time<-min(exp_design$time)
  if(is.null(max_time))
    max_time<-max(exp_design$time)
  for(i in 1:nrow(all_combs))
  {
    gr<-as.character(all_combs[i,"group"])

    timeax <- seq(min_time, max(c(period,max_time)), length.out = 200)
    newdata <- data.frame(time = timeax, inphase = cos(2 * pi * timeax/period),
                          outphase = sin(2 * pi * timeax/period),group=gr)

    to_be_predited<-rbind(to_be_predited,newdata)
  }

  to_be_predited$Y.hat<-predict(object = fit,newdata = to_be_predited,re.form=NA,level = 0)

  library(ggplot2)
  library(ggsci)

  if(class(fit)=="lm")
    raw_data<-data.frame(measure=results[[1]]@fit$model[,"measure"],exp_design)
  else
    raw_data<-data.frame(measure=fit@frame[,"measure"],exp_design)
  time<-raw_data$time

  plot_tmp<-ggplot(data=raw_data,aes(x=time,y=measure,group=group,color=group,shape=group))

  if(plot_points)
    plot_tmp <- plot_tmp+geom_point(data=raw_data,aes(x=time,y=measure,group=group,color=group,shape=group),position = position_dodge(2),alpha=0.2)

  if(plot_trend)
    plot_tmp<-plot_tmp+geom_line(data=to_be_predited, aes_string(x = "time", y = "Y.hat", col = "group"))

  if(plot_smooth)
    plot_tmp<-plot_tmp+geom_smooth(data=raw_data,aes(x=time,y=measure,group=group,color=group),linetype = "dashed",se = F,alpha=0.2)

  plot_tmp<-plot_tmp+theme_classic2()+
    scale_color_aaas()+scale_x_continuous( labels = unique(as.character(time)), breaks = unique(time))+ylab("Expression")+scale_y_continuous()+
    theme(axis.text = element_text(color = "black"))+ggtitle(plot_title)+xlab(xlab)+ylab(ylab)

  if(length(unique(object@exp_design[,"group"]))<2)
    plot_tmp<-plot_tmp+theme(legend.position="none")
  plot_tmp
}

