#' Compares circadian rhythm
#'
#' This functions compare two circadian rhythm fitted using mixedcirc_detect function using cross-correlation.
#'
#' @param x A class of mixedcirc_fit
#' @param y A class of mixedcirc_fit
#' @param type It can be one of "original", "fitted", "simulate". original: original data will used. fitted values will be used. simulate: data will be simulated based on the models and will be used to do cross correlation. Default: original
#' @param merge_groups If TRUE, for multi-group experiments, groups are NOT correlated separately. Default: FALSE
#' @param lag.max maximum lag at which to calculate the cross-correlation. Will be automatically limited as in ccf. Default: NULL
#' @param level confidence level, from 0 to 1. Default is 0.95, that is, 95% confidence.
#' @B number of bootstrap simulations to obtain empirical critical values. Default is 1000.
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
#' A statistical components of cross correlation.
#'
#' @details
#' In case of RRBS, the data is assumed to be log2. The M-values are calculated as Mathylated-Unmethylated.
#' In any case, the data in x and y will be matched using all the columns in the exp_design (except the measurement)
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
mixedcirc_compare<- function(x,y,type=c( "original", "fitted", "simulate")[1],merge_groups=FALSE, period=24,
                             lag.max=NULL,
                             level = 0.95,
                             B = 1000,
                             xlab="x",ylab="y",
                              min_time=NULL,
                              max_time=NULL) {

  if(missing(x))
    stop("x has to be provided")


  if(missing(y))
    stop("y has to be provided")


  if(!is(x,"mixedcirc_fit"))
    stop("x has to be a mixedcirc_fit object")

  if(!is(y,"mixedcirc_fit"))
    stop("y has to be a mixedcirc_fit object")

  if(is.null(min_time)) { min_time <- min(exp_design$time) }
  if(is.null(max_time)) { max_time <- max(exp_design$time) }

  object <- x

  object2 <- y

  data_for_cff<-NULL


  type_cff <- match.arg(type,c("original","fitted","simulate"))


  if(type_cff%in%c("original","fitted"))
  {
    type_of_analysis<-object@type
    fit<-object@fit

    pr_rows <- c()
    if(class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    } else {
      pr_rows <- rownames(fit@frame)
    }

    exp_design <- object@exp_design[rownames(object@exp_design)%in%pr_rows,]
    if(class(fit) == "lm" & type_cff=="original") {
      raw_data<-data.frame(measure=fit$model[,"measure"],exp_design)
    } else {
      raw_data<-data.frame(measure=fit@frame[,"measure"],exp_design)
    }

    if(class(fit) == "lm"& type_cff=="fitted") {
      raw_data<-data.frame(measure=fitted(fit),exp_design)
    } else {
      raw_data<-data.frame(measure=fitted(fit),exp_design)
    }


    if(type_of_analysis=="RRBS")
    {
      raw_data_2<-raw_data
      raw_data<-raw_data[raw_data$scaler==1,]
      raw_data$measure<-raw_data_2[raw_data_2$scaler==1,]$measure-raw_data_2[raw_data_2$scaler==0,]$measure
    }

    raw_data_x<-raw_data

    ###### For y

    type_of_analysis<-object2@type
    fit<-object2@fit

    pr_rows <- c()
    if(class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    } else {
      pr_rows <- rownames(fit@frame)
    }

    exp_design <- object2@exp_design[rownames(object2@exp_design)%in%pr_rows,]
    if(class(fit) == "lm" & type_cff=="original") {
      raw_data<-data.frame(measure=fit$model[,"measure"],exp_design)
    } else if( type_cff=="original"){
      raw_data<-data.frame(measure=fit@frame[,"measure"],exp_design)
    }

    if(class(fit) == "lm"& type_cff=="fitted") {
      raw_data<-data.frame(measure=fitted(fit),exp_design)
    } else if( type_cff=="fitted"){
      raw_data<-data.frame(measure=fitted(fit),exp_design)
    }


    if(type_of_analysis=="RRBS")
    {
      raw_data_2<-raw_data
      raw_data<-raw_data[raw_data$scaler==1,]
      raw_data$measure<-raw_data_2[raw_data_2$scaler==1,]$measure-raw_data_2[raw_data_2$scaler==0,]$measure
    }

    raw_data_y<-raw_data
    data_for_cff<-merge.data.frame(raw_data_x[,c("time" ,"group","rep","measure")],raw_data_y[,c("time" ,"group","rep","measure")],by = c("time","group","rep"))
  }



  if(type_cff%in%c("simulate"))
  {
    type_of_analysis<-object@type
    fit<-object@fit

    pr_rows <- c()
    if(class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    } else {
      pr_rows <- rownames(fit@frame)
    }

    exp_design <- object@exp_design[rownames(object@exp_design)%in%pr_rows,]

    to_be_predited <- c()
    if(is.null(min_time)) { min_time <- min(exp_design$time) }
    if(is.null(max_time)) { max_time <- max(exp_design$time) }
    replicate_id<-NULL
    if(type_of_analysis=="RRBS")
    {
      set.seed(1)
      replicate_id<-sample(fit@frame$replicate_id,size = 1)
    }
    all_combs <- expand.grid(group=unique(object@exp_design[,"group"]))

    for(i in 1:nrow(all_combs)) {
      gr<-as.character(all_combs[i,"group"])

      timeax <- seq(min_time, max(c(period,max_time)), length.out = 200)
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
    to_be_predited$measure<-to_be_predited$Y.hat

    raw_data_x<-to_be_predited



    ###### For y

    type_of_analysis<-object2@type
    fit<-object2@fit

    pr_rows <- c()
    if(class(fit) == "lm") {
      pr_rows <- rownames(fit$model)
    } else {
      pr_rows <- rownames(fit@frame)
    }

    exp_design <- object2@exp_design[rownames(object2@exp_design)%in%pr_rows,]

    to_be_predited <- c()
    if(is.null(min_time)) { min_time <- min(exp_design$time) }
    if(is.null(max_time)) { max_time <- max(exp_design$time) }
    replicate_id<-NULL
    if(type_of_analysis=="RRBS")
    {
      set.seed(1)
      replicate_id<-sample(fit@frame$replicate_id,size = 1)
    }
    all_combs <- expand.grid(group=unique(object@exp_design[,"group"]))

    for(i in 1:nrow(all_combs)) {
      gr<-as.character(all_combs[i,"group"])

      timeax <- seq(min_time, max(c(period,max_time)), length.out = 200)
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
    to_be_predited$measure<-to_be_predited$Y.hat

    raw_data_y<-to_be_predited

    data_for_cff<-merge.data.frame(raw_data_x[,c("time" ,"group","measure")],raw_data_y[,c("time" ,"group","measure")],by = c("time","group"))
  }

  if(length(unique(data_for_cff$group))==1)
    merge_groups <- TRUE
  list_output_values<-list()
if(merge_groups==TRUE)
{
  output<-funtimes::ccf_boot(data_for_cff[,"measure.x"],data_for_cff[,"measure.y"],lag.max = lag.max,level = level,B = B,plot = "none")
  list_output_values[["merged_information"]]<-output
  sig_groups<-rep("N.S",nrow(output))
  sig_groups[output$pS<(1-level)]<-paste0("<",1-level)
  plt<-ggplot(output)+geom_point(aes(Lag,rS,color=sig_groups,fill=sig_groups))+
    geom_segment(aes(x = Lag, y = 0, xend = Lag, yend = rS))+
    geom_ribbon(aes(ymin=lowerS, ymax=upperS,x=Lag,y=rS), linetype=2, alpha=0.1)+theme_classic2()+scale_color_aaas()+xlab("Lag")+
    ylab(paste0("Spearman", " correlation of ", xlab, "(t + Lag)",
                " and ", ylab, "(t)\n", "with ", level * 100,
                "% bootstrapped confidence region"))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))

  list_output_values[["plot"]]<-plt
  return(list_output_values)
}else{


  for(gr in unique(data_for_cff$group))
  {
    output<-funtimes::ccf_boot(data_for_cff[data_for_cff$group==gr,"measure.x"],data_for_cff[data_for_cff$group==gr,"measure.y"],lag.max = lag.max,level = level,B = B,plot = "none")


    list_output_values[[paste0(gr,"_information")]]<-output
    sig_groups<-rep("N.S",nrow(output))
    sig_groups[output$pS<(1-level)]<-paste0("<",1-level)

    assign(x = paste0("var_",gr,"_output"),value = cbind(output,sig_groups=sig_groups))
    list_output_values[[paste0(gr,"_plot")]]<- ggplot(get(paste0("var_",gr,"_output")))+geom_point(aes(Lag,rS,color=sig_groups,fill=sig_groups))+
      geom_segment(aes(x = Lag, y = 0, xend = Lag, yend = rS))+
      geom_ribbon(aes(ymin=lowerS, ymax=upperS,x=Lag,y=rS), linetype=2, alpha=0.1)+theme_classic2()+scale_color_aaas()+xlab("Lag")+
      ylab(paste0("Spearman", " correlation of ", xlab, "(t + Lag)",
                  " and ", ylab, "(t)\n", "with ", level * 100,
                  "% bootstrapped confidence region"))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))

  }
  return(list_output_values)
}

}
