#' Removes linear trend from the data
#'
#' This functions performs trend removing using mixed models.
#'
#' @param data_input A numerical matrix or data.frame (N*P) or DGEList where in the rows are samples (N) and the columns are variables (P). If DGEList is provided, regression weights will be estimated!
#' @param time A vector of length N, showing circadian time of each sample
#' @param group A character vector of length N. If performing differential circadian rhythm analysis, group is a factor, showing grouping of the samples. Analysis of two groups is supported at this stage! See details!
#' @param id A vector of length N showing identity of each *unique* sample. See details
#' @param lm_method The regression method to use. At this stage, `lm` and `lme` are supported! If lm is selected, normal regression will be performed. Default: "lme"
#' @param obs_weights Regression weights. Default: NULL. See details
#' @param RRBS If `TRUE`, the data is assumed to be RRBS methylation data. if TRUE, obs_weights must be set. Default FALSE
#' @param replicate_id If `RRBS` is set to `TRUE`, This has to be a factor showing identity of each unique replicate.
#' @param remove_trend_separate_groups If TRUE, the detrending is performed separately on each group (default:TRUE)
#' @param force_weight_estimation If TRUE, variance-mean trend weight estimation will be performed regardless of the data input type (default: FALSE)
#' @param ncores number of cores
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @param ... additionl arguments to the regression function
#' @export
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_remove_trend(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,verbose = TRUE)
#'
#' @return
#' A data.frame
#'
#' @details
#' For each variable we use the following mode:
#' In this part we do rhythmicity analysis on individual variables using the following model:
#' measure ~0 + time
#' The residulas will be outputed
#'
#' `obs_weights` is a matrix of size N*P where each colum shows the weights for all the observations for that particular variable.
#'
#' If `RRBS` is set to `TRUE`, we assume that the data is RRBS methylation. In this case, the regression will be change to suit this type of analysis.
#' In BS-seq methylation analysis, each DNA sample generates two counts, a count of methylated reads and a count of unmethylated reads, for each genomic locus for each sample.
#' The samples are assumed to be ordered as methylated and then unmethylated. For example, given the samples are A1_1, A1_2,A2_1, and B1_1.
#' The rows in `data_input` are assumed to be ordered as  A1_1_methylated,A1_1_unmethylated, A1_2_methylated,A1_2_unmethylated,A2_1_methylated,A2_1_unmethylated,B1_1_methylated,B1_1_unmethylated.
#' In this setting, `replicate_id` must show the unique identify of each replicate (in contrast to unique biological sample). This means for example above,
#' `replicate_id` would be A1_1,A1_2,A2_1,B1_1.
#' Please note that if `RRBS` is set to `TRUE`, `obs_weights` must be set. We provide a starting function to estimate these weights
#' (\code{\link{mixedcirc_rrbs_voom}}) for class of `methylBaseDB` from `methylKit` package. Alternatively, one can use `voomWithDreamWeights` with correct formula.
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
#' @import parallel

mixedcirc_remove_trend <- function(data_input=NULL,time=NULL,group=NULL,id=NULL,
                                   lm_method=c("lm","lme")[2],
                                   obs_weights=NULL,RRBS=FALSE,
                                   replicate_id=NULL,remove_trend_separate_groups=TRUE,force_weight_estimation=FALSE,ncores=1,verbose=FALSE,...){



  registerDoFuture()
  if(ncores>1)
  {
    plan(multisession,workers = ncores)
  }else{
    plan(sequential)
  }

  if(RRBS==TRUE)
  {
    if(is.null(obs_weights))
    {
      stop("For RRBS data, obs_weights must be set")
    }
    if(is.null(replicate_id))
    {
      stop("For RRBS data, replicate_id must be set")
    }

  }

  if(verbose)cat("Checking inputs ...\n")
  # checking inputs
  if(is.null(data_input))
    stop("data_input must be a data frame or matrix")

  if((!is.matrix(data_input) | !is.matrix(data_input)) & !is(data_input,"DGEList"))
    stop("data_input must be a data frame or matrix")

  if(is.null(time))
    stop("time must be a vector")

  if(RRBS==TRUE)
  {
    if(length(time)!=length(replicate_id))
    {
      stop("time must have the same length as replicate_id")
    }
    time <- rep(time,each=2)
  }



  if(!is(data_input,"DGEList"))
  {
    if(length(time)!=nrow(data_input))
      stop("The length of *time* is not equal to the number of rows of *data_input*")
  }else{
    if(length(time)!=ncol(data_input))
      stop("The length of *time* is not equal to the number of rows of *data_input*")
  }


  if(!is.null(group))
  {

    if(RRBS==TRUE)
    {
      if(length(group)!=length(replicate_id))
      {
        stop("time must have the same length as replicate_id")
      }
      group <- rep(group,each=2)
    }


    if(!is(data_input,"DGEList"))
    {
      if(length(group)!=nrow(data_input))
        stop("The length of *group* is not equal to the number of rows of *data_input*")
    }else{
      if(length(group)!=ncol(data_input))
        stop("The length of *group* is not equal to the number of rows of *data_input*")
    }



    if(length(unique(group))>2)
      stop("We are only supporting two groups at this stage!")

    multiple_groups <- TRUE
  }else{
    warning("No group information was provided! Skipping differential analysis!")
    multiple_groups <- FALSE
    group<-"dummy"
  }

  if(!is.null(id))
  {
    if(RRBS==TRUE)
    {
      if(length(id)!=length(replicate_id))
      {
        stop("time must have the same length as replicate_id")
      }
      id <- rep(id,each=2)
    }

    if(!is(data_input,"DGEList"))
    {
      if(length(id)!=nrow(data_input))
        stop("The length of *id* is not equal to the number of rows of *data_input*")
    }else{
      if(length(id)!=ncol(data_input))
        stop("The length of *id* is not equal to the number of rows of *data_input*")
    }


  }else{
    warning("No id information was provided! Switching to normal linear regression!")
    lm_method <- "lm"
    id<-"dummy"
  }


  lm_method <- match.arg(lm_method,c("lm","lme"))

  if(!is(data_input,"DGEList"))
    eset <- data_input[,,drop=F]
  else
    eset <- NA


  if(verbose)cat("Building experiment ...\n")
  if(RRBS==TRUE)
  {
    exp_design<-cbind.data.frame(time=as.numeric(time),group=as.factor(group),rep=as.character(id),
                                 replicate_id=rep(replicate_id,each=2),scaler=rep(c(1,0),(length(id)/2)))
  }else{
    exp_design<-cbind.data.frame(time=as.numeric(time),group=as.factor(group),rep=as.character(id))
  }


  if(!is.factor(exp_design$group))stop("Group must be a factor!")


  if(!is.null(obs_weights) & !is(data_input,"DGEList"))
  {
    if(!all(dim(obs_weights)==dim(data_input)))
      stop("obs_weights must be the same size as data_input!")
  }


  if(multiple_groups==TRUE)
  {
    if(verbose)cat("Performing multigroup detrending ...\n")

    group_id <- base::levels(exp_design$group)


    if(RRBS==TRUE)
    {

      design <- stats::model.matrix(~0+ replicate_id+ time:scaler ,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0+replicate_id+  time:scaler,
                                      data = exp_design)

    }else{

      design <- stats::model.matrix(~0+ time,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0+  time,
                                      data = exp_design)
    }



    design_t<-design

    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design) <- gsub(":", "_", colnames(design))

    if(is(data_input,"DGEList")| force_weight_estimation==TRUE){
      if(verbose)cat("Estimating weights for the input variables ...\n")

      if(lm_method=="lm")
      {
        if(verbose)cat("lm_method is lm, voom will be used ...\n")
        voom_res<-limma::voom(data_input, design, plot=FALSE)
      }else{
        if(RRBS==TRUE)
        {

          formula<-~0+ replicate_id + time:scaler +(1 | rep)
        }else{
          formula<-~0+  time +(1 | rep)
        }

        voom_res<-variancePartition::voomWithDreamWeights(data_input,formula = formula,data = exp_design,BPPARAM = BiocParallel::SnowParam(workers = ncores))
      }
      obs_weights <-t(voom_res$weights)
      eset<-t(voom_res$E)
    }


    chunks <- parallel::splitIndices(ncol(eset), min(ncol(eset),  ncores))
    if(verbose)cat("Spliting data to",length(chunks)," chunks...\n")

    res<-foreach(chIndx=chunks)%dopar%
      {

        outputs_fn<-foreach(i=chIndx)%do%{


          if(verbose)cat("Processing variable",i,"...\n")

          data_grouped<-cbind(measure=as.numeric(eset[,i]),exp_design)
          if(!remove_trend_separate_groups){


            if(verbose)cat("Fitting the model for variable",i,"...\n")
            if(RRBS==TRUE)
            {

              suppressMessages({
                model_ln<-switch(lm_method,
                                 lm = lm(measure ~0+replicate_id+ time:scaler ,data=data_grouped,weights = obs_weights[,i],...),
                                 lme = lme4::lmer(measure~0+ replicate_id+ time:scaler +(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
              })
            }else{
              suppressMessages({
                model_ln<-switch(lm_method,
                                 lm = lm(measure ~0+ time ,data=data_grouped,weights = obs_weights[,i],...),
                                 lme = lme4::lmer(measure~0+  time +(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
              })
            }


            data_grouped$measure[!is.na( data_grouped$measure)]<-residuals(model_ln)
          }else{
            # single_rhythm A
            if(RRBS==TRUE)
            {
              suppressMessages({
                model_ln_A<-switch(lm_method,
                                   lm = lm(measure ~0+replicate_id+ time:scaler ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...),
                                   lme = lmerTest::lmer(measure~0+ replicate_id+ time:scaler +(1 | rep) ,
                                                        control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),
                                                        data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...))
              })
            }else{
              suppressMessages({
                model_ln_A<-switch(lm_method,
                                   lm = lm(measure ~0+ time  ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...),
                                   lme = lmerTest::lmer(measure~0+  time +(1 | rep) , control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),
                                                        data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...))
              })
            }


            data_grouped$measure[data_grouped$group==group_id[1] & !is.na(data_grouped$measure)]<-residuals(model_ln_A)

            # single_rhythm B
            if(RRBS==TRUE)
            {
              suppressMessages({
                model_ln_B<-switch(lm_method,
                                   lm = lm(measure ~0+replicate_id + time:scaler ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...),
                                   lme = lmerTest::lmer(measure~0+ replicate_id + time:scaler +(1 | rep),
                                                        control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...))
              })
            }else{
              suppressMessages({
                model_ln_B<-switch(lm_method,
                                   lm = lm(measure ~0+ time  ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...),
                                   lme = lmerTest::lmer(measure~0+  time +(1 | rep),
                                                        control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...))
              })
            }

            data_grouped$measure[data_grouped$group==group_id[2] & !is.na(data_grouped$measure)]<-residuals(model_ln_B)

          }

          data_grouped$measure
        }
        outputs_fn
      }




    future:::ClusterRegistry("stop")

    results_out<-(do.call("cbind",lapply(rapply(res, enquote, how="unlist"), eval)))
    colnames(results_out)<-colnames(eset)
    rownames(results_out)<-rownames(eset)
    return(results_out)
  }else{
    if(verbose)cat("Performing single group inference","...\n")

    if(verbose)cat("Preparing design matrix ...\n")

    if(RRBS==TRUE)
    {

      design <- stats::model.matrix(~0+ replicate_id + time:scaler,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0+replicate_id+  time:scaler,
                                      data = exp_design)

    }else{

      design <- stats::model.matrix(~0+ time ,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0+  time,
                                      data = exp_design)
    }


    design_t<-design

    colnames(design) <- gsub(":", "_", colnames(design))

    if(is(data_input,"DGEList")| force_weight_estimation==TRUE){
      if(verbose)cat("Estimating weights for the input variables ...\n")

      if(lm_method=="lm")
      {
        if(verbose)cat("lm_method is lm, voom will be used ...\n")
        voom_res<-limma::voom(data_input, design, plot=FALSE)
      }else{

        if(RRBS==TRUE)
        {

          formula<-~0+ replicate_id  + time:scaler+(1 | rep)
        }else{
          formula<-~0+  time +(1 | rep)
        }


        voom_res<-variancePartition::voomWithDreamWeights(data_input,formula = formula,data = exp_design,BPPARAM = BiocParallel::SnowParam(workers = ncores))
      }
      obs_weights <-t(voom_res$weights)
      eset<-t(voom_res$E)
    }

    chunks <- parallel::splitIndices(ncol(eset), min(ncol(eset),  ncores))
    if(verbose)cat("Spliting data to",length(chunks)," chunks...\n")

    res<-foreach(chIndx=chunks)%dopar%
      {

        outputs_fn<-foreach(i=chIndx)%do%{



          data_grouped<-cbind(measure=as.numeric(eset[,i]),exp_design)


          if(RRBS==TRUE)
          {
            suppressMessages({
              model_ln<-switch(lm_method,
                               lm = lm(measure ~0+replicate_id + time:scaler ,data=data_grouped,weights = obs_weights[,i],...),
                               lme = lme4::lmer(measure~0+ replicate_id +time:scaler +(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
            })
          }else{
            suppressMessages({
              model_ln<-switch(lm_method,
                               lm = lm(measure ~0+time ,data=data_grouped,weights = obs_weights[,i],...),
                               lme = lme4::lmer(measure~0+ time +(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
            })
          }

          data_grouped$measure[!is.na(data_grouped$measure)]<-residuals(model_ln)
          data_grouped$measure
        }

        outputs_fn
      }

    future:::ClusterRegistry("stop")

  }
  results_out<-(do.call("cbind",lapply(rapply(res, enquote, how="unlist"), eval)))
  colnames(results_out)<-colnames(eset)
  rownames(results_out)<-rownames(eset)
  return(results_out)

}



