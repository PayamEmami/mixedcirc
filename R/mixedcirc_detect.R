#' Perform differential circadian rhythm analysis using mixed models
#'
#' This functions performs differential circadian rhythm analysis using mixed models.
#'
#' @param data_input A numerical matrix or data.frame (N*P) or DGEList where in the rows are samples (N) and the columns are variables (P). If DGEList is provided, regression weights will be estimated!
#' @param time A vector of length N, showing circadian time of each sample
#' @param group A character vector of length N. If performing differential circadian rhythm analysis, group is a factor, showing grouping of the samples. Analysis of two groups is supported at this stage! See details!
#' @param id A vector of length N showing identity of each *unique* sample. See details
#' @param period Period of circadian rhythm. Default: 24
#' @param lm_method The regression method to use. At this stage, `lm` and `lme` are supported! If lm is selected, normal regression will be performed. Default: "lme"
#' @param f_test Type of f-test for calculating p-value of the rhythm. Possible values are "multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger". Default: Satterthwaite
#' @param abs_phase Whether to return absolute phase or not. Default: TRUE
#' @param obs_weights Regression weights. Default: NULL. See details
#' @param RRBS If `TRUE`, the data is assumed to be RRBS methylation data. if TRUE, obs_weights must be set. Default FALSE
#' @param replicate_id If `RRBS` is set to `TRUE`, This has to be a factor showing identity of each unique replicate.
#' @param separate_tests If TRUE (defualt), separate models will be fit for each group to estimate rhythmicity otherwise, groupwise rhythmicity will be based on the a global model.
#' @param ncores number of cores
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @param ... additionl arguments to the regression function
#' @export
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#'
#' @return
#' A class of mixedcirc_fit_list.
#'
#' @details
#' For each variable we use the following mode:
#' In this part we do rhythmicity analysis on individual variables using the following model:
#' measure ~0 + group + group:in + group:out
#' Where `in` is defined as cos(2 * pi * time / period) and `out` as sin(2 * pi * time / period).
#' The global inference is tested on group:in == 0 or group:out == 0
#' The between group difference is tested as differences between in and out across different groups. No p-value adjustment is performed!
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

mixedcirc_detect <- function(data_input=NULL,time=NULL,group=NULL,id=NULL,
                             period=24,lm_method=c("lm","lme")[2],
                             f_test=c("multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger")[3],
                             abs_phase=TRUE,obs_weights=NULL,RRBS=FALSE,replicate_id=NULL,separate_tests=TRUE,ncores=1,verbose=FALSE,...){



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

  if(!is.numeric(period))
    stop("*period* must be numeric!")

  lm_method <- match.arg(lm_method,c("lm","lme"))
  f_test <- match.arg(f_test,c("multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger"))

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

  back_transform_method <- "card"

  if(f_test%in%c("Satterthwaite", "Kenward-Roger") & lm_method=="lm")
  {
    stop("Datasets without hierarchical structures do not work with Satterthwaite and Kenward-Roger! Change f_test to multcomp_f or multcomp_chi!")
  }


  if(multiple_groups==TRUE)
  {
    if(verbose)cat("Performing multigroup inference ...\n")

    group_id <- base::levels(exp_design$group)
    if(verbose)cat("Creating design matrix for ...\n")
    exp_design <- base::cbind(exp_design,
                              inphase = cos(2 * pi * exp_design$time / period),
                              outphase = sin(2 * pi * exp_design$time / period))


    if(RRBS==TRUE)
    {

      design <- stats::model.matrix(~0+ replicate_id+ group:scaler + group:inphase:scaler + group:outphase:scaler,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0 +replicate_id+  inphase:scaler + outphase:scaler,
                                      data = exp_design)

    }else{

      design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0 +  inphase + outphase,
                                      data = exp_design)
    }



    design_t<-design

    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design) <- gsub(":", "_", colnames(design))

    if(is(data_input,"DGEList")){
      if(verbose)cat("Estimating weights for the input variables ...\n")

      if(lm_method=="lm")
      {
        if(verbose)cat("lm_method is lm, voom will be used ...\n")
        voom_res<-limma::voom(data_input, design, plot=FALSE)
      }else{
        if(RRBS==TRUE)
        {

          formula<-~ 0+replicate_id + group:scaler + group:inphase:scaler + group:outphase:scaler+(1 | rep)
        }else{
          formula<-~ 0 + group + group:inphase + group:outphase+(1 | rep)
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

      if(verbose)cat("Fitting the model for variable",i,"...\n")
      if(RRBS==TRUE)
      {
        suppressMessages({
        model_ln<-switch(lm_method,
                         lm = lm(measure ~0 +replicate_id+ group:scaler + group:inphase:scaler + group:outphase:scaler ,data=data_grouped,weights = obs_weights[,i],...),
                         lme = lme4::lmer(measure~ 0 +replicate_id+ group:scaler + group:inphase:scaler + group:outphase:scaler+(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
})
      }else{
        suppressMessages({
        model_ln<-switch(lm_method,
                         lm = lm(measure ~0 + group + group:inphase + group:outphase ,data=data_grouped,weights = obs_weights[,i],...),
                         lme = lme4::lmer(measure~ 0 + group + group:inphase + group:outphase+(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
})
      }

      if(verbose)cat("Performing f-test for variable",i,"...\n")
      ## calculate f-test
      cof<-matrix(0,nrow = ncol(design_t),ncol = ncol(design_t))
      colnames(cof)<-rownames(cof)<-c(colnames(design_t))

      conts<-c()
      for(x in colnames(cof)[grep("phase", colnames(cof))])conts<-c(conts,paste(add_sym(x)," == 0",sep = ""))


      g <- multcomp::glht(model_ln, linfct = conts)

      f_test_results<-NULL
      if(f_test=="multcomp_chi")
      {
        f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
      }else if(f_test == "multcomp_f")
      {
        f_test_results<-multcomp:::summary.glht(g, test = Ftest())
      }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
      {
        f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
      }else{
        stop("Wrong f test method!")
      }


      f_p_value<-f_test_results$test$pvalue

      if(separate_tests)
      {
        # single_rhythm A
        if(RRBS==TRUE)
        {
          suppressMessages({
          model_ln_A<-switch(lm_method,
                             lm = lm(measure ~0 +replicate_id+ inphase:scaler + outphase:scaler ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...),
                             lme = lmerTest::lmer(measure~ 0 +replicate_id+ inphase:scaler + outphase:scaler+(1 | rep) ,
                                                  control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),
                                                  data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...))
})
        }else{
          suppressMessages({
          model_ln_A<-switch(lm_method,
                             lm = lm(measure ~0 + inphase + outphase ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...),
                             lme = lmerTest::lmer(measure~ 0 + inphase + outphase+(1 | rep) , control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),
                                                  data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights[data_grouped$group==group_id[1],i],...))
})
        }

        cof_s<-matrix(0,nrow = ncol(design_s),ncol = ncol(design_s))
        colnames(cof_s)<-rownames(cof_s)<-c(colnames(design_s))

        conts_g<-c()
        for(x in colnames(cof_s)[grep("phase", colnames(cof_s))])conts_g<-c(conts_g,paste(add_sym(x)," == 0",sep = ""))

        g <- multcomp::glht(model_ln_A, linfct = conts_g)

        f_test_results<-NULL
        if(f_test=="multcomp_chi")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
        }else if(f_test == "multcomp_f")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Ftest())
        }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
        {
          f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln_A,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
        }else{
          stop("Wrong f test method!")
        }

        f_p_value_A<-f_test_results$test$pvalue


        # single_rhythm B
        if(RRBS==TRUE)
        {
          suppressMessages({
          model_ln_B<-switch(lm_method,
                             lm = lm(measure ~0+replicate_id + inphase:scaler + outphase:scaler ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...),
                             lme = lmerTest::lmer(measure~ 0+replicate_id + inphase:scaler + outphase:scaler+(1 | rep),
                                                  control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...))
})
        }else{
          suppressMessages({
          model_ln_B<-switch(lm_method,
                             lm = lm(measure ~0 + inphase + outphase ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...),
                             lme = lmerTest::lmer(measure~ 0 + inphase + outphase+(1 | rep),
                                                  control = lmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e4)),data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights[data_grouped$group==group_id[2],i],...))
})
        }


        g <- multcomp::glht(model_ln_B, linfct = conts_g)

        f_test_results<-NULL
        if(f_test=="multcomp_chi")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
        }else if(f_test == "multcomp_f")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Ftest())
        }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
        {
          f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln_B,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
        }else{
          stop("Wrong f test method!")
        }

        f_p_value_B<-f_test_results$test$pvalue

      }else{

        # single_rhythm A

        model_ln_A<-model_ln
        cof_s<-matrix(0,nrow = ncol(design_s),ncol = ncol(design_s))
        colnames(cof_s)<-rownames(cof_s)<-c(colnames(design_s))
        conts_g<-c()
        for(x in colnames(cof_s)[grep("phase", colnames(cof_s))])conts_g<-c(conts_g,paste(add_sym(x)," == 0",sep = ""))

        conts_g<-conts[grepl(paste0(group_id[1],":"),x = conts,fixed = T)]

        g <- multcomp::glht(model_ln_A, linfct = conts_g)

        f_test_results<-NULL
        if(f_test=="multcomp_chi")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
        }else if(f_test == "multcomp_f")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Ftest())
        }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
        {
          f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln_A,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
        }else{
          stop("Wrong f test method!")
        }

        f_p_value_A<-f_test_results$test$pvalue


        # single_rhythm B


        model_ln_B<-model_ln
        cof_s<-matrix(0,nrow = ncol(design_s),ncol = ncol(design_s))
        colnames(cof_s)<-rownames(cof_s)<-c(colnames(design_s))

        conts_g<-c()
        for(x in colnames(cof_s)[grep("phase", colnames(cof_s))])conts_g<-c(conts_g,paste(add_sym(x)," == 0",sep = ""))

        conts_g<-conts[grepl(paste0(group_id[2],":"),x = conts,fixed = T)]



        g <- multcomp::glht(model_ln_B, linfct = conts_g)

        f_test_results<-NULL
        if(f_test=="multcomp_chi")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
        }else if(f_test == "multcomp_f")
        {
          f_test_results<-multcomp:::summary.glht(g, test = Ftest())
        }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
        {
          f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln_B,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
        }else{
          stop("Wrong f test method!")
        }

        f_p_value_B<-f_test_results$test$pvalue
      }


      ## prepare contrasts
      contrasts<-c(
        paste0(add_sym(colnames(design_t)[grepl(pattern = "inphase",x = colnames(design_t))&grepl(pattern = group_id[1],x = colnames(design_t))]),"-",
               add_sym(colnames(design_t)[grepl(pattern = "inphase",x = colnames(design_t))&grepl(pattern = group_id[2],x = colnames(design_t))]),
               "==0"),
        paste0(add_sym(colnames(design_t)[grepl(pattern = "outphase",x = colnames(design_t))&grepl(pattern = group_id[1],x = colnames(design_t))]),"-",
               add_sym(colnames(design_t)[grepl(pattern = "outphase",x = colnames(design_t))&grepl(pattern = group_id[2],x = colnames(design_t))]),
               "==0"))



      g_diff <- multcomp::glht(model_ln, linfct =contrasts)

      f_test_results_diff<-NULL
      if(f_test=="multcomp_chi")
      {
        f_test_results_diff<-multcomp:::summary.glht(g_diff, test = Chisqtest())
      }else if(f_test == "multcomp_f")
      {
        f_test_results_diff<-multcomp:::summary.glht(g_diff, test = Ftest())
      }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
      {
        f_test_results_diff<-list(test=list(pvalue=lmerTest::contest(model_ln,L=g_diff$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
      }else{
        stop("Wrong f test method!")
      }
      f_p_value_diff<-f_test_results_diff$test$pvalue

      ## back transform the coeffcients

      # extract coeffcients
      if(verbose)cat("Extracting and estimating coefficients",i,"...\n")
      ext_cof<-g_diff$linfct
      diag(ext_cof)<-1
      g_coef <- multcomp::glht(model_ln, linfct = ext_cof)
      coeffs<-g_coef$coef
      coeffs<-data.frame(t(coeffs),check.names = F)
      names(coeffs) <- gsub("group", "", names(coeffs))
      names(coeffs) <- gsub(":", "_", names(coeffs))




      if (any(base::grepl(paste0(group_id[1], "_"), names(coeffs)))) {
        rhy_params <- coeffs[, base::paste(group_id[1], c(ifelse(RRBS,yes = "scaler_inphase",no = "inphase"),
                                                          ifelse(RRBS,yes = "scaler_outphase",no = "outphase")), sep = "_")]
        mesor_A<-coeffs[,1]
        if(back_transform_method=="compareRhythms")
        {
          amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
          phases_A <- base::atan2(rhy_params[, 2], rhy_params[,
                                                              1])%%(2 * pi)
        }else if(back_transform_method=="cosinor"){
          amps_A <- sqrt(base::rowSums(rhy_params^2))
          phases_A<-base::atan(rhy_params[, 2]/ rhy_params[,
                                                           1])
        }else if(back_transform_method=="card"){
          mesor <- coeffs[1]
          amps_A <- sqrt(base::rowSums(rhy_params^2))

          sb <- sign(rhy_params[1])
          sg <- sign(rhy_params[2])

          theta <- atan(abs(rhy_params[, 2]/rhy_params[,1]))

          if ((sb == 1 | sb == 0) & sg == 1) {
            phi <- -theta
          }else if (sb == -1 & (sg == 1 | sg == 0)) {
            phi <- theta - pi
          }else if ((sb == -1 | sb == 0) & sg == -1) {
            phi <- -theta - pi
          }else if (sb == 1 & (sg == -1 | sg == 0)) {
            phi <- theta - (2 * pi)
          }

          phases_A<-phi
        }else if(back_transform_method=="edpclau"){
          mesor <- coeffs[1]
          amps_A <- sqrt(base::rowSums(rhy_params^2))

          cos_coeff <- (rhy_params[1])
          sin_coeff <- (rhy_params[2])

          acrophase <- atan(abs(rhy_params[, 2]/rhy_params[,1]))

          if (cos_coeff < 0 & sin_coeff >= 0) {
            acrophase <-   acrophase + pi
          }
          if (cos_coeff < 0 & sin_coeff < 0) {
            acrophase <- pi + acrophase
          }

          if (cos_coeff >= 0 & sin_coeff < 0) {
            acrophase <- 2*pi + acrophase
          }

          # if (cos_coeff >= 0 & sin_coeff >= 0) {
          #   acrophase <- acrophase + pi
          # }

          phases_A<-acrophase

        }else{
          amps_A <- 0
          phases_A <- 0
        }
      }else {
        amps_A <- 0
        phases_A <- 0
      }
      if (any(base::grepl(paste0(group_id[2], "_"), colnames(coeffs)))) {
        rhy_params <- coeffs[, base::paste(group_id[2], c(ifelse(RRBS,yes = "scaler_inphase",no = "inphase"),
                                                          ifelse(RRBS,yes = "scaler_outphase",no = "outphase")), sep = "_")]


        mesor_B<-coeffs[,2]
        if(back_transform_method=="compareRhythms")
        {
          amps_B <- 2 * sqrt(base::rowSums(rhy_params^2))
          phases_B <- base::atan2(rhy_params[, 2], rhy_params[,
                                                              1])%%(2 * pi)
        }else if(back_transform_method=="cosinor"){
          amps_B <- sqrt(base::rowSums(rhy_params^2))
          phases_B<-base::atan(rhy_params[, 2]/ rhy_params[,
                                                           1])
        }else if(back_transform_method=="card"){
          amps_B <- sqrt(base::rowSums(rhy_params^2))

          sb <- sign(rhy_params[1])
          sg <- sign(rhy_params[2])

          theta <- atan(abs(rhy_params[, 2]/rhy_params[,1]))

          if ((sb == 1 | sb == 0) & sg == 1) {
            phi <- -theta
          }else if (sb == -1 & (sg == 1 | sg == 0)) {
            phi <- theta - pi
          }else if ((sb == -1 | sb == 0) & sg == -1) {
            phi <- -theta - pi
          }else if (sb == 1 & (sg == -1 | sg == 0)) {
            phi <- theta - (2 * pi)
          }

          phases_B<-phi
        }else if(back_transform_method=="edpclau"){
          mesor <- coeffs[1]
          amps_B <- sqrt(base::rowSums(rhy_params^2))

          cos_coeff <- (rhy_params[1])
          sin_coeff <- (rhy_params[2])

          acrophase <- atan(abs(rhy_params[, 2]/rhy_params[,1]))

          if (cos_coeff < 0 & sin_coeff >= 0) {
            acrophase <-   acrophase + pi
          }
          if (cos_coeff < 0 & sin_coeff < 0) {
            acrophase <- pi + acrophase
          }

          if (cos_coeff >= 0 & sin_coeff < 0) {
            acrophase <- 2*pi + acrophase
          }
          phases_B<-acrophase

        }else{
          amps_B <- 0
          phases_B <- 0
        }
      }else {
        amps_B <- 0
        phases_B <- 0
      }

      if(verbose)cat("Preparing output for variable",i,"...\n")
      dt_out<-data.frame(mesor_A=mesor_A,mesor_B=mesor_B,amps_A=amps_A,amps_B=amps_B,phases_A=phases_A,phases_B=phases_B,global_p_value=f_p_value,diff_p_value=f_p_value_diff,f_p_value_A=f_p_value_A,f_p_value_B=f_p_value_B)
      rownames(dt_out)<-colnames(eset)[i]
      if(abs_phase){dt_out<-abs(dt_out)}
      colnames(dt_out)<-gsub(pattern = "_A",replacement = paste("_",group_id[1],sep = ""),x = colnames(dt_out),fixed = T)
      colnames(dt_out)<-gsub(pattern = "_B",replacement = paste("_",group_id[2],sep = ""),x = colnames(dt_out),fixed = T)
      if(verbose)cat("Variable",i,"finished\n")
      #list(results=dt_out,fit=model_ln,exp_design=exp_design)
     new("mixedcirc_fit",results = dt_out,fit=model_ln,exp_design=exp_design,type=ifelse(RRBS,yes="RRBS",no = "expression"))
      }
        outputs_fn
      }

    future:::ClusterRegistry("stop")

    res_tmp<-new("mixedcirc_fit_list",results = lapply(rapply(res, enquote, how="unlist"), eval))
    return((res_tmp))


  }else{
    if(verbose)cat("Performing single group inference","...\n")

    if(verbose)cat("Preparing design matrix ...\n")
    exp_design <- base::cbind(exp_design,
                              inphase = cos(2 * pi * exp_design$time / period),
                              outphase = sin(2 * pi * exp_design$time / period))


    if(RRBS==TRUE)
    {

      design <- stats::model.matrix(~0+ replicate_id + inphase:scaler + outphase:scaler,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0 +replicate_id+  inphase:scaler + outphase:scaler,
                                      data = exp_design)

    }else{

      design <- stats::model.matrix(~0 + inphase + outphase,
                                    data = exp_design)

      design_s <- stats::model.matrix(~0 +  inphase + outphase,
                                      data = exp_design)
    }


    design_t<-design

    colnames(design) <- gsub(":", "_", colnames(design))

    if(is(data_input,"DGEList")){
      if(verbose)cat("Estimating weights for the input variables ...\n")

      if(lm_method=="lm")
      {
        if(verbose)cat("lm_method is lm, voom will be used ...\n")
        voom_res<-limma::voom(data_input, design, plot=FALSE)
      }else{

        if(RRBS==TRUE)
        {

          formula<-~ 0+replicate_id  + inphase:scaler + outphase:scaler+(1 | rep)
        }else{
          formula<-~ 0 + inphase + outphase+(1 | rep)
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
      if(verbose)cat("Fitting the model for variable",i,"...\n")
      suppressMessages({
      model_ln<-switch(lm_method,
                       lm = lm(measure ~  inphase + outphase ,data=data_grouped,weights = obs_weights[,i],...),
                       lme = lme4::lmer(measure~ 0 + inphase + outphase+(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
})

      if(RRBS==TRUE)
      {
        suppressMessages({
        model_ln<-switch(lm_method,
                         lm = lm(measure ~0 +replicate_id+ scaler + inphase:scaler +outphase:scaler ,data=data_grouped,weights = obs_weights[,i],...),
                         lme = lme4::lmer(measure~ 0 +replicate_id +inphase:scaler + outphase:scaler+(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
})
      }else{
        suppressMessages({
        model_ln<-switch(lm_method,
                         lm = lm(measure ~0 + inphase + outphase ,data=data_grouped,weights = obs_weights[,i],...),
                         lme = lme4::lmer(measure~ 0 + inphase + outphase+(1 | rep) ,data=data_grouped,weights = obs_weights[,i],...))
})
      }

      if(verbose)cat("Calculating f-test for variable",i,"...\n")
      ## calculate f-test
      cof<-matrix(0,nrow = ncol(design_t),ncol = ncol(design_t))
      colnames(cof)<-rownames(cof)<-c(colnames(design_t))

      conts<-c()
      for(x in colnames(cof)[grep("phase", colnames(cof))])conts<-c(conts,paste(x," == 0",sep = ""))

      g <- multcomp::glht(model_ln, linfct = conts)

      f_test_results<-NULL
      if(f_test=="multcomp_chi")
      {
        f_test_results<-multcomp:::summary.glht(g, test = Chisqtest())
      }else if(f_test == "multcomp_f")
      {
        f_test_results<-multcomp:::summary.glht(g, test = Ftest())
      }else if(f_test%in%c("Satterthwaite", "Kenward-Roger"))
      {
        f_test_results<-list(test=list(pvalue=lmerTest::contest(model_ln,L=g$linfct,joint = TRUE,ddf=f_test)[,"Pr(>F)"]))
      }else{
        stop("Wrong f test method!")
      }

      f_p_value<-f_test_results$test$pvalue


      ## back transform the coeffcients

      # extract coeffcients
      if(verbose)cat("Extracting and estimating coefficients for variable ",i,"...\n")
      ext_cof<-cof
      diag(ext_cof)<-1
      g_coef <- multcomp::glht(model_ln)
      coeffs<-g_coef$coef
      coeffs<-data.frame(t(coeffs),check.names = F)
      names(coeffs) <- gsub(":", "_", names(coeffs))


      if (all(base::is.element(c("inphase", "outphase"), colnames(coeffs)))) {


        rhy_params <- coeffs[, c(ifelse(RRBS,yes = "scaler_inphase",no = "inphase"), ifelse(RRBS,yes = "scaler_outphase",no = "outphase"))]
        amps <- 2 * sqrt(base::rowSums(rhy_params^2))
        phases <- base::atan2(rhy_params[, 2], rhy_params[, 1])%%(2 *
                                                                    pi)

        if(back_transform_method=="compareRhythms")
        {
          amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
          phases_A <- base::atan2(rhy_params[, 2], rhy_params[,
                                                              1])%%(2 * pi)
        }else if(back_transform_method=="cosinor"){
          amps_A <- sqrt(base::rowSums(rhy_params^2))
          phases_A<-base::atan(rhy_params[, 2]/ rhy_params[,
                                                           1])%%(2 * pi)
        }else if(back_transform_method=="card"){
          mesor <- coeffs[1]
          amps_A <- sqrt(base::rowSums(rhy_params^2))

          sb <- sign(rhy_params[1])
          sg <- sign(rhy_params[2])

          theta <- atan(abs(rhy_params[, 2]/rhy_params[,1]))

          if ((sb == 1 | sb == 0) & sg == 1) {
            phi <- -theta
          }else if (sb == -1 & (sg == 1 | sg == 0)) {
            phi <- theta - pi
          }else if ((sb == -1 | sb == 0) & sg == -1) {
            phi <- -theta - pi
          }else if (sb == 1 & (sg == -1 | sg == 0)) {
            phi <- theta - (2 * pi)
          }

          phases_A<-phi
        }else{
          amps_A <- 0
          phases_A <- 0
        }
      }else {
        amps_A <- 0
        phases_A <- 0
      }

      if(verbose)cat("Preparing output for variable ",i,"...\n")
      dt_out<-data.frame(amps_A=amps_A,phases_A=phases_A,pglobal_p_value=f_p_value)
      rownames(dt_out)<-colnames(eset)[i]
      if(abs_phase){dt_out<-abs(dt_out)}

      if(verbose)cat("Variable ",i," finished...\n")
      #list(results=dt_out,fit=model_ln,exp_design=exp_design)
      new("mixedcirc_fit",results = dt_out,fit=model_ln,exp_design=exp_design,type=ifelse(RRBS,yes="RRBS",no = "expression"))

        }

        outputs_fn
      }

    future:::ClusterRegistry("stop")

    res_tmp<-new("mixedcirc_fit_list",results = lapply(rapply(res, enquote, how="unlist"), eval))
    return((res_tmp))

  }

}



