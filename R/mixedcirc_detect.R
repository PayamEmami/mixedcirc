#' Perform differential circadian rhythm analysis using mixed models
#'
#' This functions performs differential circadian rhythm analysis using mixed models.
#'
#' @param data_input A numerical matrix or data.frame (N*P) where in the rows are samples (N) and the columns are variables (P)
#' @param time A vector of length N, showing circadian time of each sample
#' @param group A character vector of length N. If performing differential circadian rhythm analysis, group is a factor, showing grouping of the samples. Analysis of two groups is supported at this stage! See details!
#' @param id A vector of length N showing identity of each *unique* sample. See details
#' @param period Period of circadian rhythm. Default: 24
#' @param lm_method The regression method to use. At this stage, `lm`,`lme`,`nlme` are supported! If lm is selected, normal regression will be performed. Default: "lme"
#' @param f_test Type of f-test for calculating p-value of the rhythm. Possible values are "multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger". Default: Satterthwaite
#' @param abs_phase Whether to return absolute phase or not. Default: TRUE
#' @param obs_weights Regression weights. Default: NULL. See details
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
#' A nested list. Each element of the first list contains a data frame (result) and a fit class. The data frame includes the estimated mesor, amplitude and phase for each of the groups.
#' The data frame also includes the p-values for each estimates as well as differences between the groups.
#' The fit object is the actual model fitted.
#'
#' @details
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
#'



mixedcirc_detect <- function(data_input=NULL,time=NULL,group=NULL,id=NULL,
                                    period=24,lm_method=c("lm","lme","nlme")[2],
                                    f_test=c("multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger")[3],
                                    abs_phase=TRUE,obs_weights=NULL,verbose=FALSE,...){
  doFuture::registerDoFuture()
  future::plan(future::multisession)

  if(verbose)cat("Checking inputs ...\n")
  # checking inputs
  if(is.null(data_input))
    stop("data_input must be a data frame or matrix")
  if(!is.matrix(data_input) | !is.matrix(data_input))
    stop("data_input must be a data frame or matrix")

  if(is.null(time))
    stop("time must be a vector")

  if(length(time)!=nrow(data_input))
    stop("The length of *time* is not equal to the number of rows of *data_input*")

  if(!is.null(group))
  {
    if(length(group)!=nrow(data_input))
      stop("The length of *group* is not equal to the number of rows of *data_input*")

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
    if(length(id)!=nrow(data_input))
      stop("The length of *id* is not equal to the number of rows of *data_input*")
  }else{
    warning("No id information was provided! Switching to normal linear regression!")
    lm_method <- "lm"
    id<-"dummy"
  }

  if(!is.numeric(period))
    stop("*period* must be numeric!")

  lm_method <- match.arg(lm_method,c("lm","lme","nlme"))
  f_test <- match.arg(f_test,c("multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger"))

  eset <- data_input[,,drop=F]

  if(verbose)cat("Building experiment ...\n")
  exp_design<-cbind.data.frame(time=as.numeric(time),group=as.factor(group),rep=as.character(id))

  if(!is.factor(exp_design$group))stop("Group must be a factor!")

  back_transform_method <- "card"
  obs_weights<-NULL
  if(multiple_groups==TRUE)
  {
    if(verbose)cat("Performing multigroup inference ...\n")
    res<-foreach(i=1:ncol(eset))%do%{#ncol(eset)

      if(verbose)cat("Processing variable",i,"...\n")
      group_id <- base::levels(exp_design$group)
      if(verbose)cat("Creating design matrix for variable",i,"...\n")
      exp_design <- base::cbind(exp_design,
                                inphase = cos(2 * pi * exp_design$time / period),
                                outphase = sin(2 * pi * exp_design$time / period))

      if ("batch" %in% colnames(exp_design)) {

        design <- stats::model.matrix(~0 + group + group:inphase + group:outphase + batch,
                                      data = exp_design)

        design_s <- stats::model.matrix(~0 +  inphase + outphase,
                                        data = exp_design)

      } else {

        design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                      data = exp_design)
        design_s <- stats::model.matrix(~0 +  inphase + outphase,
                                        data = exp_design)
      }

      design_t<-design

      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design) <- gsub(":", "_", colnames(design))

      data_grouped<-cbind(measure=as.numeric(eset[,i]),exp_design)

      if(verbose)cat("Fitting the model for variable",i,"...\n")
      model_ln<-switch(lm_method,
                       lm = lm(measure ~0 + group + group:inphase + group:outphase ,data=data_grouped,weights = obs_weights,...),
                       lme = lme4::lmer(measure~ 0 + group + group:inphase + group:outphase+(1 | rep) ,data=data_grouped,weights = obs_weights,...),
                       nlme = nlme::lme(measure~ 0 + group + group:inphase + group:outphase,random=~1 | rep,data=data_grouped,weights = obs_weights,...))

      if(verbose)cat("Performing f-test for variable",i,"...\n")
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

      # single_rhythm A

      model_ln_A<-switch(lm_method,
                         lm = lm(measure ~0 + inphase + outphase ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights,...),
                         lme = lme4::lmer(measure~ 0 + inphase + outphase+(1 | rep) ,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights,...),
                         nlme = nlme::lme(measure~ 0 + inphase + outphase,random=~1 | rep,data=data_grouped[data_grouped$group==group_id[1],],weights = obs_weights,...))

      cof_s<-matrix(0,nrow = ncol(design_s),ncol = ncol(design_s))
      colnames(cof_s)<-rownames(cof_s)<-c(colnames(design_s))

      conts_g<-c()
      for(x in colnames(cof_s)[grep("phase", colnames(cof_s))])conts_g<-c(conts_g,paste(x," == 0",sep = ""))

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
      model_ln_B<-switch(lm_method,
                         lm = lm(measure ~0 + inphase + outphase ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights,...),
                         lme = lme4::lmer(measure~ 0 + inphase + outphase+(1 | rep) ,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights,...),
                         nlme = nlme::lme(measure~ 0 + inphase + outphase,random=~1 | rep,data=data_grouped[data_grouped$group==group_id[2],],weights = obs_weights,...))


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


      ## prepare for contrast

      contrasts <- c(paste0(group_id, "_inphase", collapse = "-"),
                     paste0(group_id, "_outphase", collapse = "-"))

      diff_rhy_contrast <- limma::makeContrasts(contrasts = contrasts,
                                                levels = design)

      g_diff <- multcomp::glht(model_ln, linfct = t(diff_rhy_contrast))

      f_test_results_diff<-NULL
      if(f_test=="chi")
      {
        f_test_results_diff<-multcomp:::summary.glht(g_diff, test = Chisqtest())
      }else if(f_test == "f")
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
      ext_cof<-cof
      diag(ext_cof)<-1
      g_coef <- multcomp::glht(model_ln, linfct = ext_cof)
      coeffs<-g_coef$coef
      coeffs<-data.frame(t(coeffs),check.names = F)
      names(coeffs) <- gsub("group", "", names(coeffs))
      names(coeffs) <- gsub(":", "_", names(coeffs))


      if (any(base::grepl(paste0(group_id[1], "_"), names(coeffs)))) {
        rhy_params <- coeffs[, base::paste(group_id[1], c("inphase",
                                                          "outphase"), sep = "_")]
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
        rhy_params <- coeffs[, base::paste(group_id[2], c("inphase",
                                                          "outphase"), sep = "_")]

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
      if(abs_phase){dt_out<-abs(dt_out)}
      colnames(dt_out)<-gsub(pattern = "_A",replacement = paste("_",group_id[1],sep = ""),x = colnames(dt_out),fixed = T)
      colnames(dt_out)<-gsub(pattern = "_B",replacement = paste("_",group_id[2],sep = ""),x = colnames(dt_out),fixed = T)
      if(verbose)cat("Variable",i,"finished\n")
      #list(results=dt_out,fit=model_ln,exp_design=exp_design)
      new("mixedcirc_fit",results = dt_out,fit=model_ln,exp_design=exp_design)

    }
    future:::ClusterRegistry("stop")
    return((res))


  }else{
    if(verbose)cat("Performing single group inference","...\n")
    res<-foreach::foreach(i=1:ncol(eset))%do%{#ncol(eset)

      if(verbose)cat("Preparing design matrix for variable",i,"...\n")
      exp_design <- base::cbind(exp_design,
                                inphase = cos(2 * pi * exp_design$time / period),
                                outphase = sin(2 * pi * exp_design$time / period))

      if ("batch" %in% colnames(exp_design)) {

        design <- stats::model.matrix(~0 + inphase + outphase,
                                      data = exp_design)
      } else {

        design <- stats::model.matrix(~0 + inphase + outphase,
                                      data = exp_design)
      }

      design_t<-design

      colnames(design) <- gsub(":", "_", colnames(design))

      data_grouped<-cbind(measure=as.numeric(eset[,i]),exp_design)
      if(verbose)cat("Fitting the model for variable",i,"...\n")
      model_ln<-switch(lm_method,
                       lm = lm(measure ~  inphase + outphase ,data=data_grouped,weights = obs_weights,...),
                       lme = lme4::lmer(measure~ 0 + inphase + outphase+(1 | rep) ,data=data_grouped,weights = obs_weights,...),
                       nlme = nlme::lme(measure~ 0 + inphase + outphase,random=~1 | rep,data=data_grouped,weights = obs_weights,...))

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


        rhy_params <- coeffs[, c("inphase", "outphase")]
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
      if(abs_phase){dt_out<-abs(dt_out)}

      if(verbose)cat("Variable ",i," finished...\n")
      #list(results=dt_out,fit=model_ln,exp_design=exp_design)
      new("mixedcirc_fit",results = dt_out,fit=model_ln,exp_design=exp_design)

    }
    future:::ClusterRegistry("stop")
    return(res)

  }

}
