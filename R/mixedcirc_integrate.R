#' Circadian rhythm integration
#'
#' This functions performs circadian rhythm data integration using canonical correlations.
#'
#' @param data_input A named list of numerical matrices or data.frames (N*P) where in the rows are samples (N) and the columns are variables (P). Each element of the list is a data set
#' @param group A character vector of length N. If performing differential circadian rhythm analysis, group is a factor, showing grouping of the samples. Analysis of two groups is supported at this stage! See details!
#' @param id A vector of length N showing identity of each *unique* sample. See details
#' @param period Period of circadian rhythm. Default: 24
#' @param ncomp Number of components. If set to 0, maximum number of possible components will be calculated (default: 2)
#' @param variables A named list where in each element (one per data set) is a vector of length ncomp. Each element of the vector must be the number of variables retained on that particular component
#' @param lm_method The regression method to use. At this stage, `lm`,`lme`,`nlme` are supported! If lm is selected, normal regression will be performed. Default: "lme"
#' @param f_test Type of f-test for calculating p-value of the rhythm. Possible values are "multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger". Default: Satterthwaite
#' @param abs_phase Whether to return absolute phase or not. Default: TRUE
#' @param obs_weights Regression weights. Default: NULL. See details
#' @param decompose decomposes the Within variation in the data set with respect to id (default: FALSE)
#' @param center Centers the columns of each data set to zero
#' @param scale scale each data set by dividing each column by its standard deviation
#' @param merge If TRUE, the groups and the type will be merge in a single matrix for regression (default: FALSE)
#' @param no_correlation If TRUE, zero covariance will be set between time and groups (default: FALSE)
#' @param no_interaction If TRUE, no interaction between groups and time is assumed. (default FALSE)
#' @param max.iter Maximum number of iteration of model (default: 100000)
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @param ... additionl arguments to the regression function
#' @export
#' @examples
#'library(mixedcirc)
#'data("circa_data")
#'data_input<-list(a=circa_data$data_matrix[,1:3],b=circa_data$data_matrix[,3:7])
#'results<-mixedcirc_integrate(data_input,time = circa_data$time,group = circa_data$group,id = circa_data$id)
#'data_matrix<-mixedcirc_getscore(ee,type = "partial",merge = T)
#'fitted_data<-mixedcirc_detect(data_matrix,time = circa_data$time,group = circa_data$group,id = circa_data$id)
#'plot(fitted_data[1])
#'
#'
#' @return
#' A class of \code{\link{mixedcirc_integration}}.
#'
#' @details
#' This method is based on SGCCA analysis
#' All the data sets will be correlated together with the linearized rhythm and the group
#' The resulting partial or average scores can be used in \code{\link{mixedcirc_detect}} to perform differential circadian rhythm analysis using mixed models
#' The loadings can be used for example to cluster the corresponding variables across different domains
#' The scores can be retrieved using \code{\link{mixedcirc_getscore}}
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



mixedcirc_integrate <- function(data_input=NULL,time=NULL,group=NULL,id=NULL,
                             period=24,ncomp=2,variables=NULL,lm_method=c("lm","lme","nlme")[2],
                             f_test=c("multcomp_f","multcomp_chi","Satterthwaite", "Kenward-Roger")[3],
                             abs_phase=TRUE,decompose=FALSE,center=TRUE,scale=FALSE,merge=FALSE,no_correlation=FALSE,
                             no_interaction=FALSE,
                             max.iter = 100000,verbose=FALSE,...){

  if(verbose)cat("Checking inputs ...\n")
  # checking inputs
  if(is.null(data_input))
    stop("data_input must be a list of one or several data frames or matrices")

  if(!is.list(data_input))
    stop("data_input must be a list of one or several data frames or matrices")

  if(is.null(names(data_input)))
    stop("data_input must be a named list (each element must have a name)")

  check_each_element<-sapply(data_input,function(x){is.matrix(x)|is.data.frame(x)})

  if(any(!check_each_element))
    stop("All elements of data_input data frame or matrix."," Check (", paste(names(data_input)[!check_each_element],sep = " and "),")")

  if(is.null(ncomp) | is.na(ncomp) | is.infinite(ncomp))
    stop("ncomp must be a positive integer!")

  if(ncomp<=0)
    stop("ncomp must be a positive integer!")

  if(ncomp%%1!=0)
    stop("ncomp must be a positive integer!")

  if(!is.null(variables))
  {
    if(is.null(names(variables)))
      stop("data_input must be a named list (each element must have a name)")

    if(!all(names(variables)%in%names(data_input)))
      stop("variables must be a named list (wih the same names as data_input)")

    check_each_element<-sapply(variables,function(x){is.vector(x)})

    if(any(!check_each_element))
      stop("All elements of variables must be data frame or matrix."," Check (", paste(names(variables)[!check_each_element],sep = " and "),")")

    check_each_element<-sapply(variables,function(x){length(x)<ncomp})
    if(any(check_each_element))
      warning("Each element of variables should have ",ncomp," elements"," selecting all variables for ", paste(names(variables)[check_each_element],sep = " and "))
  }


  nrows_unique<-unique(sapply(data_input,nrow))

  if(length(nrows_unique)>1)
    stop("data_input must contain data frames or matrices with the exact same number of rows!")

  if(is.null(time))
    stop("time must be a vector")

 if(length(time)!=nrows_unique)
      stop("The length of *time* is not equal to the number of rows of *data_input*")



  if(!is.null(group))
  {

    if(length(group)!=nrows_unique)
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

   if(length(id)!=nrows_unique)
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

  if(verbose)cat("Building experiment ...\n")
  exp_design<-cbind.data.frame(time=as.numeric(time),group=as.factor(group),rep=as.character(id))

  if(!is.factor(exp_design$group))stop("Group must be a factor!")


  if(!is.null(id) & decompose){
    if(verbose)
    cat("Performing the Within variation decomposition ...\n")
    data_input<-lapply(data_input,mixOmics::withinVariation)
  }
  if(verbose)
  cat("Performing scaling ...\n")
  data_input<-lapply(data_input,function(x){scale(x,center=center,scale=scale)})
  if(verbose)
  cat("Building experimental design ...\n")


  back_transform_method <- "card"
  if(multiple_groups==TRUE)
  {


    group_id <- base::levels(exp_design$group)
    exp_design <- base::cbind(exp_design,
                              inphase = cos(2 * pi * exp_design$time / period),
                              outphase = sin(2 * pi * exp_design$time / period))

      design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                    data = exp_design)



      if(merge)
      {
        Y_time<-design
        Y_group<-NULL
        if(no_interaction)
        {
          Y_time<-exp_design
        }


      }else{
        Y_time<-design[,-c(1,2)]
        Y_group<-design[,c(1,2)]
        if(no_interaction)
        {
          Y_time<-exp_design
        }
      }


  }else{
    exp_design <- base::cbind(exp_design,
                              inphase = cos(2 * pi * exp_design$time / period),
                              outphase = sin(2 * pi * exp_design$time / period))

      design <- stats::model.matrix(~0 + inphase + outphase,
                                    data = exp_design)
      Y_group <- NULL
      Y_time<-design
    }

  if(verbose)
  cat("Building design matrix ...\n")
  data_input$Y_time <- Y_time
  if(!is.null(Y_group))
    data_input$Y_group <- Y_group


  design_matrix<-diag(length(data_input))
  rownames(design_matrix)<-names(data_input)
  colnames(design_matrix)<-names(data_input)

  design_matrix[,"Y_time"]<-1
  design_matrix["Y_time",]<-1
  if(!is.null(Y_group))
  {
    design_matrix[,"Y_group"]<-1
    design_matrix["Y_group",]<-1
  }

  if(no_correlation)
  {
    design_matrix["Y_time","Y_group"]<-0
    design_matrix["Y_group","Y_time"]<-0
  }

  if(is.null(variables))
  {
    variables<-lapply(data_input,function(x){rep(ncol(x),ncomp)})
  }else{

    data_input$Y_time <- Y_time

    variables$Y_time<-rep(ncol(Y_time),ncomp)
    if(!is.null(Y_group))
      variables$Y_group<-rep(ncol(Y_group),ncomp)
  }

  if(verbose)
  cat("Running the model ...\n")
  sgcca_model<-mixOmics::wrapper.sgcca(data_input,ncomp = ncomp,scale = F,design = design_matrix,mode = "canonical",
                                      keepX  =variables,max.iter = 100000)
  if(verbose)
  cat("Extracting data ...\n")
  partial_scores<-sgcca_model$variates[-which(names(sgcca_model$variates)%in%c("Y_time","Y_group"))]
  partial_loadings<-sgcca_model$variates[-which(names(sgcca_model$loadings)%in%c("Y_time","Y_group"))]
  if(verbose)
  cat("Calculating the average score ...\n")
  block_object<-sgcca_model
  X_blocks <- with(block_object, names$blocks[-which(names$block %in%
                                                       c("Y_group","Y_time"))])

  .get_average_variates <- function(object, X_blocks, weighted = FALSE) {
    arrays <- object$variates
    arrays <- arrays[X_blocks]
    if (weighted) {

      weights <- object$weights
    }
    else {
      weights <- matrix(rep(1, length(X_blocks) * object$ncomp[1]),
                        nrow = length(X_blocks))
      dimnames(weights) <- list(X_blocks, paste0("comp",
                                                 seq_len(object$ncomp[1])))
      weights <- as.data.frame(weights)
    }
    block_names <- mixOmics:::.name_list(names(arrays))
    weighted_arrays <- lapply(block_names, function(x) {
      variates <- arrays[[x]]
      weights <- diag(weights[x, ])
      weighted_variates <- variates %*% weights
      dimnames(weighted_variates) <- dimnames(variates)
      weighted_variates
    })
    wtd_sum <- Reduce(f = "+", weighted_arrays)
    sweep(wtd_sum, MARGIN = 2, colSums(weights), FUN = "/")
  }
  average_blocks<-"average"
  for (average_block in average_blocks) {
      block_object$variates[[average_block]] <- .get_average_variates(object = block_object,
                                                                      X_blocks = X_blocks, weighted = FALSE)
      block_object$variates$avg<- block_object$variates[[average_block]]
    block_object$names$blocks <- c(block_object$names$blocks,
                                   average_block)
    block_object$ncomp[average_block] <- block_object$ncomp[1]
    block_object$prop_expl_var[average_block] <- 0
  }
  average_score<-block_object$variates$avg



  output<-new("mixedcirc_integration",partial = partial_scores, average=average_score, loadings=partial_loadings,model=sgcca_model)
  if(verbose)
  cat("Finished!\n")

return(output)
}

#' Retrieved scores from mixedcirc_integration
#'
#' This functions performs circadian rhythm data integration using canonical correlations.
#'
#' @param input A mixedcirc_integration object
#' @param type what scores to retrieve (partial, average, loadings)
#' @param dataset name of the data set to extract the data from
#' @param merge If TRUE, the scores or loadings will be merged into a single data frame (default: FALSE). This has no effect if \code{type} is average!
#' @param verbose Enables printing logs
#' @export
#' @examples
#'library(mixedcirc)
#'data("circa_data")
#'data_input <- list(a=circa_data$data_matrix[,1:3],b=circa_data$data_matrix[,3:7])
#'results <- mixedcirc_integrate(data_input,time = circa_data$time,group = circa_data$group,id = circa_data$id)
#'score_matrix <- mixedcirc_getscore(results,type="average")
#' @return
#' A data frame of values
#'
#' @details
#' This method is used to extract scores and loadings from a mixedcirc_integration object
#'
mixedcirc_getscore<-function(input=NULL,type=NULL,dataset=NULL,merge=FALSE,verbose=FALSE)
{
  if(verbose)
    cat("Checking parameters ...!\n")
if(!is(input,"mixedcirc_integration"))
  stop("input must be mixedcirc_integrate!")

  if(is.null(type))stop("type need to be provided!")
  if(!type%in%c("partial","average","loadings"))
    stop("type must be one of partial, average or loadings!")

  if(!merge){
    if( is.null(dataset))stop("dataset needs to be provided or set the merge to TRUE")
  }else{
    dataset <- names(input@partial)[1]
  }

  if(!dataset%in%c(names(input@partial),"average"))
    stop("dataset must be one of ", paste(c(names(input@partial),"average"),collapse = ", "))

  if(verbose)
    cat("Extracting data ...!\n")
  if(dataset=="average")
    return(as.matrix(input@average))

  if(merge){
    return(as.matrix(dplyr::bind_cols(slot(input,type))))
  }else{
    return(as.matrix(slot(input,type)[[dataset]]))
  }
}

