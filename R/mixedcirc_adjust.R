#' Perform p-value adjustment on mixedcirc_fit_list
#'
#' This functions performs p-value adjustment on mixedcirc_fit_list
#'
#' @param input An object of mixedcirc_fit_list. It has to contain more than one variable
#' @param method A vector of character strings containing the names of the multiple testing procedures for which adjusted p-values are to be computed. This vector should include any of the following: "holm","hochberg","hommel","bonferroni", "BH", "BY","fdr" "none"
#' @param pattern Columns of results with this pattern will be selected for doing adjustment. Default: p_value
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @param ... additional arguments passed to mt.rawp2adjp from multtest package
#' @export
#' @examples
#' data("circa_data")
#'
#'results<-mixedcirc_detect(data_input = circa_data$data_matrix,
#'time = circa_data$time,group = circa_data$group,id = circa_data$id,period = 24,verbose = TRUE)
#'results<-mixedcirc_adjust(results,method="Bonferroni")
#'
#' @return
#' A class of mixedcirc_fit_list in which the adjusted p-values have been added in the results
#'
#' @details
#' The p-values from the matching columns will be adjusted one column at the time using the selected method.
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
#' @import multtest

mixedcirc_adjust <- function(input, method="BH", pattern="p_value", verbose=FALSE, ...) {

  if (verbose) {
    cat("Checking inputs ...\n")
  }
  # checking inputs
  if (is(input,"mixedcirc_fit_list.")) {
    stop("input must be a data frame or matrix")
  }

  if(!method%in%c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                  "fdr", "none")) {
    stop("method must be one of Bonferroni, holm, hochberg, hommel, bonferroni, BH, BY,
                  fdr, none")
  }

  if(!is.character(pattern)) {
    stop("pattern must be a string (character)")
  }


  if(verbose) {
    cat("gathering the data ...\n")
  }

  data_gathered <- input@results
  data_gathered <- do.call(rbind,lapply(data_gathered,function(x){ x@results }))

  selected_pvalues <- grepl(pattern = pattern,x = colnames(data_gathered))

  if(sum(selected_pvalues)==0) {
    stop("Did not find any column matching the provided pattern")
  }

  selected_data<-data_gathered[ ,selected_pvalues,drop=F]

  if (verbose) {
    cat("Performing adjustment ...\n")
  }

  adjusted_pvalues<-apply(selected_data,MARGIN = 2,function(x){
    matrix(p.adjust(p = x,method = method))
  })[,,drop=F]

  colnames(adjusted_pvalues)<-paste(colnames(selected_data),"_adjusted_",method,sep="")

  if (verbose) {
    cat("Preparing the results ...\n")
  }

  for (i in 1:length(input)) {
    input@results[[i]]@results<-cbind(input@results[[i]]@results[,,drop=F],adjusted_pvalues[i,,drop=F])
  }

  if(verbose) {
    cat("Finished\n")
  }

return((input))
}

