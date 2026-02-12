# most of the functions here have been adapted from VarianceParition package


#' Transform RRBS data for linear mixed modelling (experimental)
#'
#' This function is used to derive weights for RRBS data using linear mixed modelling.
#'
#' @param counts a class of methylBaseDB
#' @param formula specifies variables for the linear (mixed) model. Must only specify covariates, e.g.: ~ a + b + (1|c). Formulas with only fixed effects also work, and fit() is run.
#' @param data data.frame with columns corresponding to formula
#' @param lib.size numeric vector containing total library sizes for each sample. This will be estimated based on RRBS data
#' @param chunk.size the number of chunk to retrieve from the data at the same time
#' @param span width of the lowess smoothing window as a proportion.
#' @param plot logical, should a plot of the mean-variance trend be displayed?
#' @param save.plot logical, should the coordinates and line of the plot be saved in the output?
#' @param quiet suppress message, default FALSE
#' @param BPPARAM parameters for parallel evaluation
#' @param ignore_missing if TRUE, the missing values will be left as they are. Otherwise, they will be imputed by zero. Default: False
#' @param ... other arguments are passed to lmer
#'
#' @return
#' A two element list each of which is methylBaseDB (same dimentions).
#' * E: log2 transformed expression
#' * W: weights
#'
#' @details
#' In order to run this function correctly, one has to create an experimental design desired for rhythm analysis.
#' The way that the regression is formulated for RRBS data is
#' ~0+ replicate_id+ group:scaler + group:in:scaler + group:out:scaler
#' this means that `data` must contain replicate_id, scaler, in and out.
#' `replicate_id` is a factor showing identity of each unique replicate.
#' `scaler` is a numeric which is one (1) for methylated and zero (0) for unmethylated data
#' `in` is cos(2 * pi * time / period)
#' `out` is sin(2 * pi * time / period)
#' `group` is an optional (if differential rhythm analysis is done) two level factor that shows the grouping of the samples.
#' If constructed correctly, the output of the function will be a list of two elemens (methylBaseDB) which can be used to extract the transformed methylated/unmethylated as well as the weights.
#' The reason for this function is to solve the problem with high memory demand of RRBS data. If you have just a few loci, use can use `voomWithDreamWeights` or `voom` essentially using the same setup.
#' By default NAs are imputed by zero. However if ignore_missing, the missing values are ignores. The later is not according to voom recommendation. Use at your own risk!
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



mixedcirc_rrbs_voom_mixed<-function (counts, formula, data, lib.size = NULL,chunk.size=100,
                               span = 0.5, plot = FALSE, save.plot = FALSE, quiet = FALSE,
                               BPPARAM = bpparam(),ignore_missing=FALSE, ...)
{


  formula = as.formula(formula)
  out <- list()
  design = NULL

  n <- counts@num.records
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(lib.size)){
    tbxFile_t<-methylKit::getDBPath(counts)
    lib.size_partial<-methylKit:::applyTbxByChunk(tbxFile = tbxFile_t,
                                                  chunk.size = chunk.size,return.type = "data.frame",FUN = function(x){
                                                    if(!ignore_missing)
                                                    {
                                                      x[is.na(x)]<-0;
                                                    }



                                                    as.data.frame(t(colSums(x[,counts@coverage.index],na.rm = TRUE)))[,,drop=F]


                                                  })


    lib.size<-colSums(lib.size_partial,na.rm = TRUE)

    # double the library one for meth and unmeth
    lib.size <- rep(lib.size,each=2)
  }

  ## do y transformation
  dir_for_file<-dirname(methylKit::getDBPath(counts))
  file_name <- paste(tools::file_path_sans_ext(basename(methylKit::getDBPath(counts))),"_mixedcirc","_libnorm","",sep = "")


  selec_sm<-
    c(rbind(counts@numCs.index,counts@numTs.index))

  select_all_info<-
    c(rbind(counts@coverage.index ,
            counts@numCs.index,counts@numTs.index))

  mBase<-counts
  slotList <- list(dbtype = mBase@dbtype,
                   sample.ids = mBase@sample.ids,assembly = mBase@assembly,
                   context = mBase@context,resolution = mBase@resolution,
                   treatment = mBase@treatment,destranded = mBase@destranded,
                   coverage.index = mBase@coverage.index,
                   numCs.index = mBase@numCs.index,
                   numTs.index = mBase@numTs.index)

  tabixHead <- methylKit:::makeTabixHeader(slotList)
  tabixHeadString <- methylKit:::.formatTabixHeader(class = "methylBaseDB",
                                                    tabixHead = tabixHead)


  y<-methylKit:::applyTbxByChunk(tbxFile = tbxFile_t,tabixHead = tabixHeadString,
                                 chunk.size = chunk.size,return.type = "tabix",dir = dir_for_file,filename = file_name,
                                 FUN = function(x){

                                   xx<-x[,-c(1:4,counts@coverage.index),drop=F];
                                   if(!ignore_missing)
                                   {
                                   xx[is.na(xx)]<-0;
                                   }
                                   xx_res<-x;
                                   x2<-x;
                                   xx<-t(log2(t(xx + 0.5)/(lib.size + 1) * 1e+06));
                                   xx_res[,-c(1:4,counts@coverage.index)]<-
                                     xx;
                                   xx_res
                                 })

  y_norm <- methylKit:::readMethylBaseDB(dbpath = y,
                                         sample.ids = slotList$sample.ids,
                                         assembly = counts@assembly,
                                         dbtype = "tabix",
                                         resolution = counts@resolution,
                                         context = counts@context,
                                         treatment = slotList$treatment,destranded = F)


  if (.isMixedModelFormula(formula)) {
    if (missing(data)) {
      stop("Must specify argument 'data'\n")
    }






    selec_sm<-
      c(rbind(counts@numCs.index,counts@numTs.index))

    select_all_info<-
      c(rbind(counts@coverage.index ,
              counts@numCs.index,counts@numTs.index))

    max_mc<-max(c(mBase@coverage.index,mBase@numCs.index,mBase@numTs.index))
    slotList <- list(dbtype = mBase@dbtype,
                     sample.ids = c(mBase@sample.ids,"sigma","resid","Ameans"),assembly = mBase@assembly,
                     context = mBase@context,resolution = mBase@resolution,
                     treatment = c(mBase@treatment,NA,NA,NA),destranded = mBase@destranded,
                     coverage.index = c(mBase@coverage.index,(max_mc)+1,(max_mc)+4,(max_mc)+7),
                     numCs.index = c(mBase@numCs.index,max_mc+2,max_mc+5,(max_mc)+8),
                     numTs.index = c(mBase@numTs.index,max_mc+3,max_mc+6,(max_mc)+9))
    #slotList$sample.ids
    tabixHead <- methylKit:::makeTabixHeader(slotList)
    tabixHeadString <- methylKit:::.formatTabixHeader(class = "methylBaseDB",
                                                      tabixHead = tabixHead)

    file_name <- paste(tools::file_path_sans_ext(basename(methylKit::getDBPath(counts))),"_mixedcirc","_stats","",sep = "")
    tbxFile_t<-methylKit::getDBPath(y_norm)



    y<-methylKit:::applyTbxByChunk(tbxFile = tbxFile_t,tabixHead = tabixHeadString,
                                   chunk.size = chunk.size,return.type = "tabix",dir = dir_for_file,filename = file_name,
                                   FUN = function(x){
                                     xx<-as.matrix(x[,-c(1:4,counts@coverage.index),drop=F]);
                                     Ameans<-rowMeans(xx,na.rm = TRUE)[,drop=F];


                                     xx_res<-x;
                                     vpList = variance_fit(xx, formula, data, showWarnings = FALSE,
                                                                                  fxn = function(fit) {
                                                                                    list(sd = attr(lme4::VarCorr(fit), "sc"), fitted.values = predict(fit))
                                                                                  }, BPPARAM = BPPARAM,ignore_na = ignore_missing);
                                     fitted.values <- lapply(vpList, function(x) x$fitted.values);
                                     fitted.values <- do.call("rbind", fitted.values);
                                     fit = list();
                                     fit$sigma <- sapply(vpList, function(x) x$sd);
                                     fit$df.residual = rep(2, length(fit$sigma));
                                     xx_res[,-c(1:4,counts@coverage.index)]<-fitted.values;
                                     cbind(xx_res,fit$sigma,fit$sigma,fit$sigma,fit$df.residual,fit$df.residual,fit$df.residual,
                                           Ameans,Ameans,Ameans)

                                   })

    y_norm <- methylKit:::readMethylBaseDB(dbpath = y,
                                           sample.ids = slotList$sample.ids,
                                           assembly = counts@assembly,
                                           dbtype = "tabix",
                                           resolution = counts@resolution,
                                           context = counts@context,
                                           treatment = slotList$treatment,destranded = F)


    tbxFile_t<-methylKit::getDBPath(y_norm)
    sig_ameans<-methylKit:::applyTbxByChunk(tbxFile = tbxFile_t,tabixHead = tabixHeadString,
                                            chunk.size = chunk.size,return.type = "data.frame",dir = dir_for_file,filename = file_name,
                                            FUN = function(x){


                                              cbind(x[,y_norm@numCs.index[y_norm@sample.ids%in%c("sigma","Ameans")],drop=F],
                                                    rowSums( x[,c(y_norm@numCs.index[!y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )],
                                                                  y_norm@numTs.index[!y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )]),drop=F],na.rm = TRUE))


                                            })

    sx <- sig_ameans[,2] + mean(log2(lib.size + 1)) - log2(1e+06)

    sy <- sqrt(sig_ameans[,1])

    allzero <- sig_ameans[,3] == 0
    if (any(allzero)) {
      sx <- sx[!allzero]
      sy <- sy[!allzero]
    }
    l <- stats::lowess(sx, sy, f = span)
    if (plot) {
      plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
           pch = 16, cex = 0.25)
      title("voom: Mean-variance trend")
      lines(l, col = "red")
    }
    suppressWarnings({
      f <- approxfun(l, rule = 2)
    })

    ## estimate weights

    file_name <- paste(tools::file_path_sans_ext(basename(methylKit::getDBPath(counts))),"_mixedcirc","_weights","",sep = "")
    slotList_weights <- list(dbtype = mBase@dbtype,
                             sample.ids = c(mBase@sample.ids),assembly = mBase@assembly,
                             context = mBase@context,resolution = mBase@resolution,
                             treatment = c(mBase@treatment),destranded = mBase@destranded,
                             coverage.index = c(mBase@coverage.index),
                             numCs.index = c(mBase@numCs.index),
                             numTs.index = c(mBase@numTs.index))
    #slotList$sample.ids
    tabixHead_weights <- methylKit:::makeTabixHeader(slotList_weights)
    tabixHeadString_weights <- methylKit:::.formatTabixHeader(class = "methylBaseDB",
                                                              tabixHead = tabixHead_weights)


    tbxFile_t<-methylKit::getDBPath(y_norm)
    weights_norm<-methylKit:::applyTbxByChunk(tbxFile = tbxFile_t,tabixHead = tabixHeadString_weights,
                                              chunk.size = chunk.size,return.type = "tabix",dir = dir_for_file,filename = file_name,
                                              FUN = function(x){

                                                fitted.cpm<-2^x[,-c(1:4,y_norm@coverage.index,y_norm@numCs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans")],
                                                                    y_norm@numTs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )])]
                                                fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
                                                fitted.logcount <- log2(fitted.count)
                                                w <- 1/f(fitted.logcount)^4
                                                dim(w) <- dim(fitted.logcount)
                                                x[,-c(1:4,y_norm@coverage.index,y_norm@numCs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans")],
                                                      y_norm@numTs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )])]<-w
                                                x[,-c(y_norm@numCs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans")],
                                                      y_norm@numTs.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )],
                                                      y_norm@coverage.index[y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans")])]

                                              })

    weights_final <- methylKit:::readMethylBaseDB(dbpath = weights_norm,
                                                  sample.ids = slotList_weights$sample.ids,
                                                  assembly = counts@assembly,
                                                  dbtype = "tabix",
                                                  resolution = counts@resolution,
                                                  context = counts@context,
                                                  treatment = slotList_weights$treatment,destranded = F)


    y_data_norm<-paste(dir_for_file,"/",tools::file_path_sans_ext(basename(methylKit::getDBPath(counts))),"_mixedcirc","_libnorm",".bgz",sep = "")

    slotList <- list(dbtype = mBase@dbtype,
                     sample.ids = mBase@sample.ids,assembly = mBase@assembly,
                     context = mBase@context,resolution = mBase@resolution,
                     treatment = mBase@treatment,destranded = mBase@destranded,
                     coverage.index = mBase@coverage.index,
                     numCs.index = mBase@numCs.index,
                     numTs.index = mBase@numTs.index)
    y_norm <- methylKit:::readMethylBaseDB(dbpath = y_data_norm  ,
                                           sample.ids = slotList$sample.ids,
                                           assembly = counts@assembly,
                                           dbtype = "tabix",
                                           resolution = counts@resolution,
                                           context = counts@context,
                                           treatment = slotList$treatment,destranded = F)

    return(list(E=y_norm,weights=weights_final))

  }


}


#' This function is used to calculate sigma of fit.
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
variance_fit<-function (exprObj, formula, data, REML = FALSE, useWeights = TRUE,
                        weightsMatrix = NULL, showWarnings = TRUE, fxn = identity,
                        colinearityCutoff = 0.999, control = lme4::lmerControl(calc.derivs = FALSE,
                                                                               check.rankX = "stop.deficient"), quiet = FALSE, BPPARAM = bpparam(), ignore_na=FALSE,
                        ...)
{
  formula = stats::as.formula(formula)
  if (ncol(exprObj) != nrow(data)) {
    stop("the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)")
  }
  if (!is(exprObj, "sparseMatrix")) {
    countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
    if (countNA > 0) {
      if(!ignore_na)
      {
        stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
      }else{

        warning("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
      }

    }
    rv = apply(exprObj, 1, var,na.rm=ignore_na)
  }
  else {
    rv = c()
    for (i in seq_len(nrow(exprObj))) {
      rv[i] = var(exprObj[i, ],na.rm=ignore_na)
    }
  }
  if (any(rv == 0)) {
    idx = which(rv == 0)
    stop(paste("Response variable", idx[1], "has a variance of 0"))
  }
  if (useWeights && is.null(weightsMatrix)) {
    useWeights = FALSE
  }
  if (useWeights && !identical(dim(exprObj), dim(weightsMatrix))) {
    stop("exprObj and weightsMatrix must be the same dimensions")
  }
  if (.isDisconnected()) {
    stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
  }
  if (!identical(colnames(exprObj), rownames(data))) {
    warning("Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently.")
  }
  form = paste("responsePlaceholder$E", paste(as.character(formula),
                                              collapse = ""))
  responsePlaceholder = iterators::nextElem(exprIter(exprObj, weightsMatrix,
                                                                         useWeights))
  possibleError <- tryCatch(lmer(eval(parse(text = form)),
                                 data = data, ..., control = control), error = function(e) e)
  if (inherits(possibleError, "error")) {
    err = grep("object '.*' not found", possibleError$message)
    if (length(err) > 0) {
      stop("Variable in formula is not found: ", gsub("object '(.*)' not found",
                                                      "\\1", possibleError$message))
    }
  }
  pb <- progress::progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",
                                   , total = nrow(exprObj), width = 60, clear = FALSE)
  timeStart = proc.time()
  mesg <- "No random effects terms specified in formula"
  method = ""
  if (inherits(possibleError, "error") && identical(possibleError$message,
                                                    mesg)) {
    fit <- lm(eval(parse(text = form)), data = data, ...)
    checkModelStatus(fit, showWarnings = showWarnings, colinearityCutoff = colinearityCutoff)
    res <- foreach(responsePlaceholder = exprIter(exprObj,
                                                                      weightsMatrix, useWeights), .packages = c("splines",
                                                                                                                "lme4")) %do% {
                                                                                                                  fit = lm(eval(parse(text = form)), data = data, weights = responsePlaceholder$weights,
                                                                                                                           na.action = stats::na.exclude, ...)
                                                                                                                  fxn(fit)
                                                                                                                }
    method = "lm"
  }
  else {
    if (inherits(possibleError, "error") && grep("the fixed-effects model matrix is column rank deficient",
                                                 possibleError$message) == 1) {
      stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
    }
    responsePlaceholder = iterators:::nextElem(exprIter(exprObj, weightsMatrix,
                                                                            useWeights))
    timeStart = proc.time()
    fitInit <- lmer(eval(parse(text = form)), data = data,
                    ..., REML = REML, control = control)
    timediff = proc.time() - timeStart
    objSize = object.size(fxn(fitInit)) * nrow(exprObj)
    if (!quiet)
      message("Memory usage to store result: >", format(objSize,
                                                        units = "auto"))
    checkModelStatus(fitInit, showWarnings = showWarnings,
                                         colinearityCutoff = colinearityCutoff)
    data2 = data.frame(data, expr = responsePlaceholder$E,
                       check.names = FALSE)
    form = paste("expr", paste(as.character(formula), collapse = ""))
    .eval_models = function(responsePlaceholder, data2, form,
                            REML, theta, fxn, control, na.action = stats::na.exclude,
                            ...) {
      data2$expr = responsePlaceholder$E
      fit = lmer(eval(parse(text = form)), data = data2,
                 ..., REML = REML, weights = responsePlaceholder$weights,
                 control = control, na.action = na.action)
      fxn(fit)
    }
    .eval_master = function(obj, data2, form, REML, theta,
                            fxn, control, na.action = stats::na.exclude, ...) {
      lapply(seq_len(nrow(obj$E)), function(j) {
        .eval_models(list(E = obj$E[j, ], weights = obj$weights[j,
        ]), data2, form, REML, theta, fxn, control,
        na.action, ...)
      })
    }
    it = iterBatch(exprObj, weightsMatrix, useWeights, n_chunks = 100)
    if (!quiet)
      message(paste0("Dividing work into ", attr(it, "n_chunks"),
                     " chunks..."))
    res <- BiocParallel::bpiterate(it, .eval_master, data2 = data2, form = form,
                                   REML = REML, theta = fitInit@theta, fxn = fxn, control = control,
                                   ..., REDUCE = c, reduce.in.order = TRUE, BPPARAM = BPPARAM)
    if (is(res, "remote_error")) {
      stop("Error evaluating fxn:\n\n", res)
    }
    method = "lmer"
  }
  if (!quiet)
    message("\nTotal:", paste(format((proc.time() - timeStart)[3],
                                     digits = 0), "s"))
  names(res) <- rownames(exprObj)
  new("VarParFitList", res, method = method)
}






# Check that lm/lmer model is valid
# Throw warning if
#	1) Intercept is ommited
#	2) Any coefficient is NA
#	3) a categorical variable is modeled as a fixed effect
setGeneric("checkModelStatus", signature="fit",
           function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE )
             standardGeneric("checkModelStatus")
)

setMethod("checkModelStatus", "lm",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE  )
          {
            # if no intercept is specified, give warning
            if( showWarnings && length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
              txt = "No Intercept term was specified in the formula"
              stop(txt)
              # warning(txt, immediate.=immediate)
            }

            # if any coefficient is NA
            if( showWarnings && any(is.na(coef(fit))) ){
              stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
            }

            # check colinearity
            score = colinearityScore(fit)
            if( score > colinearityCutoff ){
              stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
            }
          }
)

setMethod("checkModelStatus", "lmerMod",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE  ){
            run_model_check_mixed( fit, showWarnings, dream, colinearityCutoff, immediate )
          })


setMethod("checkModelStatus", "glmerMod",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE  ){
            run_model_check_mixed( fit, showWarnings, dream, colinearityCutoff, immediate )
          })

#' @importFrom aod negbin
setMethod("checkModelStatus", "negbin",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE  ){
            # run_model_check( fit, showWarnings, dream, colinearityCutoff, immediate )
          })

#' @importFrom lme4 isSingular
run_model_check_mixed = function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999, immediate=FALSE  ){
  # if no intercept is specified, give warning
  if( !dream && showWarnings && length(which(colnames(fit@pp$X) == "(Intercept)")) == 0 ){
    txt = "No Intercept term was specified in the formula."
    stop(txt)
    # warning(txt, immediate.=immediate)
  }

  # if any coefficient is NA
  if( ( showWarnings | dream) && any(is.na(coef(fit))) ){
    stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
  }

  # # check colinearity
  # ###################
  # score = colinearityScore(fit)

  # if( score > colinearityCutoff ){
  # 	stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
  # }

  # Check condition number of covariance matrix
  condNum = kappa(cov2cor(as.matrix(vcov(fit)))) # exact=TRUE

  if( condNum > 1e8 ){
    stop(paste0("Condition number (", format(condNum, digits=1), ") is very high.\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
  }

  # check that factors are random and continuous variables are fixed
  ###################################################################

  # remove backticks with gsub manually
  # solve issue that backticks are conserved is some but not all parts of lmer()

  # Simplified testing of random versus fixed effects
  # allows (A|B) only where A is continuous

  # variables fit by regression
  testVar = attr(attr(fit@frame, "terms"), "term.labels")
  testVar = gsub("`", "", testVar)

  # get type for each variable
  # keep only tested variables
  varType = attr(attr(fit@frame, "terms"), "dataClasses")[-1]
  varType = varType[testVar]

  # random effects
  randVar = names(fit@flist)

  # fixed effects
  # starting with all variables, remove random variables
  fixedVar = setdiff(testVar, randVar)

  for( i in 1:length(varType) ){

    # if factor is not random
    if( (showWarnings && ! dream) && varType[i] %in% c("factor", "character") && (! names(varType)[i] %in% randVar) ){
      txt = paste("Categorical variables modeled as fixed effect:", paste(names(varType)[i], collapse=', '), "\nMust model either _all_ or _no_ categorical variables as random effects here")
      stop(txt)
    }

    # If numeric/double is not fixed
    if( (showWarnings && ! dream) && varType[i] %in% c("numeric", "double") && (!names(varType)[i] %in% fixedVar) ){
      stop(paste("Continuous variable cannot be modeled as a random effect:", names(varType)[i]))
    }
  }

  # show convergance message, if model is not singular
  if( showWarnings && !is.null(fit@optinfo$conv$lme4$messages) && ! isSingular(fit)){
    stop(fit@optinfo$conv$lme4$messages)
  }
}







.isMixedModelFormula <- function(formula) {
  !is.null(findbars(as.formula(formula)))
}

isDisconnected<-function (){
  i = NULL
  possibleError <- tryCatch(suppressWarnings(foreach(i = seq_len(2)) %dopar%
                                               {
                                                 i
                                               }), error = function(e) e)
  return(isTRUE(inherits(possibleError, "error") && identical(possibleError$message,
                                                              "invalid connection")))
}


# requared to that iterator return NULL after the last element
icount2 = function (count){
  if (missing(count))
    count <- NULL
  else if (!is.numeric(count) || length(count) != 1)
    stop("count must be a numeric value")
  i <- 0L
  nextEl <- function(){
    if( is.null(i) )
      (i <<- NULL)
    else if (is.null(count) || i < count)
      (i <<- i + 1L)
    else
      (i <<- NULL)
  }
  it <- list(nextElem = nextEl)
  class(it) <- c("abstractiter", "iter")
  it
}


# Iterator over genes
#' @importFrom iterators icount
exprIter = function( exprObj, weights, useWeights = TRUE, scale=TRUE, iterCount = "icount"){

  n_features = nrow(exprObj)

  if( iterCount == 'icount2'){
    xit <- icount2( n_features )
  }else{
    xit <- icount( n_features )
  }

  nextEl <- function() {
    j <- nextElem(xit)

    if( is.null(j) || j > n_features){
      res = NULL
    }else{
      if( useWeights && !is.null(weights) ){

        w = weights[j,]

        # scale weights to have mean of 1, otherwise it affects the residual variance too much
        # scale should be false when signa(fit) needs to be evaluted
        if(scale){
          w = w / mean(w)
        }
      }else{
        w = NULL
      }

      res = list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
    }
    res
  }
  it <- list(nextElem = nextEl)
  class(it) <- c("abstractiter", "iter")
  it
}



# exprIterOrig = function( exprObj, weights, useWeights = TRUE, scale=TRUE){

# 	n_features = nrow(exprObj)
# 	xit <- icountn( n_features )

#     nextEl <- function() {
#     	j <- nextElem(xit)

#     	if( useWeights && !is.null(weights) ){
# 			# scale weights to have mean of 1, otherwise it affects the residual variance too much
#     		if(scale){
#     			w = weights[j,] /  mean(weights[j,])
#     		}else{
#     			w = weights[j,]
#     		}
#     	}else{
#     		w = NULL
# 		}

#        	list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
#     }
#     it <- list(nextElem = nextEl)
#     class(it) <- c("abstractiter", "iter")
#     it
# }


#' @importFrom BiocParallel bpworkers
iterBatch <- function(exprObj, weights, useWeights = TRUE, scale=TRUE, n_chunks = nrow(exprObj) / 500, min_chunk_size = 20, BPPARAM = NULL ) {
  # Adjust number of chunks upward to the next multiple of number of
  # workers in BPPARAM, if this can be determined. If any errors are
  # encountered, just continue without adjusting.
  tryCatch(
    if (is(BPPARAM, "BiocParallelParam")) {
      n_workers <- bpworkers(BPPARAM)
      if (!is.null(n_workers) && is.numeric(n_workers) && n_workers >= 1) {
        chunks_per_worker <- ceiling(n_chunks / n_workers)
        n_chunks <- chunks_per_worker * n_workers
      }
    },
    error = function(...) NULL
  )

  # Don't split into chunks smaller than min_chunk_size
  max_allowed_chunks <- floor(nrow(exprObj) / min_chunk_size)
  n_chunks = min(n_chunks, max_allowed_chunks)
  # Make sure we have at least 1 chunk (since we can get 0 if
  # min_chunk_size > nrow)
  n_chunks <- max(n_chunks, 1)

  # specify chunks
  idx <- parallel::splitIndices(nrow(exprObj), min(nrow(exprObj), n_chunks))
  i <- 0L

  f = function() {
    if (i == length(idx)){
      return(NULL)
    }
    i <<- i + 1L
    E = exprObj[ idx[[i]],, drop = FALSE ]

    if( useWeights && !is.null(weights) ){
      # scale weights to have mean of 1, otherwise it affects the residual variance too much
      if(scale){
        # w = weights[j,] /  mean(weights[j,])
        w = weights[idx[[i]],,drop=FALSE]
        # for each row, devide by row mean
        w = w / rowMeans(w)
      }else{
        # w = weights[j,]
        w = weights[idx[[i]],,drop=FALSE]
      }
    }else{
      w = matrix(1, nrow(E), ncol(E))
    }

    list(E = E, weights = w )
  }

  # get number of chunks
  attr( f, "n_chunks") = length(idx)
  f
}



#' Transform RNA-Seq Data Ready for Linear Mixed Modelling with \code{dream()}
#'
#' Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear mixed modelling with \code{dream()}.   This method is the same as \code{limma::voom()}, except that it allows random effects in the formula
#'
#' @param counts a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object. Counts must be non-negative and NAs are not permitted.
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work, and \code{lmFit()} followed by contrasts.fit() are run.
#' @param data \code{data.frame} with columns corresponding to formula
#' @param lib.size numeric vector containing total library sizes for each sample.  Defaults to the normalized (effective) library sizes in \code{counts} if \code{counts} is a \code{DGEList} or to the columnwise count totals if \code{counts} is a matrix.
#' @param normalize.method the microarray-style normalization method to be applied to the logCPM values (if any).  Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.  Any normalization factors found in \code{counts} will still be used even if \code{normalize.method="none"}.
#' @param span width of the lowess smoothing window as a proportion.
#' @param weights Can be a numeric matrix of individual weights of same dimensions as the \code{counts}, or a numeric vector of sample weights with length equal to \code{ncol(counts)}
#' @param plot logical, should a plot of the mean-variance trend be displayed?
#' @param save.plot logical, should the coordinates and line of the plot be saved in the output?
#' @param quiet suppress message, default FALSE
#' @param BPPARAM parameters for parallel evaluation
#' @param      ... other arguments are passed to \code{lmer}.
#'
#' @return
#' An \code{EList} object just like the result of \code{limma::voom()}
#'
#' @details Adapted from \code{vomm()} in \code{limma} v3.40.2
#' @seealso \code{limma::voom()}
#' @examples
#' # library(variancePartition)
#' library(edgeR)
#' library(BiocParallel)
#'
#' data(varPartDEdata)
#'
#' # normalize RNA-seq counts
#' dge = DGEList(counts = countMatrix)
#' dge = calcNormFactors(dge)
#'
#' # specify formula with random effect for Individual
#' form <- ~ Disease + (1|Individual)
#'
#' # compute observation weights
#' vobj = voomWithDreamWeights( dge[1:20,], form, metadata)
#'
#' # fit dream model
#' res = dream( vobj, form, metadata)
#' res = eBayes(res)
#'
#' # extract results
#' topTable(res, coef="Disease1", number=3)
#'
#' @importFrom lme4 VarCorr
#' @importFrom stats approxfun predict as.formula
#' @importFrom limma asMatrixWeights
#' @export
voomWithDreamWeights <- function(counts, formula, data, lib.size=NULL, normalize.method="none", span=0.5, weights = NULL, plot=FALSE, save.plot=FALSE, quiet=FALSE, BPPARAM=SerialParam(),...){

  formula = as.formula( formula )

  # only retain columns used in the formula
  data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]

  # check if variables in formula has NA's
  hasNA = hasMissingData(formula, data)

  if( any(hasNA) ){
    txt = paste("Variables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\n  Missing data is not supported")
    stop( txt )
  }

  out <- list()

  design = NULL

  #	Check counts
  if(is(counts,"DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    # if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
    if(is.null(lib.size)) lib.size <- counts$samples$lib.size*counts$samples$norm.factors
    counts <- counts$counts
  } else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
    if(isExpressionSet) {
      if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
      if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    } else if( is.null(dim(counts)) ){
      stop("counts is type '", class(counts), "' and can't be converted to matrix unambiguously")
    } else {
      counts <- as.matrix(counts)
    }
  }

  n <- nrow(counts)
  if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")

  # #	Check design
  # if(is.null(design)) {
  # 	design <- matrix(1,ncol(counts),1)
  # 	rownames(design) <- colnames(counts)
  # 	colnames(design) <- "GrandMean"
  # }

  # Check lib.size
  if(is.null(lib.size)) lib.size <- colSums(counts)

  #	Fit linear model to log2-counts-per-million
  y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
  y <- normalizeBetweenArrays(y,method=normalize.method)

  # Fit regression model
  #---------------------

  # use pre-specified weights, if available

  if( .isMixedModelFormula( formula ) ){

    if( missing(data) ){
      stop("Must specify argument 'data'\n")
    }

    if( !is.null(weights) ){

      # if weights is a per-sample vector
      if( length(weights) == ncol(y) ){
        # convert weights vector to matrix
        weights = asMatrixWeights(weights, dim(y))
      }
    }

    # fit linear mixed model
    vpList = fitVarPartModel( y, formula, data, weightsMatrix=weights, showWarnings=FALSE, ...,fxn = function(fit){
      # extract
      # 1) sqrt residual variance (i.e. residual standard deviation)
      # 2) fitted values
      list( sd = attr(VarCorr(fit), 'sc'),
            fitted.values = predict(fit) )
    }, BPPARAM=BPPARAM, quiet=quiet )

    fit = list()
    fit$sigma <- sapply( vpList, function(x) x$sd)
    fit$df.residual = rep(2, length(fit$sigma)) # check this

    # extract fitted values
    fitted.values <- lapply( vpList, function(x) x$fitted.values)
    fitted.values <- do.call("rbind", fitted.values )

  }else{

    #if( ! quiet) message("Fixed effect model, using limma directly...")

    design = model.matrix(formula, data)
    fit <- lmFit(y,design,weights=weights,...)

    if(fit$rank < ncol(design)) {
      j <- fit$pivot[1:fit$rank]
      fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
    } else {
      fitted.values <- fit$coef %*% t(fit$design)
    }
  }

  if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

  #		If no replication found, set all weight to 1
  NWithReps <- sum(fit$df.residual > 0L)
  if(NWithReps < 2L) {
    if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
    if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    if( !is.null(design) ) out$design <- design
    if(is.null(out$targets))
      out$targets <- data.frame(lib.size=lib.size)
    else
      out$targets$lib.size <- lib.size
    return(new("EList",out))
  }

  # Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)

  # get residual standard deviation
  if( is(fit, "MArrayLM2") ){
    # fit is result of dream()
    sy <- sqrt(attr(fit, "varComp")$resid)
  }else{
    # fit is result of lmFit()
    sy <- sqrt(fit$sigma)
  }

  allzero <- rowSums(counts)==0
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- stats::lowess(sx,sy,f=span)
  if(plot) {
    plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
    title("voom: Mean-variance trend")
    lines(l,col="red")
  }

  #	Make interpolating rule
  #	Special treatment of zero counts is now removed;
  #	instead zero counts get same variance as smallest gene average.
  #	l$x <- c(0.5^0.25, l$x)
  #	l$x <- c(log2(0.5), l$x)
  #	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
  #	var0 <- max(var0,1e-6)
  #	l$y <- c(var0, l$y)
  suppressWarnings({
    f <- approxfun(l, rule=2)
  })

  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
  fitted.logcount <- log2(fitted.count)

  #	Apply trend to individual observations
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)

  #	Output
  out$E <- y
  out$weights <- w
  if( !is.null(design) ) out$design <- design
  if(is.null(out$targets))
    out$targets <- data.frame(lib.size=lib.size)
  else
    out$targets$lib.size <- lib.size

  if(save.plot) {
    out$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )")
    out$voom.line <- l
  }

  # Check max value of precision weights
  maxValue = max(out$weights)
  if( maxValue > 1e8){
    txt = paste0("The maximum precision weight is ", format(maxValue, scientific=TRUE), ", suggesting a poor smoothing fit\non the mean-variance plot for large expression values. Such large weights can\nhave unexpected effects downstream.  Consider examining the mean-variance plot\nand reducing the span parameter.")
    warning(txt)
  }

  new("EList",out)
}


## default control for lmer fits
vpcontrol <- lme4::lmerControl(calc.derivs = FALSE,
                               check.rankX = "stop.deficient",
                               check.conv.singular =
                                 lme4::.makeCC("ignore", tol = 1e-4))


#' Fit linear (mixed) model
#'
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  See Details.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from voom().  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from \code{voom()}.  Used only if useWeights is TRUE
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param fxn apply function to model fit for each gene.  Defaults to identify function so it returns the model fit itself
#' @param control control settings for \code{lmer()}
#' @param BPPARAM parameters for parallel evaluation
#' @param quiet suppress message, default FALSE
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#'
#' @return
#' \code{list()} of where each entry is a model fit produced by \code{lmer()} or \code{lm()}
#'
#' @importFrom MASS ginv
# @importFrom RSpectra eigs_sym
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom graphics abline axis hist image layout lines mtext par plot plot.new rect text title
#' @importFrom stats anova as.dendrogram as.dist cancor coef cov2cor density dist fitted.values hclust lm median model.matrix order.dendrogram quantile reorder residuals sd terms var vcov pt qt
#' @importFrom scales rescale
#' @importFrom iterators nextElem
#' @import Rdpack
#'
#' @details
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is
#'
#' \code{fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)}
#'
#' If there are no random effects, so formula is \code{~ a + b + c}, a 'standard' linear model is used:
#'
#' \code{fit <- lm( exprObj[j,] ~ a + b + c, data=data)}
#'
#' In both cases, \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run in parallel.  Parallel processing must be enabled before calling this function.  See below.
#'
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lm/lmer.
#'
#' Since this function returns a list of each model fit, using this function is slower and uses more memory than \code{fitExtractVarPartModel()}.
#'
#' \code{REML=FALSE} uses maximum likelihood to estimate variance fractions.  This approach produced unbiased estimates, while \code{REML=TRUE} can show substantial bias.  See Vignette "3) Theory and practice of random effects and REML"
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' library(BiocParallel)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue)
#'
#' # Step 1: fit linear mixed model on gene expression
#' # If categorical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'
#' # violin plot of contribution of each variable to total variance
#' # also sort columns
#' plotVarPart( sortCols( varPart ) )
#'
#' # Advanced:
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results <- fitVarPartModel( geneExpr, form, info )
#'
#' # Step 2: extract variance fractions
#' varPart <- extractVarPart( results )
#'
#' # Note: fitVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form <- ~ (1|sex) + (1|type) + score
#' info2 <- Biobase::pData(sample.ExpressionSet)
#' results2 <- fitVarPartModel( sample.ExpressionSet, form, info2 )
#'
# # Parallel processing using multiple cores with reduced memory usage
# param <- SnowParam(4, "SOCK", progressbar=TRUE)
# results2 <- fitVarPartModel( sample.ExpressionSet, form, info2, BPPARAM=param)
#'
#' @export
#' @docType methods
#' @rdname fitVarPartModel-method
setGeneric("fitVarPartModel", signature="exprObj",
           function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(),...)
             standardGeneric("fitVarPartModel")
)

# internal driver function
#' @importFrom BiocParallel SerialParam bpiterate bplapply bpok
#' @importFrom methods is new
#' @importFrom lme4 lmer
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom progress progress_bar
#' @importFrom utils object.size
#' @import foreach
#' @import lme4
.fitVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, colinearityCutoff=.999,control = vpcontrol, quiet=quiet, BPPARAM=SerialParam(), ...){

  # convert to data.frame
  data = as.data.frame(data)

  # if( ! is(exprObj, "sparseMatrix")){
  # 	exprObj = as.matrix( exprObj )
  # }
  formula = stats::as.formula( formula )

  # only retain columns used in the formula
  data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]
  data = droplevels(data)

  # check that variables in the formula are all in the data
  idx = unique(all.vars(formula)) %in% colnames(data)
  if( any(!idx) ){
    txt = paste(unique(all.vars(formula))[!idx], collapse=', ')
    stop("Variable in formula not found in data: ", txt)
  }

  # check dimensions of reponse and covariates
  if( ncol(exprObj) != nrow(data) ){
    stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
  }

  # check if variables in formula has NA's
  hasNA = hasMissingData(formula, data)

  if( any(hasNA) ){
    warning(paste("Variables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\nSamples with missing data will be dropped.\n"), immediate.=TRUE, call.=FALSE)

    # drop samples with missing data in formula variables
    idx = sapply(all.vars(formula), function(v) {
      which(is.na(data[[v]]))
    })
    idx = unique(unlist(idx))

    data = droplevels(data[-idx,,drop=FALSE])
    exprObj = exprObj[,-idx,drop=FALSE]
    exprObjMat = as.matrix( exprObj )
  }

  # check if all genes have variance
  if( ! is(exprObj, "sparseMatrix")){
    # check if values are NA
    countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
    if( countNA > 0 ){
      stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
    }

    rv = apply( exprObj, 1, var)
  }else{
    rv = c()
    for( i in seq_len(nrow(exprObj)) ){
      rv[i] = var( exprObj[i,])
    }
  }
  if( any( rv == 0) ){
    idx = which(rv == 0)
    stop(paste("Response variable", idx[1], 'has a variance of 0'))
  }

  # if weightsMatrix is not specified, set useWeights to FALSE
  if( useWeights && is.null(weightsMatrix) ){
    # warning("useWeights was ignored: no weightsMatrix was specified")
    useWeights = FALSE
  }

  # if useWeights, and (weights and expression are the same size)
  if( useWeights && !identical( dim(exprObj), dim(weightsMatrix)) ){
    stop( "exprObj and weightsMatrix must be the same dimensions" )
  }
  if( .isDisconnected() ){
    stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
  }

  # If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
  if( ! identical(colnames(exprObj), rownames(data)) ){
    warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
  }

  # add response (i.e. exprObj[j,]) to formula
  # Use term 'responsePlaceholder' to store the value of reponse j (exprObj[j,])
  # This is not an ideal way to structure the evaluation of response j.
  # 	The formula is evaluated a different scope, (i.e. within lmer()), and there is no need to pass the
  # 	entire exprObj object into that scope.  With lexical scope, I found it was possible that
  # 	the value of exprObj[j,] could be different when evaluated in the lower scope
  # This placeholder term addresses that issue
  form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

  # run lmer() to see if the model has random effects
  # if less run lmer() in the loop
  # else run lm()
  responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
  possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,...,control=control ), error = function(e) e)

  # detect error when variable in formula does not exist
  if( inherits(possibleError, "error") ){
    err = grep("object '.*' not found", possibleError$message)
    if( length(err) > 0 ){
      stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
    }
  }

  pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",,
                         total = nrow(exprObj), width= 60, clear=FALSE)

  timeStart = proc.time()

  mesg <- "No random effects terms specified in formula"
  method = ''
  if( isTRUE(inherits(possibleError, "error") && identical(possibleError$message, mesg)) ){

    # fit the model for testing
    fit <- lm( eval(parse(text=form)), data=data,...)

    # check that model fit is valid, and throw warning if not
    checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

    resList <- foreach(responsePlaceholder=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4") ) %do% {
      # fit linear mixed model
      fit = lm( eval(parse(text=form)), data=data, weights=responsePlaceholder$weights,na.action=stats::na.exclude,...)

      # apply function
      fxn( fit )
    }

    method = "lm"

  }else{

    if( isTRUE(inherits(possibleError, "error") &&  grep('the fixed-effects model matrix is column rank deficient', possibleError$message) == 1) ){
      stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
    }

    # fit first model to initialize other model fits
    # this make the other models converge faster
    responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

    timeStart = proc.time()
    fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )

    timediff = proc.time() - timeStart

    # check size of stored objects
    objSize = object.size( fxn(fitInit) ) * nrow(exprObj)

    if( !quiet ) message("Memory usage to store result: >", format(objSize, units = "auto"))

    # check that model fit is valid, and throw warning if not
    checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

    # specify gene explicitly in data
    # required for downstream processing with lmerTest
    data2 = data.frame(data, expr=responsePlaceholder$E, check.names=FALSE)
    form = paste( "expr", paste(as.character( formula), collapse=''))

    # Define function for parallel evaluation
    .eval_models = function(responsePlaceholder, data2, form, REML, theta, fxn, control, na.action=stats::na.exclude,...){

      # modify data2 for this gene
      data2$expr = responsePlaceholder$E

      # fit linear mixed model
      fit = lmer( eval(parse(text=form)), data=data2, ..., REML=REML, weights=responsePlaceholder$weights, control=control,na.action=na.action)

      # apply function
      fxn( fit )
    }

    .eval_master = function( obj, data2, form, REML, theta, fxn, control, na.action=stats::na.exclude,... ){

      # use only 1 OpenMP thread for linear algebra
      omp_set_num_threads(1)

      lapply(seq_len(nrow(obj$E)), function(j){
        .eval_models( list(E=obj$E[j,], weights=obj$weights[j,]), data2, form, REML, theta, fxn, control, na.action,...)
      })
    }

    # Evaluate function
    ###################

    it = iterBatch(exprObj, weightsMatrix, useWeights, n_chunks = 100, BPPARAM = BPPARAM)

    if( !quiet ) message(paste0("Dividing work into ",attr(it, "n_chunks")," chunks..."))

    resList <- bpiterate( it, .eval_master,
                          data2=data2, form=form, REML=REML, theta=fitInit@theta, fxn=fxn, control=control,...,
                          BPPARAM=BPPARAM)

    # if there is an error in evaluating fxn (usually in parallel backend)
    if( !bpok(list(resList)) ){
      stop("Error evaluating fxn:\n\n", resList)
    }
    # It can also return a list of errors, or a list where only some elements are errors
    if( !all(bpok(resList)) ){
      first_error <- resList[[which(!bpok(resList))[1]]]
      stop("Error evaluating fxn:\n\n", first_error)
    }

    # If no errors, then it's safe to concatenate all the results together.
    resList <- do.call(c, resList)

    method = "lmer"
  }

  # pb$update( responsePlaceholder$max_iter / responsePlaceholder$max_iter )
  if( !quiet ) message("\nTotal:", paste(format((proc.time() - timeStart)[3], digits = 1, scientific = FALSE), "s"))
  # set name of each entry
  names(resList) <- rownames( exprObj )

  new( "VarParFitList", resList, method=method )
}

## matrix
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,matrix-method
setMethod("fitVarPartModel", "matrix",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitVarPartModel(exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM,...)
          }
)

# data.frame
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,data.frame-method
setMethod("fitVarPartModel", "data.frame",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitVarPartModel( as.matrix(exprObj), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
          }
)

## EList
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,EList-method
setMethod("fitVarPartModel", "EList",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitVarPartModel( as.matrix(exprObj$E), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM,...)
          }
)

## ExpressionSet
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,ExpressionSet-method
#' @importFrom Biobase ExpressionSet exprs
setMethod("fitVarPartModel", "ExpressionSet",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitVarPartModel( as.matrix(exprs(exprObj)), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
          }
)

# sparseMatrix
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,sparseMatrix-method
#' @importFrom Matrix sparseMatrix
setMethod("fitVarPartModel", "sparseMatrix",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitVarPartModel( exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
          }
)

#' Fit linear (mixed) model, report variance fractions
#'
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Report fraction of variance attributable to each variable
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.   See Details.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from \code{voom()}.  Used only if \code{useWeights} is TRUE
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param control control settings for \code{lmer()}
#' @param quiet suppress message, default FALSE
#' @param BPPARAM parameters for parallel evaluation
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#'
#' @return
#' list() of where each entry is a model fit produced by \code{lmer()} or \code{lm()}
#'
#' @details
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' If there are no random effects, so formula is ~ a + b + c, a 'standard' linear model is used:
#'
#' \code{fit <- lm( exprObj[j,] ~ a + b + c, data=data)}
#'
#' In both cases, \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run in parallel.  Parallel processing must be enabled before calling this function.  See below.
#'
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to \code{lm}/\code{lmer}.
#'
#' \code{REML=FALSE} uses maximum likelihood to estimate variance fractions.  This approach produced unbiased estimates, while \code{REML=TRUE} can show substantial bias.  See Vignette "3) Theory and practice of random effects and REML"
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' library(BiocParallel)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue)
#'
#' # Step 1: fit linear mixed model on gene expression
#' # If categorical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#'
#' # Note: fitExtractVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form <- ~ (1|sex) + (1|type) + score
#' info2 <- Biobase::pData(sample.ExpressionSet)
#' varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2 )
#'
# # Parallel processing using multiple cores with reduced memory usage
# param = SnowParam(4, "SOCK", progressbar=TRUE)
# varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2, BPPARAM = param)
#'
#'
#' @export
#' @docType methods
#' @rdname fitExtractVarPartModel-method
#' @importFrom BiocParallel SerialParam bpiterate bplapply bpok
setGeneric("fitExtractVarPartModel", signature="exprObj",
           function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
             standardGeneric("fitExtractVarPartModel")
)

# internal driver function
#' @importFrom methods is new
#' @importFrom lme4 lmer
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom progress progress_bar
#' @import foreach
.fitExtractVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE, colinearityCutoff=.999, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(),...){

  # convert to data.frame
  data = as.data.frame(data)

  # exprObj = as.matrix( exprObj )
  formula = stats::as.formula( formula )

  # only retain columns used in the formula
  data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]
  data = droplevels(data)

  # check that variables in the formula are all in the data
  idx = unique(all.vars(formula)) %in% colnames(data)
  if( any(!idx) ){
    txt = paste(unique(all.vars(formula))[!idx], collapse=', ')
    stop("Variable in formula not found in data: ", txt)
  }

  # check dimensions of reponse and covariates
  if( ncol(exprObj) != nrow(data) ){
    stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
  }

  # check if variables in formula has NA's
  hasNA = hasMissingData(formula, data)

  if( any(hasNA) ){
    warning(paste("Variables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\nSamples with missing data will be dropped.\n"), immediate.=TRUE, call.=FALSE)

    # drop samples with missing data in formula variables
    idx = sapply(all.vars(formula), function(v) {
      which(is.na(data[[v]]))
    })
    idx = unique(unlist(idx))

    data = droplevels(data[-idx,,drop=FALSE])
    exprObj = exprObj[,-idx,drop=FALSE]
    exprObjMat = as.matrix( exprObj )
  }

  if( ! is(exprObj, "sparseMatrix")){
    # check if values are NA
    countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
    if( countNA > 0 ){
      stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
    }

    rv = apply( exprObj, 1, var)
  }else{
    # if exprObj is a sparseMatrix, this method will compute row-wise
    # variances with using additional memory
    rv = c()
    for( i in seq_len(nrow(exprObj)) ){
      rv[i] = var( exprObj[i,])
    }
  }
  if( any( rv == 0) ){
    idx = which(rv == 0)
    stop(paste("Response variable", idx[1], 'has a variance of 0'))
  }

  # if weightsMatrix is not specified, set useWeights to FALSE
  if( useWeights && is.null(weightsMatrix) ){
    # warning("useWeights was ignored: no weightsMatrix was specified")
    useWeights = FALSE
  }

  # if useWeights, and (weights and expression are the same size)
  if( useWeights && !identical( dim(exprObj), dim(weightsMatrix)) ){
    stop( "exprObj and weightsMatrix must be the same dimensions" )
  }
  if( .isDisconnected() ){
    stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
  }

  # add response (i.e. exprObj[,j] to formula
  form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

  # run lmer() to see if the model has random effects
  # if yes run lmer() in the loop
  # else run lm()
  responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
  possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data, control=control,... ), error = function(e) e)

  # detect error when variable in formula does not exist
  if( inherits(possibleError, "error") ){
    err = grep("object '.*' not found", possibleError$message)
    if( length(err) > 0 ){
      stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
    }
  }

  if( !quiet) pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",
                                     total = nrow(exprObj), width= 60, clear=FALSE)

  if( ! .isMixedModelFormula( formula ) ){

    # fit the model for testing
    fit <- lm( eval(parse(text=form)), data=data,...)

    # check that model fit is valid, and throw warning if not
    checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

    testValue = calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )

    timeStart = proc.time()

    varPart <- foreach(responsePlaceholder=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4") ) %do% {

      # fit linear mixed model
      fit = lm( eval(parse(text=form)), data=data, weights=responsePlaceholder$weights,na.action=stats::na.exclude,...)

      calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )
    }

    modelType = "anova"

  }else{

    if( inherits(possibleError, "error") && grep('the fixed-effects model matrix is column rank deficient', possibleError$message) == 1 ){
      stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
    }

    # fit first model to initialize other model fits
    # this make the other models converge faster
    responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

    timeStart = proc.time()
    fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control)
    timediff = proc.time() - timeStart

    # check that model fit is valid, and throw warning if not
    checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

    timeStart = proc.time()

    # Define function for parallel evaluation
    .eval_models = function(responsePlaceholder, data, form, REML, theta, control, na.action=stats::na.exclude,...){
      # fit linear mixed model
      fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=responsePlaceholder$weights, control=control,na.action=na.action)

      calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )
    }

    .eval_master = function( obj, data, form, REML, theta, control, na.action=stats::na.exclude,... ){

      # use only 1 OpenMP thread for linear algebra
      omp_set_num_threads(1)

      lapply(seq_len(nrow(obj$E)), function(j){
        .eval_models( list(E=obj$E[j,], weights=obj$weights[j,]), data, form, REML, theta, control, na.action,...)
      })
    }

    # Evaluate function
    ####################

    it = iterBatch(exprObj, weightsMatrix, useWeights, n_chunks = 100, BPPARAM = BPPARAM)

    if( !quiet) message(paste0("Dividing work into ",attr(it, "n_chunks")," chunks..."))

    varPart <- bpiterate( it, .eval_master,
                          data=data, form=form, REML=REML, theta=fitInit@theta, control=control,...,
                          BPPARAM=BPPARAM)

    # if there is an error in evaluating fxn (usually in parallel backend)
    if( !bpok(list(varPart)) ){
      stop("Error evaluating fxn:\n\n", varPart)
    }
    # It can also return a list of errors, or a list where only some elements are errors
    if( !all(bpok(varPart)) ){
      first_error <- varPart[[which(!bpok(varPart))[1]]]
      stop("Error evaluating fxn:\n\n", first_error)
    }

    # If no errors, then it's safe to concatenate all the results together.
    varPart <- do.call(c, varPart)

    modelType = "linear mixed model"
  }

  if(!quiet) message("\nTotal:", paste(format((proc.time() - timeStart)[3], digits = 1, scientific = FALSE), "s"))

  varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
  colnames(varPartMat) <- names(varPart[[1]])
  rownames(varPartMat) <- rownames(exprObj)

  res <- new("varPartResults", varPartMat, type=modelType, method="Variance explained (%)")

  return( res )
}

# matrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,matrix-method
setMethod("fitExtractVarPartModel", "matrix",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitExtractVarPartModel(exprObj, formula, data,
                                    REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                                    BPPARAM=BPPARAM, ...)
          }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod("fitExtractVarPartModel", "data.frame",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitExtractVarPartModel( as.matrix(exprObj), formula, data,
                                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                                     BPPARAM=BPPARAM, ...)
          }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod("fitExtractVarPartModel", "EList",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitExtractVarPartModel( as.matrix(exprObj$E), formula, data,
                                     REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights,  showWarnings=showWarnings, control=control, quiet=quiet,
                                     BPPARAM=BPPARAM, ...)
          }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod("fitExtractVarPartModel", "ExpressionSet",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitExtractVarPartModel( as.matrix(exprs(exprObj)), formula, data,
                                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                                     BPPARAM=BPPARAM,...)
          }
)

# sparseMatrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,sparseMatrix-method
setMethod("fitExtractVarPartModel", "sparseMatrix",
          function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = vpcontrol, quiet=FALSE, BPPARAM=SerialParam(), ...)
          {
            .fitExtractVarPartModel( exprObj, formula, data,
                                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                                     BPPARAM=BPPARAM, ...)
          }
)



