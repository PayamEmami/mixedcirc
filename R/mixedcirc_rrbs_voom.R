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
#' @param ... other arguments are passed to lmer
#' @examples

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
#' By default NAs are imputed by zero.
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



mixedcirc_rrbs_voom<-function (counts, formula, data, lib.size = NULL,chunk.size=100,
                               span = 0.5, plot = FALSE, save.plot = FALSE, quiet = FALSE,
                               BPPARAM = bpparam(), ...)
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
                                                    x[is.na(x)]<-0


                                                    as.data.frame(t(colSums(x[,counts@coverage.index])))[,,drop=F]


                                                  })


    lib.size<-colSums(lib.size_partial)

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
                                   xx[is.na(xx)]<-0
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


  if (variancePartition:::.isMixedModelFormula(formula)) {
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
                                     xx<-x[,-c(1:4,counts@coverage.index),drop=F];
                                     Ameans<-rowMeans(xx,na.rm = TRUE)[,drop=F]


                                     xx_res<-x;
                                     vpList = variancePartition:::fitVarPartModel(xx, formula, data, showWarnings = FALSE,
                                                                                  fxn = function(fit) {
                                                                                    list(sd = attr(lme4::VarCorr(fit), "sc"), fitted.values = predict(fit))
                                                                                  }, BPPARAM = BPPARAM)
                                     fitted.values <- lapply(vpList, function(x) x$fitted.values)
                                     fitted.values <- do.call("rbind", fitted.values)
                                     fit = list()
                                     fit$sigma <- sapply(vpList, function(x) x$sd)
                                     fit$df.residual = rep(2, length(fit$sigma))
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
                                                                  y_norm@numTs.index[!y_norm@sample.ids%in%c("sigma"  ,   "resid"    , "Ameans" )]),drop=F]))


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
