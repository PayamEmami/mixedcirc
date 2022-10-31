
add_sym<-
  function (x, include.backtick = "as.needed", dat = NULL)
  {
    original.format.dt <- FALSE
    len.x <- length(x)
    if (include.backtick == "all") {
      w <- 1:len.x
    }
    if (include.backtick == "as.needed") {
      if (is.null(dat)) {
        w <- which(x != make.names(names = x))
      }
      if (!is.null(dat)) {
        if (is.data.frame(dat)) {
          original.format.dt <- is.data.table(x = dat)
        }
        data.table::setDT(dat)
        requires.backtick <- logical(length = len.x)
        for (i in 1:len.x) {
          value.exists <- !is.null(tryCatch(expr = unique(eval(parse(text = sprintf("dat[, `%s`]",
                                                                                    eval((eval(x[i]))))))), error = function(e) return(NULL)))
          if (value.exists == TRUE & x[i] != make.names(x[i])) {
            requires.backtick[i] <- TRUE
          }
        }
        w <- which(requires.backtick == TRUE)
      }
    }
    if (is.data.frame(x = dat)) {
      if (original.format.dt == F) {
        setDF(x = dat)
      }
    }
    if (length(w) > 0) {
      x[w] <- sprintf("`%s`", x[w])
    }
    return(x)
  }



#' Serially apply a function on the list regions of tabix files
#'
#' The function reads a list of regions of a tabix file and applies a function on them.
#' The function (FUN argument) should apply on data.frames and return a data frame
#' as a result. The function is serially applied to chunks (means no parallelization).
#' However, the function FUN itself can be a parallelized function
#' and related arguments could be passed on to the function via ... argument.
#'
#'
#' @param tbxFile tabix file to read. a TabixFile object
#' @param ranges a list of GRanges object specifying the regions
#' @param chunk.size number of rows to be taken as a chunk, default: 1e6
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes a data.frame and returns a
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN.
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be
#' a path, just a file name.
#' @param tabixHead optional header to prepend to the file
#'
#' @return either a path to a tabix or text file, or a data frame or data.table
#' @noRd
applyTbxByOverlap_list<-function(tbxFile,ranges,chunk.size=1e6,dir,filename,
                                 return.type=c("tabix","data.frame","data.table"),
                                 FUN,...,tabixHead=NULL){

  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)

  # open tabix file with given chunk size
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile)

  }
  #   else {
  #     if(Rsamtools::isOpen(tbxFile)){close(tbxFile)}# close if already open
  #     Rsamtools::yieldSize(tbxFile) <-  chunk.size
  #   }


  if(return.type =="tabix"){

    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chunk.num,region.split,tbxFile,dir,filename,FUN,...){
      data <- try(expr = data <- methylKit:::getTabixByOverlap(
        tbxFile,granges = region.split[[chunk.num]],
        return.type="data.frame"),silent = TRUE)

      if( class(data)== "try-error") {

        #         warning( paste("No records found in range between",
        # min(IRanges::end(region.split[[chunk.num]])),
        #                        "and",max(IRanges::end(region.split[[chunk.num]])),
        #                        "for Chromosome",
        # unique(as.character(region.split[[chunk.num]]@seqnames))))

      } else {

        res=FUN(data,...)

        # for tabix
        outfile= file.path(path.expand( dir),paste(chunk.num,filename,sep="_"))
        methylKit:::.write.table.noSci(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                                       sep="\t")
      }
    }

    # attach a random string to the file name
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")

    # apply function to chunks
    res=lapply(1:length(ranges),myFunc,ranges,tbxFile,dir,filename2,FUN,...)

    # collect & cat temp files,then make tabix
    path <- methylKit:::catsub2tabix(dir,pattern=filename2,filename,sort = TRUE,
                                     tabixHead = tabixHead)

    return(gsub(".tbi","",path))

  } else if(return.type=="data.frame"){

    # create a custom function that contains the function
    # to be applied
    myFunc2<-function(chunk.num,region.split,tbxFile,FUN,...){
      data <- try(expr = data <- methylKit:::getTabixByOverlap(
        tbxFile,granges = region.split[[chunk.num]],
        return.type="data.frame"),silent = TRUE)
      if( !(class(data)== "try-error") ) {
        res=FUN(data,...)
      }
    }

    res=lapply(1:length(ranges),myFunc2,ranges,tbxFile,FUN,...)

    # collect and return
    data.frame(rbindlist(res))
  }else{

    myFunc3<-function(chunk.num,region.split,tbxFile,FUN,...){
      data <- try(expr = data <- methylKit:::getTabixByOverlap(
        tbxFile,
        granges = region.split[[chunk.num]],
        return.type="data.table"),silent = TRUE)

      if( !(class(data)[1] == "try-error") ) { ## class of data.table is both "data.table" and "data.frame
        res=FUN(data,...)
      }
    }

    res=lapply(1:length(ranges),myFunc3,ranges,tbxFile,FUN,...)


    # collect and return
    rbindlist(res)
  }

}
