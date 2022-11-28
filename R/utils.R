
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




# applyTbxByChunk
#' Serially apply a function on chunks of two tabix files
#'
#' The function reads the same chunks of two tabix files and applies a function on them.
#' The function (FUN argument) should apply on data.frames and
#' return a data frame
#' as a result. The function is serially applied to chunks
#' (means no parallelization).
#' However, the function FUN itself can be a parallelized function
#' and related arguments could be passed on to the function via ... argument.
#'
#' @param tbxFile1 tabix file to read. a TabixFile object
#' @param tbxFile2 tabix file to read. a TabixFile object
#' @param chunk.size number of rows to be taken as a chunk, default: 1e6
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes two data.frames and returns a
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN.
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be
#' a path, just a file name.
#' @param tabixHead optional header to prepend to the file
#'
#' @return either a path to a tabix or text file, or a data frame or data.table
#' @noRd

applyTbxByChunk_multiple <-
  function (tbxFile1,tbxFile2, chunk.size = 1e+06, dir, filename, return.type = c("tabix",
                                                                                  "data.frame", "data.table", "text"), FUN, ..., tabixHead = NULL,
            textHeader = NULL)
  {
    return.type <- match.arg(return.type)
    FUN <- match.fun(FUN)
    if (class(tbxFile1) != "TabixFile" | class(tbxFile2) != "TabixFile") {
      tbxFile1 <- Rsamtools::TabixFile(tbxFile1, yieldSize = chunk.size)
      tbxFile2 <- Rsamtools::TabixFile(tbxFile2, yieldSize = chunk.size)
    }
    else {
      if (Rsamtools::isOpen(tbxFile1)) {
        close(tbxFile1)
      }

      Rsamtools::yieldSize(tbxFile1) <- chunk.size

      if (Rsamtools::isOpen(tbxFile2)) {
        close(tbxFile2)
      }
      Rsamtools::yieldSize(tbxFile2) <- chunk.size
    }
    recs = Rsamtools::countTabix(tbxFile1)[[1]]
    chunk.num = ceiling(recs/chunk.size)
    open(tbxFile1)
    open(tbxFile2)
    if (return.type == "tabix") {
      myFunc <- function(chunk.num, tbxFile1,tbxFile2, dir, filename,
                         FUN, ...) {
        data1 = methylKit:::getTabixByChunk(tbxFile1, chunk.size = NULL,
                                            return.type = "data.frame")
        data2 = methylKit:::getTabixByChunk(tbxFile2, chunk.size = NULL,
                                            return.type = "data.frame")
        res = FUN(data1,data2, ...)
        outfile = file.path(path.expand(dir), paste(chunk.num,
                                                    filename, sep = "_"))
        methylKit:::.write.table.noSci(res, outfile, quote = FALSE, col.names = FALSE,
                                       row.names = FALSE, sep = "\t")
      }
      rndFile = paste(sample(c(0:9, letters, LETTERS), 9, replace = TRUE),
                      collapse = "")
      filename2 = paste(rndFile, filename, sep = "_")
      res = lapply(1:chunk.num, myFunc, tbxFile1,tbxFile2, dir, filename2,
                   FUN, ...)
      path <- methylKit:::catsub2tabix(dir, pattern = filename2, filename,
                                       sort = TRUE, tabixHead = tabixHead)
      return(gsub(".tbi", "", path))
    }
    else if (return.type == "text") {
      myFunc2 <- function(chunk.num, tbxFile1,tbxFile2, dir, filename,
                          FUN, ...) {
        data1 = methylKit:::getTabixByChunk(tbxFile1, chunk.size = NULL,
                                            return.type = "data.frame")
        data2 = methylKit:::getTabixByChunk(tbxFile2, chunk.size = NULL,
                                            return.type = "data.frame")
        res = FUN(data1,data2, ...)
        outfile = file.path(path.expand(dir), paste(chunk.num,
                                                    filename, sep = "_"))
        methylKit:::.write.table.noSci(res, outfile, quote = FALSE, col.names = FALSE,
                                       row.names = FALSE, sep = "\t")
      }
      rndFile = paste(sample(c(0:9, letters, LETTERS), 9, replace = TRUE),
                      collapse = "")
      filename2 = paste(rndFile, filename, sep = "_")
      res = lapply(1:chunk.num, myFunc2, tbxFile1,tbxFile2, dir, filename2,
                   FUN, ...)
      outfile = file.path(path.expand(dir), filename)
      if (file.exists(outfile)) {
        message("overwriting ", outfile)
        unlink(outfile)
      }
      con = file(outfile, open = "a", blocking = TRUE)
      if (!is.null(textHeader))
        write(file = con, x = textHeader, ncolumns = length(textHeader),
              sep = "\t")
      for (file in gtools::mixedsort(list.files(path = dir,
                                                pattern = filename2, full.names = TRUE))) {
        file.append(outfile, file)
      }
      close(con)
      unlink(list.files(path = dir, pattern = filename2, full.names = TRUE))
      return(outfile)
    }
    else if (return.type == "data.frame") {
      myFunc3 <- function(chunk.num, tbxFile1,tbxFile2, FUN, ...) {
        data1 = methylKit:::getTabixByChunk(tbxFile1, chunk.size = NULL,
                                            return.type = "data.frame")
        data2 = methylKit:::getTabixByChunk(tbxFile2, chunk.size = NULL,
                                            return.type = "data.frame")
        FUN(data1,data2, ...)
      }
      res = lapply(1:chunk.num, myFunc3, tbxFile1,tbxFile2, FUN, ...)
      data.frame(data.table::rbindlist(res))
    }
    else {
      myFunc4 <- function(chunk.num, tbxFile1,tbxFile2, FUN, ...) {
        data1 = methylKit:::getTabixByChunk(tbxFile1, chunk.size = NULL,
                                            return.type = "data.table")
        data2 = methylKit:::getTabixByChunk(tbxFile2, chunk.size = NULL,
                                            return.type = "data.table")
        FUN(data1,data2, ...)
      }
      res = lapply(1:chunk.num, myFunc4, tbxFile1,tbxFile2, FUN, ...)
      data.table::rbindlist(res)
    }
  }
