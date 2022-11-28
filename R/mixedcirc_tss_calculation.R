#' @title shallow copy a \code{data.table}
#'
#' @description Convenience function to shallow copy a \code{data.table}
#' (until this function is exported in the \code{data.table} package). For
#' internal use only.
#' @param x A \code{data.table}.
#' @param cols Character vector of column names (from \code{x}).
#' @param reset_class logical (default \code{FALSE}). If \code{TRUE}, resets
#' the class to \code{"data.table", "data.frame"}.
#' @return A shallow copied \code{data.table}.
#' @examples
#' \dontrun{
#' # For internal use only
#' library(data.table)
#' x <- data.table(a=1:2, b=3:4)
#' setattr(x, 'class', c("tmp", class(x)))
#'
#' y <- gread:::shallow(x) # only copies column pointers
#' class(y) # class(x) is retained
#' }
shallow <- function(x, cols = names(x), reset_class = FALSE) {
  stopifnot(is.data.table(x), all(cols %in% names(x)))
  ans = vector("list", length(cols))
  setattr(ans, 'names', data.table::copy(cols))
  for (col in cols)
    ans[[col]] = x[[col]]
  setDT(ans)
  class = if (!reset_class) data.table::copy(class(x))
  else c("data.table", "data.frame")
  setattr(ans, 'class', class)
  ans[]
}


#' @title Construct introns from GRanges object
#'
#' @description This function generates intronic coordinates by extracting
#' all the \code{exons} from a \code{GRanges} object.
#'
#' @param x An object of class \code{GRanges}. It has to have \code{exon}
#' feature present.
#' @param update If \code{TRUE} (default), \code{x} is updated by reference
#' and returned invisibily, else just the \code{intron} coordinates are
#' returned.
#' @return An object of class \code{GRanges} with updated \code{intron}
#' coordinates or just the intron coordinates depending on \code{update}.
#' They all inherit from \code{GRanges}.
#' @examples
#' gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
#' gff = gffToGRanges(gff.file)
#' construct_introns(gff)
#' @import data.table
construct_introns <- function(x, update=TRUE) {
  # stopifnot(is.ra(x) || is.gff(x))
  x = data.table::as.data.table(x)
  stopifnot("type" %chin% names(x),
            "transcript_id" %chin% names(x),
            update %in% c(TRUE, FALSE))

  # to please R CMD CHECK
  type=seqnames=transcript_id=NULL
  x_class = copy(class(x))
  exons = x[type == "exon"]
  if (!nrow(exons)) stop("type == 'exon' returned 0 rows.")
  setorderv(exons, c("transcript_id", "seqnames", "start", "end", "strand"))
  introns = exons[, .(seqnames = seqnames[1L],
                      start = end[seq_len(.N-1L)]+1L,
                      end = start[seq_len(.N-1L)+1L]-1L,
                      type = "intron",
                      strand = strand[1L]), by=transcript_id]
  introns = na.omit(introns, cols = c("start", "end"))
  check = introns[, any(start > end), by=transcript_id]
  if (length(ids <- which(check[["V1"]])))
    stop("Exons for transcript ids [", paste(ids, collaspse=" "),
         "] have start > end")
  exons[, c(setdiff(names(introns), "transcript_id")) := NULL]
  exons = unique(shallow(exons, reset_class=TRUE), by="transcript_id")
  ecols = names(exons)
  introns[exons, (ecols) := mget(ecols), on="transcript_id"]
  # reset all exon related cols
  introns[, grep("^exon", names(introns), value=TRUE) := NA]
  colorder = names(x)
  if (update) {
    x = rbind(x, introns)
    setorderv(x, c("seqnames", "start", "end"))
    x = x[, .SD, by="transcript_id"]
    setattr(x, 'class', x_class)
    setcolorder(x, colorder)
    ans = x
  } else {
    setcolorder(introns, colorder)
    setattr(introns, 'class', c( "data.table", "data.frame"))
    ans = introns
  }
  return(as(setDF(ans), "GRanges"))
}

# ---------------------------------------------------------------------------- #
#' Function for reading exon intron and promoter structure from a given bed or GFF file
#'
#' @param location location of the bed file with 12 or more columns.
#'                 The file can end in \code{.gz}, \code{.bz2}, \code{.xz}, or \code{.zip}
#'                 and/or start with \code{http://} or \code{ftp://}. If the file is not compressed
#'                 it can also start with \code{https://} or \code{ftps://}. or location of GFF file
#' @param gff If TRUE, location is assumed to be GFF file
#' @param id This is supposed to be a column of GFF file showing the symbol or id of each genes
#' @param remove.unusual remove the chromomesomes with unsual names, mainly random chromsomes etc (only in bed file)
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom get only the unique promoters, promoter boundaries will not have
#'                    a gene name if you set this option to be TRUE
#' @usage readTranscriptFeatures(location,remove.unusual=TRUE,
#'                               up.flank=1000,down.flank=1000,unique.prom=TRUE)
#' @return a \code{\link{GRangesList}} containing locations of exon/intron/promoter/TSS
#' @note  one bed track per file is only accepted, the bed files with multiple tracks will cause en error.
#'
#' @details This function is an extension of readTranscriptFeatures function from genomation package that now supports GFF file.
#' The idea behind GFF features is that they are calculated on gene levels in contrast to transcript levels. Therefore the parameter id must be supplied
#' together with GFF file so that the function knows what genes to use. For each unique genes, the function collapse the coordinate to contain all the genomics features.
#' These new coordinates are used to calculate TSS and promoter regions.
#'
#'
#' @examples
#' my.bed12.file = system.file("extdata/chr21.refseq.hg19.bed", package = "genomation")
#' my.bed12.file
#' feats = mixed_circ_tss_calculation(my.bed12.file)
#' names(feats)
#' sapply(feats, head)
#'
#' gff.file = system.file('extdata/chr21.refseq.hg19.gtf', package='genomation')
#' feats = mixed_circ_tss_calculation(gff.file,gff=TRUE)
#' names(feats)
#'
#' @import data.table
#' @import GenomicRanges
#' @import genomation
mixed_circ_tss_calculation<-function(location,gff=FALSE,id="gene_id",remove.unusual=TRUE,
                       up.flank=1000,
                       down.flank=1000,
                       unique.prom=TRUE){

  if(!gff)
  {
    # readBed6
    message('Reading the table...\r')
    bed=genomation:::readTableFast(location,header=FALSE,skip="auto")
    if(remove.unusual)
      bed=bed[grep("_", as.character(bed[,1]),invert=TRUE),]

    # introns
    message('Calculating intron coordinates...\r')
    introns    = genomation:::convertBed2Introns(bed)
    # exons
    message('Calculating exon coordinates...\r')
    exons    = genomation:::convertBed2Exons(bed)

    # get the locations of TSSes
    message('Calculating TSS coordinates...\r')
    tss=bed
    #  + strand
    tss[tss$V6=="+",3] = tss[tss$V6=="+",2]
    #  - strand
    tss[tss$V6=="-",2]=tss[tss$V6=="-",3]

    tssg = GenomicRanges:::GRanges(seqnames=as.character(tss$V1),
                   ranges=IRanges:::IRanges(start=tss$V2, end=tss$V3),
                   strand=as.character(tss$V6),
                   score=rep(0,nrow(tss)),
                   name=tss$V4)

    message('Calculating promoter coordinates...\r')
    # get the locations of promoters
    # + strand
    bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
    bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank

    #  - strand
    bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
    bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank


    if(! unique.prom){
      prom.df = (bed[,c(1,2,3,4,6)])
      prom = GenomicRanges:::GRanges(seqnames=as.character(prom.df$V1),
                     ranges = IRanges:::IRanges(start=prom.df$V2, end=prom.df$V3),
                     strand = as.character(prom.df$V6),
                     score=rep(0,nrow(prom.df)),
                     name=prom.df$V4)
    }else{
      prom.df = unique(bed[,c(1,2,3,6)])
      prom = GenomicRanges::GRanges(seqnames=as.character(prom.df$V1),
                     ranges=IRanges:::IRanges(start=prom.df$V2, end=prom.df$V3),
                     strand=as.character(prom.df$V6),
                     score=rep(0,nrow(prom.df)),
                     name=rep(".",nrow(prom.df)) )
    }

    message('Outputting the final GRangesList...\r\n')
    GenomicRanges:::GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
  }else{

    # read the gff file
    message('Reading the GFF file...\r')
    gff_range<-genomation::gffToGRanges(location)
    if(!id%in%colnames( GenomicRanges::mcols(gff_range)))
      stop(paste0("Did not find ",id," in ",location,"!"))

    message('Calculating intron coordinates...\r')
    introns<-construct_introns(gff_range,update = FALSE)
    introns<-as.data.frame(introns)
    introns = GenomicRanges::GRanges(seqnames=as.character(introns$seqnames),
                                   ranges=IRanges::IRanges(start=introns$start, end=introns$end),
                                   strand=as.character(introns$strand),
                                   score=rep(0,nrow(introns)),
                                   name=introns[, id])

    message('Calculating exon coordinates...\r')
    exons <- gff_range[gff_range$type=="exon"]

    exons<-exons[,c("score",id)]

    exons<-as.data.frame(exons)
    exons = GenomicRanges::GRanges(seqnames=as.character(exons$seqnames),
                                  ranges=IRanges::IRanges(start=exons$start, end=exons$end),
                                  strand=as.character(exons$strand),
                                  score=rep(0,nrow(exons)),
                                  name=exons[, id])

    message(paste0("Calculating TSS coordinates based on ",id,"...\r"))
    splited_ids<-rtracklayer::split(gff_range,all_ids<-GenomicRanges::mcols(gff_range)[,id])
    merged_coordniates <-unlist(range(splited_ids))

    tss<-data.table::as.data.table(merged_coordniates)
    #  + strand
    tss[tss$strand=="+","end"] = tss[tss$strand=="+","start"]
    #  - strand
    tss[tss$strand=="-","start"] = tss[tss$strand=="-","end"]

    tssg = GenomicRanges::GRanges(seqnames=as.character(tss$seqnames),
                   ranges=IRanges::IRanges(start=tss$start, end=tss$end),
                   strand=as.character(tss$strand),
                   score=rep(0,nrow(tss)),
                   name=names(merged_coordniates))

    message(paste0("Calculating promoter coordinates based on ",id,"...\r"))

    bed<-data.table::as.data.table(merged_coordniates)
    # get the locations of promoters
    # + strand
    bed[bed$strand=="+","end"]=bed[bed$strand=="+","start"]+down.flank
    bed[bed$strand=="+","start"]=bed[bed$strand=="+","start"]-up.flank

    #  - strand
    bed[bed$strand=="-","start"]=bed[bed$strand=="-","end"]-down.flank
    bed[bed$strand=="-","end"]=bed[bed$strand=="-","end"]+up.flank



    if(! unique.prom){
      prom.df = (bed[,c(1,2,3,4,5)])
      prom = GenomicRanges::GRanges(seqnames=as.character(prom.df$seqnames),
                     ranges = IRanges::IRanges(start=prom.df$start, end=prom.df$end),
                     strand = as.character(prom.df$strand),
                     score=rep(0,nrow(prom.df)),
                     name=names(merged_coordniates))
    }else{
      prom.df = unique(bed[,c(1,2,3,5)])
      prom = GenomicRanges::GRanges(seqnames=as.character(prom.df$seqnames),
                     ranges=IRanges::IRanges(start=prom.df$start, end=prom.df$end),
                     strand=as.character(prom.df$strand),
                     score=rep(0,nrow(prom.df)),
                     name=rep(".",nrow(prom.df)) )
    }

    GenomicRanges::GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
  }

}





