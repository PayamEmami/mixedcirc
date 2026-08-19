# Optional genomic annotation compatibility helpers ------------------------

.mc_require_annotation_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "This genomic annotation workflow requires: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Construct introns from exon annotations
#'
#' Constructs the gaps between consecutive exons belonging to the same
#' transcript. It supports the BED/GFF annotation workflow for assigning loci
#' to intronic regions without modifying objects by reference.
#'
#' @param x A `GenomicRanges::GRanges` object whose metadata contain `type` and
#'   `transcript_id` columns.
#' @param update Logical. If `TRUE`, append the constructed introns to `x`;
#'   otherwise return only the introns.
#'
#' @return A `GenomicRanges::GRanges` object.
#' @export
construct_introns <- function(x, update = TRUE) {
  .mc_require_annotation_packages(c("GenomicRanges", "IRanges", "S4Vectors"))
  if (!methods::is(x, "GRanges")) {
    stop("x must be a GRanges object", call. = FALSE)
  }
  if (!is.logical(update) || length(update) != 1L || is.na(update)) {
    stop("update must be TRUE or FALSE", call. = FALSE)
  }

  metadata <- as.data.frame(S4Vectors::mcols(x), stringsAsFactors = FALSE)
  required <- c("type", "transcript_id")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop(
      "x metadata must contain: ", paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  exon_indices <- which(as.character(metadata$type) == "exon")
  if (!length(exon_indices)) {
    stop("type == 'exon' returned 0 rows", call. = FALSE)
  }
  exons <- x[exon_indices]
  transcript_ids <- as.character(S4Vectors::mcols(exons)$transcript_id)
  if (anyNA(transcript_ids) || any(!nzchar(transcript_ids))) {
    stop("Every exon must have a non-missing transcript_id", call. = FALSE)
  }

  split_exons <- split(exons, transcript_ids)
  intron_ranges <- lapply(names(split_exons), function(transcript_id) {
    transcript <- split_exons[[transcript_id]]
    sequence_names <- unique(as.character(GenomicRanges::seqnames(transcript)))
    strands <- unique(as.character(GenomicRanges::strand(transcript)))
    if (length(sequence_names) != 1L || length(strands) != 1L) {
      stop(
        "Transcript ", transcript_id,
        " spans multiple sequences or strands",
        call. = FALSE
      )
    }
    transcript <- transcript[order(IRanges::start(transcript), IRanges::end(transcript))]
    if (length(transcript) < 2L) return(NULL)

    intron_start <- IRanges::end(transcript)[-length(transcript)] + 1L
    intron_end <- IRanges::start(transcript)[-1L] - 1L
    if (any(intron_start > intron_end)) {
      stop(
        "Overlapping or out-of-order exons in transcript ", transcript_id,
        call. = FALSE
      )
    }

    result <- GenomicRanges::GRanges(
      seqnames = rep(sequence_names, length(intron_start)),
      ranges = IRanges::IRanges(start = intron_start, end = intron_end),
      strand = rep(strands, length(intron_start))
    )
    template <- as.data.frame(
      S4Vectors::mcols(transcript[1L]),
      stringsAsFactors = FALSE
    )
    intron_metadata <- template[rep(1L, length(result)), , drop = FALSE]
    intron_metadata$type <- "intron"
    exon_columns <- grep("^exon", names(intron_metadata), value = TRUE)
    for (column in exon_columns) intron_metadata[[column]] <- NA
    S4Vectors::mcols(result) <- S4Vectors::DataFrame(intron_metadata)
    result
  })
  intron_ranges <- Filter(Negate(is.null), intron_ranges)
  if (length(intron_ranges)) {
    introns <- do.call(c, unname(intron_ranges))
  } else {
    introns <- x[0]
  }

  if (!update) return(introns)
  combined <- c(x, introns)
  combined[order(
    as.character(GenomicRanges::seqnames(combined)),
    IRanges::start(combined),
    IRanges::end(combined)
  )]
}

.mc_standard_annotation <- function(ranges, ids) {
  S4Vectors::mcols(ranges) <- S4Vectors::DataFrame(
    score = rep(0, length(ranges)),
    name = as.character(ids)
  )
  ranges
}

#' Read exon, intron, promoter, and transcription-start annotations
#'
#' Reads a BED12 or GFF/GTF annotation and returns the feature groups used by
#' the mixedcirc RRBS annotation workflow.
#'
#' @param location Path or supported URL for a BED12, GFF, or GTF file.
#' @param gff Logical; treat `location` as GFF/GTF rather than BED12.
#' @param id Metadata column identifying genes in a GFF/GTF file.
#' @param remove.unusual Logical; remove unusual chromosome names for BED input.
#' @param up.flank Bases upstream of each transcription start included in a
#'   promoter.
#' @param down.flank Bases downstream of each transcription start included in a
#'   promoter.
#' @param unique.prom Logical; collapse duplicate promoter coordinates.
#'
#' @return A `GenomicRanges::GRangesList` named `exons`, `introns`,
#'   `promoters`, and `TSSes`.
#' @export
mixed_circ_tss_calculation <- function(location,
                                       gff = FALSE,
                                       id = "gene_id",
                                       remove.unusual = TRUE,
                                       up.flank = 1000,
                                       down.flank = 1000,
                                       unique.prom = TRUE) {
  logical_arguments <- list(
    gff = gff,
    remove.unusual = remove.unusual,
    unique.prom = unique.prom
  )
  invalid_logical <- !vapply(
    logical_arguments,
    function(value) is.logical(value) && length(value) == 1L && !is.na(value),
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(names(logical_arguments)[which(invalid_logical)[1L]],
         " must be TRUE or FALSE", call. = FALSE)
  }
  flanks <- list(up.flank = up.flank, down.flank = down.flank)
  invalid_flank <- !vapply(
    flanks,
    function(value) is.numeric(value) && length(value) == 1L &&
      is.finite(value) && value >= 0 && value %% 1 == 0,
    logical(1)
  )
  if (any(invalid_flank)) {
    stop(names(flanks)[which(invalid_flank)[1L]],
         " must be one non-negative integer", call. = FALSE)
  }
  up.flank <- as.integer(up.flank)
  down.flank <- as.integer(down.flank)

  if (!gff) {
    .mc_require_annotation_packages("genomation")
    return(genomation::readTranscriptFeatures(
      location,
      remove.unusual = remove.unusual,
      up.flank = up.flank,
      down.flank = down.flank,
      unique.prom = unique.prom
    ))
  }

  .mc_require_annotation_packages(c(
    "GenomicRanges", "IRanges", "S4Vectors", "rtracklayer"
  ))
  annotation <- rtracklayer::import(location)
  metadata <- as.data.frame(S4Vectors::mcols(annotation), stringsAsFactors = FALSE)
  if (!id %in% names(metadata)) {
    stop("Did not find ", id, " in ", location, call. = FALSE)
  }
  if (!"type" %in% names(metadata)) {
    stop("The annotation must contain a type column", call. = FALSE)
  }
  gene_ids <- as.character(metadata[[id]])
  if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop("The selected id column contains missing gene identifiers", call. = FALSE)
  }

  exon_indices <- which(as.character(metadata$type) == "exon")
  if (!length(exon_indices)) {
    stop("The annotation contains no exon records", call. = FALSE)
  }
  exons <- .mc_standard_annotation(annotation[exon_indices], gene_ids[exon_indices])

  introns_raw <- construct_introns(annotation, update = FALSE)
  intron_ids <- S4Vectors::mcols(introns_raw)[[id]]
  introns <- .mc_standard_annotation(introns_raw, intron_ids)

  collapsed_ids <- unique(gene_ids)
  collapsed_ranges <- lapply(collapsed_ids, function(gene_id) {
    gene_annotation <- annotation[gene_ids == gene_id]
    sequence_names <- unique(as.character(GenomicRanges::seqnames(gene_annotation)))
    strands <- unique(as.character(GenomicRanges::strand(gene_annotation)))
    if (length(sequence_names) != 1L || length(strands) != 1L) {
      stop(
        "Gene ", gene_id, " spans multiple sequences or strands",
        call. = FALSE
      )
    }
    GenomicRanges::GRanges(
      seqnames = sequence_names,
      ranges = IRanges::IRanges(
        min(IRanges::start(gene_annotation)),
        max(IRanges::end(gene_annotation))
      ),
      strand = strands
    )
  })
  gene_ranges <- do.call(c, unname(collapsed_ranges))
  names(gene_ranges) <- collapsed_ids
  gene_strands <- as.character(GenomicRanges::strand(gene_ranges))
  starts <- IRanges::start(gene_ranges)
  ends <- IRanges::end(gene_ranges)
  tss_position <- ifelse(gene_strands == "-", ends, starts)

  transcription_starts <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gene_ranges),
    ranges = IRanges::IRanges(tss_position, tss_position),
    strand = GenomicRanges::strand(gene_ranges)
  )
  transcription_starts <- .mc_standard_annotation(
    transcription_starts, collapsed_ids
  )

  promoter_start <- ifelse(
    gene_strands == "-",
    ends - down.flank,
    starts - up.flank
  )
  promoter_end <- ifelse(
    gene_strands == "-",
    ends + up.flank,
    starts + down.flank
  )
  promoter_start <- pmax(1L, promoter_start)
  promoters <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gene_ranges),
    ranges = IRanges::IRanges(promoter_start, promoter_end),
    strand = GenomicRanges::strand(gene_ranges)
  )
  if (unique.prom) {
    promoters <- unique(promoters)
    promoter_ids <- rep(".", length(promoters))
  } else {
    promoter_ids <- collapsed_ids
  }
  promoters <- .mc_standard_annotation(promoters, promoter_ids)

  GenomicRanges::GRangesList(
    exons = exons,
    introns = introns,
    promoters = promoters,
    TSSes = transcription_starts
  )
}
