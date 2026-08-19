library(mixedcirc)

annotation_packages <- c("GenomicRanges", "IRanges", "S4Vectors", "rtracklayer")
if (!all(vapply(annotation_packages, requireNamespace, logical(1), quietly = TRUE))) {
  message("Skipping genomic annotation tests because optional packages are unavailable")
} else {
  annotation <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(
      start = c(100, 200, 500, 600),
      end = c(149, 249, 549, 649)
    ),
    strand = c("+", "+", "-", "-")
  )
  S4Vectors::mcols(annotation) <- S4Vectors::DataFrame(
    type = rep("exon", 4),
    gene_id = rep(c("gene_A", "gene_B"), each = 2),
    transcript_id = rep(c("tx_A", "tx_B"), each = 2),
    exon_number = rep(1:2, 2)
  )

  introns <- construct_introns(annotation, update = FALSE)
  combined <- construct_introns(annotation, update = TRUE)
  stopifnot(
    methods::is(introns, "GRanges"),
    length(introns) == 2L,
    identical(IRanges::start(introns), c(150L, 550L)),
    identical(IRanges::end(introns), c(199L, 599L)),
    identical(as.character(S4Vectors::mcols(introns)$type), rep("intron", 2)),
    identical(as.character(S4Vectors::mcols(introns)$transcript_id), c("tx_A", "tx_B")),
    all(is.na(S4Vectors::mcols(introns)$exon_number)),
    length(combined) == 6L
  )

  gff_path <- tempfile(fileext = ".gff3")
  rtracklayer::export(annotation, gff_path, format = "gff3")
  features <- mixed_circ_tss_calculation(
    gff_path,
    gff = TRUE,
    id = "gene_id",
    up.flank = 50,
    down.flank = 20,
    unique.prom = FALSE
  )
  stopifnot(
    methods::is(features, "GRangesList"),
    identical(names(features), c("exons", "introns", "promoters", "TSSes")),
    length(features$exons) == 4L,
    length(features$introns) == 2L,
    identical(IRanges::start(features$TSSes), c(100L, 649L)),
    identical(IRanges::start(features$promoters), c(50L, 629L)),
    identical(IRanges::end(features$promoters), c(120L, 699L)),
    identical(as.character(S4Vectors::mcols(features$promoters)$name), c("gene_A", "gene_B"))
  )

  if (requireNamespace("genomation", quietly = TRUE)) {
    bed_path <- system.file(
      "extdata", "chr21.refseq.hg19.bed",
      package = "genomation"
    )
    if (nzchar(bed_path)) {
      bed_features <- suppressMessages(mixed_circ_tss_calculation(bed_path))
      stopifnot(
        methods::is(bed_features, "GRangesList"),
        identical(names(bed_features), c("exons", "introns", "promoters", "TSSes"))
      )
    }
  }
}
