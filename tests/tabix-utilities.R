suppressPackageStartupMessages(library(mixedcirc))

tabix_packages <- c(
  "data.table", "GenomicRanges", "IRanges", "methylKit", "Rsamtools"
)
if (!all(vapply(tabix_packages, requireNamespace, logical(1), quietly = TRUE))) {
  message("Skipping Tabix utility tests because optional packages are unavailable")
} else {
  locus_count <- 7L
  methyl_data <- data.frame(
    chr = rep("chr1", locus_count),
    start = seq.int(100L, by = 10L, length.out = locus_count),
    end = seq.int(100L, by = 10L, length.out = locus_count),
    strand = rep("+", locus_count),
    coverage1 = 20:26,
    numCs1 = 5:11,
    numTs1 = rep(15L, locus_count),
    coverage2 = 30:36,
    numCs2 = 10:16,
    numTs2 = rep(20L, locus_count)
  )
  make_database <- function(tag) {
    database_directory <- tempfile(paste0("mixedcirc_tabix_", tag, "_"))
    dir.create(database_directory)
    methylKit:::makeMethylBaseDB(
      df = methyl_data,
      dbpath = database_directory,
      dbtype = "tabix",
      sample.ids = c("sample_1", "sample_2"),
      assembly = "test",
      context = "CpG",
      resolution = "base",
      treatment = c(0, 1),
      coverage.index = c(5L, 8L),
      numCs.index = c(6L, 9L),
      numTs.index = c(7L, 10L),
      destranded = FALSE,
      suffix = paste0("_", tag)
    )
  }
  first_database <- make_database("first")
  second_database <- make_database("second")
  first_path <- methylKit::getDBPath(first_database)
  second_path <- methylKit::getDBPath(second_database)
  output_directory <- tempfile("mixedcirc_tabix_outputs_")
  dir.create(output_directory)

  add_chunks <- function(first, second) {
    first[, 5:10] <- first[, 5:10] + second[, 5:10]
    first
  }
  check_added_data <- function(value) {
    stopifnot(
      nrow(value) == locus_count,
      identical(as.integer(value[[5L]]), 2L * methyl_data$coverage1),
      identical(as.integer(value[[10L]]), 2L * methyl_data$numTs2)
    )
  }

  chunk_frame <- mixedcirc:::applyTbxByChunk_multiple(
    first_path,
    second_path,
    chunk.size = 3L,
    dir = output_directory,
    filename = "chunk_frame",
    return.type = "data.frame",
    FUN = add_chunks
  )
  chunk_table <- mixedcirc:::applyTbxByChunk_multiple(
    Rsamtools::TabixFile(first_path),
    Rsamtools::TabixFile(second_path),
    chunk.size = 3L,
    dir = output_directory,
    filename = "chunk_table",
    return.type = "data.table",
    FUN = add_chunks
  )
  check_added_data(chunk_frame)
  check_added_data(chunk_table)
  stopifnot(
    is.data.frame(chunk_frame),
    !data.table::is.data.table(chunk_frame),
    data.table::is.data.table(chunk_table)
  )

  text_path <- mixedcirc:::applyTbxByChunk_multiple(
    first_path,
    second_path,
    chunk.size = 3L,
    dir = output_directory,
    filename = "chunk_output.txt",
    return.type = "text",
    FUN = add_chunks,
    textHeader = names(methyl_data)
  )
  text_lines <- readLines(text_path)
  stopifnot(
    file.exists(text_path),
    length(text_lines) == locus_count + 1L,
    identical(strsplit(text_lines[[1L]], "\t", fixed = TRUE)[[1L]], names(methyl_data))
  )

  internals <- mixedcirc:::.mc_methylkit_internals()
  tabix_header <- internals$.formatTabixHeader(
    class = "methylBaseDB",
    tabixHead = internals$makeTabixHeader(
      mixedcirc:::.mc_methylbase_slots(first_database)
    )
  )
  chunk_tabix <- mixedcirc:::applyTbxByChunk_multiple(
    first_path,
    second_path,
    chunk.size = 3L,
    dir = output_directory,
    filename = "chunk_tabix",
    return.type = "tabix",
    FUN = add_chunks,
    tabixHead = tabix_header
  )
  stopifnot(
    file.exists(chunk_tabix),
    file.exists(paste0(chunk_tabix, ".tbi")),
    unname(Rsamtools::countTabix(chunk_tabix)[[1L]]) == locus_count
  )

  ranges <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(95L, 145L), end = c(115L, 165L))
  )
  range_list <- list(ranges[1L], ranges[2L])
  overlap_frame <- mixedcirc:::applyTbxByOverlap_list(
    first_path,
    range_list,
    dir = output_directory,
    filename = "overlap_frame",
    return.type = "data.frame",
    FUN = identity
  )
  overlap_table <- mixedcirc:::applyTbxByOverlap_list(
    Rsamtools::TabixFile(first_path),
    ranges,
    dir = output_directory,
    filename = "overlap_table",
    return.type = "data.table",
    FUN = identity
  )
  stopifnot(
    is.data.frame(overlap_frame),
    !data.table::is.data.table(overlap_frame),
    data.table::is.data.table(overlap_table),
    nrow(overlap_frame) == 4L,
    nrow(overlap_table) == 4L,
    identical(as.integer(overlap_frame[[2L]]), c(100L, 110L, 150L, 160L))
  )

  overlap_tabix <- mixedcirc:::applyTbxByOverlap_list(
    first_path,
    range_list,
    dir = output_directory,
    filename = "overlap_tabix",
    return.type = "tabix",
    FUN = identity,
    tabixHead = tabix_header
  )
  stopifnot(
    file.exists(overlap_tabix),
    file.exists(paste0(overlap_tabix, ".tbi")),
    unname(Rsamtools::countTabix(overlap_tabix)[[1L]]) == 4L
  )

  # Record-count mismatches are detected from the files themselves.
  short_data <- methyl_data[-1L, , drop = FALSE]
  short_directory <- tempfile("mixedcirc_tabix_short_")
  dir.create(short_directory)
  short_database <- methylKit:::makeMethylBaseDB(
    df = short_data,
    dbpath = short_directory,
    dbtype = "tabix",
    sample.ids = c("sample_1", "sample_2"),
    assembly = "test",
    context = "CpG",
    resolution = "base",
    treatment = c(0, 1),
    coverage.index = c(5L, 8L),
    numCs.index = c(6L, 9L),
    numTs.index = c(7L, 10L),
    destranded = FALSE,
    suffix = "_short"
  )
  mismatch_error <- tryCatch(
    mixedcirc:::applyTbxByChunk_multiple(
      first_path,
      methylKit::getDBPath(short_database),
      chunk.size = 3L,
      dir = output_directory,
      filename = "mismatch",
      return.type = "data.frame",
      FUN = add_chunks
    ),
    error = identity
  )
  stopifnot(
    inherits(mismatch_error, "error"),
    grepl("same number of records", conditionMessage(mismatch_error), fixed = TRUE)
  )
}
