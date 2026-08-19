library(mixedcirc)

if (!requireNamespace("methylKit", quietly = TRUE)) {
  message("Skipping chunked RRBS test because methylKit is unavailable")
} else {
  set.seed(2101)
  locus_count <- 36L
  sample_count <- 8L
  sample_ids <- paste0("sample_", seq_len(sample_count))
  sample_times <- rep(c(0, 6, 12, 18), length.out = sample_count)

  coverage <- matrix(
    sample(45:180, locus_count * sample_count, replace = TRUE),
    nrow = locus_count,
    ncol = sample_count
  )
  baseline <- stats::plogis(stats::rnorm(locus_count, -0.2, 0.7))
  probability <- outer(
    baseline,
    0.7 * cos(2 * pi * sample_times / 24),
    function(base, rhythm) stats::plogis(stats::qlogis(base) + rhythm)
  )
  methylated <- matrix(
    stats::rbinom(
      locus_count * sample_count,
      size = as.vector(coverage),
      prob = as.vector(probability)
    ),
    nrow = locus_count,
    ncol = sample_count
  )
  unmethylated <- coverage - methylated

  methyl_data <- data.frame(
    chr = rep("chr1", locus_count),
    start = seq.int(100L, by = 10L, length.out = locus_count),
    end = seq.int(100L, by = 10L, length.out = locus_count),
    strand = rep("+", locus_count)
  )
  for (sample_index in seq_len(sample_count)) {
    methyl_data[[paste0("coverage", sample_index)]] <- coverage[, sample_index]
    methyl_data[[paste0("numCs", sample_index)]] <- methylated[, sample_index]
    methyl_data[[paste0("numTs", sample_index)]] <- unmethylated[, sample_index]
  }

  coverage_indices <- seq.int(5L, by = 3L, length.out = sample_count)
  methylated_indices <- coverage_indices + 1L
  unmethylated_indices <- coverage_indices + 2L
  make_database <- function(suffix) {
    directory <- tempfile(paste0("mixedcirc_rrbs_", suffix, "_"))
    dir.create(directory)
    methylKit:::makeMethylBaseDB(
      df = methyl_data,
      dbpath = directory,
      dbtype = "tabix",
      sample.ids = sample_ids,
      assembly = "test",
      context = "CpG",
      resolution = "base",
      treatment = rep(c(0, 1), each = sample_count / 2L),
      coverage.index = coverage_indices,
      numCs.index = methylated_indices,
      numTs.index = unmethylated_indices,
      destranded = FALSE,
      suffix = paste0("_", suffix)
    )
  }

  design <- data.frame(
    replicate_id = factor(rep(sample_ids, each = 2L)),
    scaler = rep(c(1, 0), sample_count),
    time = rep(sample_times, each = 2L)
  )
  design$inphase <- cos(2 * pi * design$time / 24)
  design$outphase <- sin(2 * pi * design$time / 24)

  mixed_database <- make_database("mixed")
  mixed_result <- suppressWarnings(mixedcirc_rrbs_voom_mixed(
    mixed_database,
    ~ scaler + scaler:inphase + scaler:outphase + (1 | replicate_id),
    design,
    chunk.size = 7,
    save.plot = TRUE,
    quiet = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  ))

  stopifnot(
    methods::is(mixed_result$E, "methylBaseDB"),
    methods::is(mixed_result$weights, "methylBaseDB"),
    identical(mixed_result$E@num.records, mixed_database@num.records),
    identical(mixed_result$weights@num.records, mixed_database@num.records),
    identical(names(mixed_result), c("E", "weights", "voom.xy", "voom.line")),
    length(mixed_result$voom.xy$x) == locus_count,
    length(mixed_result$voom.line$x) >= 2L
  )

  original <- methylKit::getData(mixed_database)
  transformed <- methylKit::getData(mixed_result$E)
  weighted <- methylKit::getData(mixed_result$weights)
  value_indices <- c(rbind(methylated_indices, unmethylated_indices))
  library_sizes <- rep(colSums(coverage), each = 2L)
  expected_log_cpm <- t(
    log2(t(as.matrix(original[, value_indices]) + 0.5) /
           (library_sizes + 1) * 1e6)
  )
  estimated_weights <- as.matrix(weighted[, value_indices])

  stopifnot(
    isTRUE(all.equal(
      unname(as.matrix(transformed[, value_indices])),
      unname(expected_log_cpm),
      tolerance = 1e-8
    )),
    all(is.finite(estimated_weights)),
    all(estimated_weights > 0),
    stats::sd(as.vector(estimated_weights)) > 0
  )

  # Fixed-effect formulas should return the same complete on-disk structure as
  # mixed-effect formulas, and sample-level library sizes should be expanded to
  # the paired observations.
  fixed_database <- make_database("fixed")
  fixed_result <- suppressWarnings(mixedcirc_rrbs_voom_mixed(
    fixed_database,
    ~ scaler + scaler:inphase + scaler:outphase,
    design,
    lib.size = colSums(coverage),
    chunk.size = 11,
    quiet = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  ))
  fixed_weights <- as.matrix(
    methylKit::getData(fixed_result$weights)[, value_indices]
  )
  stopifnot(
    identical(names(fixed_result), c("E", "weights")),
    all(is.finite(fixed_weights)),
    all(fixed_weights > 0),
    stats::sd(as.vector(fixed_weights)) > 0
  )

  invalid_library_error <- tryCatch(
    mixedcirc_rrbs_voom_mixed(
      fixed_database,
      ~ scaler,
      design,
      lib.size = 1:3
    ),
    error = conditionMessage
  )
  stopifnot(grepl("one value per sample", invalid_library_error, fixed = TRUE))
}
