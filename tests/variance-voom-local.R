suppressPackageStartupMessages(library(mixedcirc))

expect_voom_error <- function(expression, pattern) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  invisible(error)
}

set.seed(920)
n_subjects <- 8L
n_replicates <- 3L
n_samples <- n_subjects * n_replicates
n_features <- 60L
metadata <- data.frame(
  id = factor(rep(seq_len(n_subjects), each = n_replicates)),
  x = rnorm(n_samples)
)
rownames(metadata) <- paste0("sample_", seq_len(n_samples))
subject_effect <- rnorm(n_subjects, sd = 0.5)[metadata$id]
means <- exp(
  outer(seq(-0.4, 0.4, length.out = n_features), metadata$x) +
    rep(subject_effect, each = n_features) + log(70)
)
counts <- matrix(
  stats::rnbinom(length(means), mu = as.vector(means), size = 10),
  nrow = n_features,
  dimnames = list(paste0("gene_", seq_len(n_features)), rownames(metadata))
)
serial <- BiocParallel::SerialParam()

# The fixed model has a directly checkable logCPM transformation and produces
# positive finite weights.
fixed <- voomWithDreamWeights(
  counts,
  ~ x,
  metadata,
  quiet = TRUE,
  save.plot = TRUE,
  BPPARAM = serial
)
library_sizes <- colSums(counts)
expected_log_cpm <- t(log2(t(counts + 0.5) / (library_sizes + 1) * 1e6))
stopifnot(
  methods::is(fixed, "EList"),
  identical(dim(fixed$E), dim(counts)),
  identical(dim(fixed$weights), dim(counts)),
  max(abs(fixed$E - expected_log_cpm)) < 1e-12,
  all(is.finite(fixed$weights)),
  all(fixed$weights > 0),
  !is.null(fixed$voom.xy),
  !is.null(fixed$voom.line)
)

# Mixed-model, automatic smoothing, and prior-weight routes are all local.
sample_weights <- seq(0.7, 1.3, length.out = n_samples)
mixed <- voomWithDreamWeights(
  counts,
  ~ x + (1 | id),
  metadata,
  weights = sample_weights,
  rescaleWeightsAfter = TRUE,
  span = "auto",
  quiet = TRUE,
  BPPARAM = serial
)
count_weights <- voomWithDreamWeights(
  counts,
  ~ x,
  metadata,
  priorWeightsAsCounts = TRUE,
  scaledByLib = TRUE,
  quiet = TRUE,
  BPPARAM = serial
)
stopifnot(
  all(is.finite(mixed$weights)),
  all(mixed$weights > 0),
  length(mixed$targets$sample.weights) == n_samples,
  all(is.finite(count_weights$weights)),
  all(count_weights$weights > 0)
)

# DGEList and ExpressionSet inputs retain their feature/sample annotations.
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
dge_result <- voomWithDreamWeights(
  dge,
  ~ x + (1 | id),
  metadata,
  quiet = TRUE,
  BPPARAM = serial
)
expression_set <- Biobase::ExpressionSet(
  assayData = counts,
  phenoData = Biobase::AnnotatedDataFrame(metadata)
)
expression_set_result <- voomWithDreamWeights(
  expression_set,
  ~ x,
  metadata,
  quiet = TRUE,
  BPPARAM = serial
)
stopifnot(
  identical(dim(dge_result$weights), dim(counts)),
  identical(dim(expression_set_result$weights), dim(counts)),
  nrow(dge_result$targets) == n_samples,
  nrow(expression_set_result$targets) == n_samples
)

# Missing metadata are dropped from every associated object.
metadata_missing <- metadata
metadata_missing$x[1] <- NA_real_
missing_result <- suppressWarnings(voomWithDreamWeights(
  counts,
  ~ x,
  metadata_missing,
  quiet = TRUE,
  BPPARAM = serial
))
stopifnot(ncol(missing_result$E) == n_samples - 1L, nrow(missing_result$targets) == n_samples - 1L)

# SOCK-parallel mixed-model weights are numerically reproducible.
parallel_backend <- BiocParallel::SnowParam(
  workers = 2,
  type = "SOCK",
  progressbar = FALSE,
  stop.on.error = TRUE
)
parallel_mixed <- voomWithDreamWeights(
  counts,
  ~ x + (1 | id),
  metadata,
  weights = sample_weights,
  rescaleWeightsAfter = TRUE,
  span = "auto",
  quiet = TRUE,
  BPPARAM = parallel_backend
)
BiocParallel::bpstop(parallel_backend)
stopifnot(
  max(abs(parallel_mixed$E - mixed$E)) < 1e-12,
  max(abs(parallel_mixed$weights - mixed$weights)) < 1e-4
)

# Input failures do not fall through into obscure model errors.
expect_voom_error(
  voomWithDreamWeights(counts[1, , drop = FALSE], ~ x, metadata, quiet = TRUE),
  "at least two features"
)
negative_counts <- counts
negative_counts[1, 1] <- -1
expect_voom_error(
  voomWithDreamWeights(negative_counts, ~ x, metadata, quiet = TRUE),
  "non-negative"
)
expect_voom_error(
  voomWithDreamWeights(counts, ~ absent, metadata, quiet = TRUE),
  "Variable in formula not found"
)
expect_voom_error(
  voomWithDreamWeights(counts, ~ x, metadata, weights = rep(1, 2), quiet = TRUE),
  "sample weights"
)
