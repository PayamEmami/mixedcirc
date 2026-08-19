suppressPackageStartupMessages(library(mixedcirc))

# The installed package must not declare, load, or require an external variance
# decomposition package.
description <- utils::packageDescription("mixedcirc")
external_package <- paste0("variance", "Partition")
stopifnot(
  !grepl(external_package, description$Imports, fixed = TRUE),
  !external_package %in% loadedNamespaces()
)

expect_error <- function(expression, pattern) {
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

set.seed(812)
n_id <- 8L
n_rep <- 4L
n <- n_id * n_rep
metadata <- data.frame(
  id = factor(rep(seq_len(n_id), each = n_rep)),
  x = rnorm(n),
  z = runif(n, -1, 1)
)
rownames(metadata) <- paste0("sample_", seq_len(n))

subject_effect <- rnorm(n_id, sd = 1.1)
expression <- rbind(
  gene_1 = 1.0 + 0.8 * metadata$x + subject_effect[metadata$id] + rnorm(n, sd = 0.5),
  gene_2 = 2.0 - 0.4 * metadata$x + 0.4 * subject_effect[metadata$id] + rnorm(n, sd = 0.8),
  gene_3 = 0.5 + 0.3 * metadata$z + 1.4 * subject_effect[metadata$id] + rnorm(n, sd = 0.4)
)
colnames(expression) <- rownames(metadata)
serial <- BiocParallel::SerialParam()

# Mixed models: compare every component with calculations from direct lme4
# fits, without using another variance-decomposition package.
mixed_formula <- ~ x + (1 | id)
observed_mixed <- fitExtractVarPartModel(
  expression,
  mixed_formula,
  metadata,
  BPPARAM = serial,
  quiet = TRUE
)
direct_mixed <- lapply(seq_len(nrow(expression)), function(index) {
  direct_data <- transform(metadata, response = expression[index, ])
  lme4::lmer(
    response ~ x + (1 | id),
    data = direct_data,
    REML = FALSE,
    control = lme4::lmerControl(
      calc.derivs = FALSE,
      check.rankX = "stop.deficient",
      check.conv.singular = list(action = "ignore", tol = 1e-4)
    )
  )
})
direct_components <- t(vapply(direct_mixed, function(fit) {
  random_variance <- as.numeric(lme4::VarCorr(fit)$id[1, 1])
  design <- lme4::getME(fit, "X")
  fixed_contribution <- design[, "x"] * lme4::fixef(fit)[["x"]]
  fixed_variance <- stats::var(fixed_contribution) *
    (nrow(design) - 1) / nrow(design)
  components <- c(
    id = random_variance,
    x = fixed_variance,
    Residuals = stats::sigma(fit)^2
  )
  components / sum(components)
}, numeric(3)))

stopifnot(
  methods::is(observed_mixed, "varPartResults"),
  identical(dim(observed_mixed), dim(direct_components)),
  all(is.finite(as.matrix(observed_mixed))),
  max(abs(rowSums(as.data.frame(observed_mixed)) - 1)) < 1e-10,
  max(abs(as.matrix(observed_mixed) - direct_components)) < 1e-7
)

# Fixed models: with one predictor the fitted fraction is ordinary R-squared.
single_fixed <- fitExtractVarPartModel(
  expression,
  ~ x,
  metadata,
  BPPARAM = serial,
  quiet = TRUE
)
direct_r_squared <- vapply(seq_len(nrow(expression)), function(index) {
  summary(stats::lm(expression[index, ] ~ metadata$x))$r.squared
}, numeric(1))
stopifnot(
  max(abs(single_fixed$x - direct_r_squared)) < 1e-12,
  max(abs(single_fixed$Residuals - (1 - direct_r_squared))) < 1e-12
)

# Multiple fixed terms must be finite, named, and exhaustive.
fixed <- fitExtractVarPartModel(
  expression,
  ~ x + z,
  metadata,
  BPPARAM = serial,
  quiet = TRUE
)
stopifnot(
  identical(colnames(fixed), c("x", "z", "Residuals")),
  all(is.finite(as.matrix(fixed))),
  max(abs(rowSums(as.data.frame(fixed)) - 1)) < 1e-10
)

# The model-returning API, custom extraction, REML, and weights all work.
weights <- matrix(
  runif(length(expression), 0.5, 2),
  nrow = nrow(expression),
  dimnames = dimnames(expression)
)
weights_before <- weights
weighted <- fitExtractVarPartModel(
  expression,
  mixed_formula,
  metadata,
  weightsMatrix = weights,
  BPPARAM = serial,
  quiet = TRUE
)
models <- fitVarPartModel(
  expression,
  mixed_formula,
  metadata,
  weightsMatrix = weights,
  BPPARAM = serial,
  quiet = TRUE
)
reml_sigma <- fitVarPartModel(
  expression,
  mixed_formula,
  metadata,
  REML = TRUE,
  fxn = stats::sigma,
  BPPARAM = serial,
  quiet = TRUE
)
extracted <- mixedcirc:::extractVarPart(models, showWarnings = FALSE)
stopifnot(
  identical(weights, weights_before),
  methods::is(models, "VarParFitList"),
  identical(models@method, "lmer"),
  length(models) == nrow(expression),
  all(vapply(models, inherits, logical(1), "lmerMod")),
  all(is.finite(unlist(reml_sigma))),
  max(abs(as.matrix(extracted) - as.matrix(weighted))) < 1e-12,
  max(abs(rowSums(as.data.frame(weighted)) - 1)) < 1e-10
)

# Every documented input container is accepted.
data_frame_result <- fitExtractVarPartModel(
  as.data.frame(expression), ~ x, metadata, BPPARAM = serial, quiet = TRUE
)
sparse_result <- fitExtractVarPartModel(
  Matrix::Matrix(expression, sparse = TRUE),
  ~ x,
  metadata,
  BPPARAM = serial,
  quiet = TRUE
)
expression_set <- Biobase::ExpressionSet(
  assayData = expression,
  phenoData = Biobase::AnnotatedDataFrame(metadata)
)
expression_set_result <- fitExtractVarPartModel(
  expression_set, ~ x, metadata, BPPARAM = serial, quiet = TRUE
)
elist <- methods::new("EList", list(E = expression, weights = weights))
elist_result <- fitExtractVarPartModel(
  elist, mixed_formula, metadata, BPPARAM = serial, quiet = TRUE
)
stopifnot(
  max(abs(as.matrix(data_frame_result) - as.matrix(single_fixed))) < 1e-12,
  max(abs(as.matrix(sparse_result) - as.matrix(single_fixed))) < 1e-12,
  max(abs(as.matrix(expression_set_result) - as.matrix(single_fixed))) < 1e-12,
  max(abs(as.matrix(elist_result) - as.matrix(weighted))) < 1e-12
)

# Missing metadata are removed consistently from expression and weights.
metadata_missing <- metadata
metadata_missing$x[2] <- NA_real_
missing_result <- suppressWarnings(fitExtractVarPartModel(
  expression,
  mixed_formula,
  metadata_missing,
  weightsMatrix = weights,
  BPPARAM = serial,
  quiet = TRUE,
  showWarnings = FALSE
))
manual_result <- fitExtractVarPartModel(
  expression[, -2, drop = FALSE],
  mixed_formula,
  droplevels(metadata[-2, , drop = FALSE]),
  weightsMatrix = weights[, -2, drop = FALSE],
  BPPARAM = serial,
  quiet = TRUE,
  showWarnings = FALSE
)
stopifnot(max(abs(as.matrix(missing_result) - as.matrix(manual_result))) < 1e-12)

# The internal fitting helper supports feature-wise missing responses.
expression_na <- expression
expression_na[1, 3] <- NA_real_
na_models <- suppressWarnings(mixedcirc:::variance_fit(
  expression_na,
  ~ x,
  metadata,
  ignore_na = TRUE,
  quiet = TRUE,
  BPPARAM = serial
))
stopifnot(length(na_models) == nrow(expression), all(vapply(na_models, inherits, logical(1), "lm")))

# SOCK parallel results must agree with serial execution.
parallel_backend <- BiocParallel::SnowParam(
  workers = 2,
  type = "SOCK",
  progressbar = FALSE,
  stop.on.error = TRUE
)
parallel_result <- fitExtractVarPartModel(
  expression,
  mixed_formula,
  metadata,
  BPPARAM = parallel_backend,
  quiet = TRUE
)
BiocParallel::bpstop(parallel_backend)
stopifnot(max(abs(as.matrix(parallel_result) - as.matrix(observed_mixed))) < 1e-7)

# Failures are explicit and feature-specific.
expect_error(
  fitVarPartModel(expression, ~ missing_value, metadata, quiet = TRUE),
  "Variable in formula not found"
)
expect_error(
  fitVarPartModel(expression, ~ x, metadata[-1, ], quiet = TRUE),
  "number of samples"
)
bad_weights <- weights
bad_weights[1, 1] <- 0
expect_error(
  fitVarPartModel(expression, ~ x, metadata, weightsMatrix = bad_weights, quiet = TRUE),
  "finite positive"
)
constant_expression <- expression
constant_expression[1, ] <- 1
expect_error(
  fitVarPartModel(constant_expression, ~ x, metadata, quiet = TRUE),
  "zero or undefined variance"
)
