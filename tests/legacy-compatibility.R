suppressPackageStartupMessages(library(mixedcirc))

expect_legacy_error <- function(expression, pattern = NULL) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) {
    stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  }
  invisible(error)
}

# Formula-name quoting works without modifying the caller's data object.
input_data <- data.frame(
  ordinary = 1:2,
  check.names = FALSE,
  "batch name" = 3:4
)
input_names <- names(input_data)
stopifnot(
  identical(
    mixedcirc:::add_sym(c("ordinary", "batch name", "absent"), dat = input_data),
    c("ordinary", "`batch name`", "absent")
  ),
  identical(
    mixedcirc:::add_sym(c("ordinary", "batch name"), include.backtick = "all"),
    c("`ordinary`", "`batch name`")
  ),
  identical(names(input_data), input_names)
)

# shallow() retains selected column pointers and optionally a custom subclass.
source_table <- data.table::data.table(a = 1:3, b = 4:6)
data.table::setattr(
  source_table,
  "class",
  c("mixedcirc_test_table", class(source_table))
)
shallow_table <- mixedcirc:::shallow(source_table, "a")
reset_table <- mixedcirc:::shallow(source_table, reset_class = TRUE)
stopifnot(
  data.table::is.data.table(shallow_table),
  identical(names(shallow_table), "a"),
  inherits(shallow_table, "mixedcirc_test_table"),
  identical(
    data.table::address(source_table[["a"]]),
    data.table::address(shallow_table[["a"]])
  ),
  identical(class(reset_table), c("data.table", "data.frame"))
)

# The formula, connection, and missing-data helpers retain their old contract.
metadata_missing <- data.frame(
  subject = factor(c("a", "a", "b", "b")),
  x = c(1, NA, 3, 4)
)
stopifnot(
  isTRUE(mixedcirc:::.isMixedModelFormula(~ x + (1 | subject))),
  isFALSE(mixedcirc:::.isMixedModelFormula(~ x)),
  identical(mixedcirc:::.isDisconnected(), FALSE),
  identical(
    mixedcirc:::hasMissingData(~ x + subject, metadata_missing),
    c(x = TRUE, subject = FALSE)
  )
)
expect_legacy_error(
  mixedcirc:::hasMissingData(~ absent, metadata_missing),
  "Variable in formula not found in data"
)

# icount2 returns NULL repeatedly after exhaustion, while the ordinary
# expression iterator signals StopIteration as required by iterators.
counter <- mixedcirc:::icount2(2)
stopifnot(
  identical(iterators::nextElem(counter), 1L),
  identical(iterators::nextElem(counter), 2L),
  is.null(iterators::nextElem(counter)),
  is.null(iterators::nextElem(counter))
)
expression <- matrix(
  1:12,
  nrow = 3,
  dimnames = list(paste0("feature_", 1:3), paste0("sample_", 1:4))
)
weights <- matrix(seq(0.5, 2, length.out = 12), nrow = 3)
iterator <- mixedcirc:::exprIter(expression, weights, iterCount = "icount2")
iterated <- lapply(seq_len(nrow(expression)), function(index) {
  iterators::nextElem(iterator)
})
stopifnot(
  is.null(iterators::nextElem(iterator)),
  identical(vapply(iterated, `[[`, integer(1), "n_iter"), 1:3),
  max(abs(vapply(iterated, function(value) mean(value$weights), numeric(1)) - 1)) <
    1e-12,
  identical(
    unname(do.call(rbind, lapply(iterated, `[[`, "E"))),
    unname(expression)
  )
)
stopping_iterator <- mixedcirc:::exprIter(
  expression[1, , drop = FALSE],
  weights[1, , drop = FALSE],
  iterCount = "icount"
)
invisible(iterators::nextElem(stopping_iterator))
stop_iteration <- expect_legacy_error(iterators::nextElem(stopping_iterator))
stopifnot(inherits(stop_iteration, "StopIteration"))

# Batch iteration covers every feature exactly once and scales weights by row.
batch_iterator <- mixedcirc:::iterBatch(
  expression,
  weights,
  n_chunks = 2,
  min_chunk_size = 1,
  BPPARAM = BiocParallel::SerialParam()
)
batches <- vector("list", attr(batch_iterator, "n_chunks"))
for (index in seq_along(batches)) batches[[index]] <- batch_iterator()
stopifnot(
  is.null(batch_iterator()),
  identical(do.call(rbind, lapply(batches, `[[`, "E")), expression),
  max(abs(rowMeans(do.call(rbind, lapply(batches, `[[`, "weights"))) - 1)) <
    1e-12
)

# Restore and exercise the original S4 generic surface and internal drivers.
stopifnot(
  methods::isGeneric("checkModelStatus"),
  methods::isGeneric("fitVarPartModel"),
  methods::isGeneric("fitExtractVarPartModel"),
  methods::is(
    methods::getMethod(
      "checkModelStatus", "ANY", where = asNamespace("mixedcirc")
    ),
    "MethodDefinition"
  ),
  methods::is(
    methods::getMethod(
      "fitVarPartModel", "ANY", where = asNamespace("mixedcirc")
    ),
    "MethodDefinition"
  ),
  methods::is(
    methods::getMethod(
      "fitExtractVarPartModel", "ANY", where = asNamespace("mixedcirc")
    ),
    "MethodDefinition"
  )
)
set.seed(410)
model_data <- data.frame(
  subject = factor(rep(seq_len(6), each = 2)),
  x = stats::rnorm(12)
)
subject_effect <- stats::rnorm(6)[model_data$subject]
model_expression <- rbind(
  feature_a = 2 + model_data$x + subject_effect + stats::rnorm(12, sd = 0.3),
  feature_b = 1 - 0.5 * model_data$x + stats::rnorm(12, sd = 0.5)
)
serial <- BiocParallel::SerialParam()
internal_models <- mixedcirc:::.fitVarPartModel(
  model_expression,
  ~ x,
  model_data,
  quiet = TRUE,
  BPPARAM = serial
)
generic_models <- fitVarPartModel(
  model_expression,
  ~ x,
  model_data,
  quiet = TRUE,
  BPPARAM = serial
)
internal_fractions <- mixedcirc:::.fitExtractVarPartModel(
  model_expression,
  ~ x,
  model_data,
  quiet = TRUE,
  BPPARAM = serial
)
generic_fractions <- fitExtractVarPartModel(
  model_expression,
  ~ x,
  model_data,
  quiet = TRUE,
  BPPARAM = serial
)
stopifnot(
  methods::is(internal_models, "VarParFitList"),
  methods::is(generic_models, "VarParFitList"),
  identical(internal_models@method, "lm"),
  length(internal_models) == nrow(model_expression),
  all(vapply(internal_models, inherits, logical(1), "lm")),
  max(abs(as.matrix(internal_fractions) - as.matrix(generic_fractions))) < 1e-12
)

ordinary_fit <- stats::lm(model_expression[1, ] ~ model_data$x)
stopifnot(identical(mixedcirc:::checkModelStatus(ordinary_fit), ordinary_fit))
expect_legacy_error(
  mixedcirc:::checkModelStatus(
    stats::lm(model_expression[1, ] ~ 0 + model_data$x)
  ),
  "No intercept term"
)
mixed_fit <- lme4::lmer(
  model_expression[1, ] ~ model_data$x + (1 | model_data$subject),
  control = lme4::lmerControl(
    calc.derivs = FALSE,
    check.conv.singular = list(action = "ignore", tol = 1e-4)
  )
)
stopifnot(
  identical(mixedcirc:::run_model_check_mixed(mixed_fit), mixed_fit),
  identical(mixedcirc:::checkModelStatus(mixed_fit), mixed_fit)
)
