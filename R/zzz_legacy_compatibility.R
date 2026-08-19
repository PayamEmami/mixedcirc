# Compatibility surface retained from the original package -----------------

# These helpers were never exported, but retaining their names prevents code
# written against mixedcirc::: internals from failing solely because the local
# variance engine and genomic annotation code were modernized.

add_sym <- function(x, include.backtick = "as.needed", dat = NULL) {
  include.backtick <- match.arg(include.backtick, c("as.needed", "all"))
  if (!is.character(x)) stop("x must be a character vector", call. = FALSE)
  quote <- if (include.backtick == "all") {
    rep(TRUE, length(x))
  } else {
    x != make.names(x)
  }
  if (!is.null(dat)) {
    if (!is.data.frame(dat)) stop("dat must be a data.frame", call. = FALSE)
    quote <- quote & x %in% names(dat)
  }
  x[quote] <- paste0("`", gsub("`", "\\`", x[quote], fixed = TRUE), "`")
  x
}

shallow <- function(x, cols = names(x), reset_class = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("shallow() requires the optional package data.table", call. = FALSE)
  }
  if (!data.table::is.data.table(x)) {
    stop("x must be a data.table", call. = FALSE)
  }
  if (!all(cols %in% names(x))) {
    stop("all cols must exist in x", call. = FALSE)
  }
  if (!is.logical(reset_class) || length(reset_class) != 1L || is.na(reset_class)) {
    stop("reset_class must be TRUE or FALSE", call. = FALSE)
  }
  result <- vector("list", length(cols))
  names(result) <- data.table::copy(cols)
  for (column in cols) result[[column]] <- x[[column]]
  data.table::setDT(result)
  if (reset_class) {
    data.table::setattr(result, "class", c("data.table", "data.frame"))
  } else {
    data.table::setattr(result, "class", class(x))
  }
  result
}

.isMixedModelFormula <- function(formula) {
  .mc_is_mixed_formula(stats::as.formula(formula))
}

.isDisconnected <- function() {
  # The original implementation tested a registered foreach cluster. The
  # replacement engine uses BiocParallel calls directly and does not retain a
  # cluster connection between calls, so there is no persistent connection to
  # become stale.
  FALSE
}

hasMissingData <- function(form, info) {
  variables <- all.vars(stats::as.formula(form))
  missing <- setdiff(variables, names(info))
  if (length(missing)) {
    stop("Variable in formula not found in data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  stats::setNames(
    vapply(variables, function(variable) anyNA(info[[variable]]), logical(1)),
    variables
  )
}

run_model_check_mixed <- function(fit,
                                  showWarnings = TRUE,
                                  dream = FALSE,
                                  colinearityCutoff = 0.999,
                                  immediate = FALSE) {
  if (!inherits(fit, "merMod")) {
    stop("fit must inherit from merMod", call. = FALSE)
  }
  .mc_check_model(
    fit,
    showWarnings = isTRUE(showWarnings) || isTRUE(dream),
    colinearityCutoff = colinearityCutoff,
    checkIntercept = !isTRUE(dream)
  )
  invisible(fit)
}

icount2 <- function(count) {
  if (missing(count)) {
    count <- NULL
  } else if (!is.numeric(count) || length(count) != 1L ||
             !is.finite(count) || count < 0 || count %% 1 != 0) {
    stop("count must be a non-negative integer", call. = FALSE)
  }
  index <- 0L
  next_element <- function() {
    if (is.null(index)) {
      NULL
    } else if (is.null(count) || index < count) {
      index <<- index + 1L
      index
    } else {
      index <<- NULL
      NULL
    }
  }
  structure(list(nextElem = next_element), class = c("abstractiter", "iter"))
}

exprIter <- function(exprObj,
                     weights,
                     useWeights = TRUE,
                     scale = TRUE,
                     iterCount = "icount") {
  expression <- as.matrix(exprObj)
  if (!is.null(weights) && !identical(dim(expression), dim(weights))) {
    stop("exprObj and weights must have identical dimensions", call. = FALSE)
  }
  if (!iterCount %in% c("icount", "icount2")) {
    stop("iterCount must be 'icount' or 'icount2'", call. = FALSE)
  }
  current <- 0L
  next_element <- function() {
    if (current >= nrow(expression)) {
      if (identical(iterCount, "icount2")) return(NULL)
      stop(structure(
        list(message = "StopIteration", call = NULL),
        class = c("StopIteration", "error", "condition")
      ))
    }
    current <<- current + 1L
    index <- current
    feature_weights <- if (isTRUE(useWeights) && !is.null(weights)) {
      as.numeric(weights[index, ])
    } else {
      NULL
    }
    if (isTRUE(scale) && !is.null(feature_weights)) {
      feature_weights <- feature_weights / mean(feature_weights)
    }
    list(
      E = expression[index, ],
      weights = feature_weights,
      n_iter = index,
      max_iter = nrow(expression)
    )
  }
  structure(list(nextElem = next_element), class = c("abstractiter", "iter"))
}

iterBatch <- function(exprObj,
                      weights,
                      useWeights = TRUE,
                      scale = TRUE,
                      n_chunks = nrow(exprObj) / 500,
                      min_chunk_size = 20,
                      BPPARAM = NULL) {
  expression <- as.matrix(exprObj)
  if (!nrow(expression)) stop("exprObj must contain at least one row", call. = FALSE)
  if (!is.null(weights) && !identical(dim(expression), dim(weights))) {
    stop("exprObj and weights must have identical dimensions", call. = FALSE)
  }
  if (!is.numeric(n_chunks) || length(n_chunks) != 1L ||
      !is.finite(n_chunks) || n_chunks <= 0) {
    stop("n_chunks must be positive", call. = FALSE)
  }
  if (!is.numeric(min_chunk_size) || length(min_chunk_size) != 1L ||
      !is.finite(min_chunk_size) || min_chunk_size < 1) {
    stop("min_chunk_size must be positive", call. = FALSE)
  }
  if (!is.null(BPPARAM) && methods::is(BPPARAM, "BiocParallelParam")) {
    workers <- max(1L, as.integer(BiocParallel::bpworkers(BPPARAM)))
    n_chunks <- ceiling(n_chunks / workers) * workers
  }
  maximum_chunks <- max(1L, floor(nrow(expression) / min_chunk_size))
  n_chunks <- max(1L, min(ceiling(n_chunks), maximum_chunks, nrow(expression)))
  indices <- parallel::splitIndices(nrow(expression), n_chunks)
  current <- 0L
  iterator <- function() {
    if (current == length(indices)) return(NULL)
    current <<- current + 1L
    selected <- indices[[current]]
    chunk <- expression[selected, , drop = FALSE]
    chunk_weights <- if (isTRUE(useWeights) && !is.null(weights)) {
      as.matrix(weights[selected, , drop = FALSE])
    } else {
      matrix(1, nrow(chunk), ncol(chunk))
    }
    if (isTRUE(scale) && isTRUE(useWeights) && !is.null(weights)) {
      chunk_weights <- chunk_weights / rowMeans(chunk_weights)
    }
    list(E = chunk, weights = chunk_weights)
  }
  attr(iterator, "n_chunks") <- length(indices)
  iterator
}

.fitVarPartModel <- function(exprObj,
                             formula,
                             data,
                             REML = FALSE,
                             useWeights = TRUE,
                             weightsMatrix = NULL,
                             showWarnings = TRUE,
                             fxn = identity,
                             colinearityCutoff = 0.999,
                             control = NULL,
                             quiet = FALSE,
                             BPPARAM = BiocParallel::SerialParam(),
                             ...) {
  global_weights <- .mc_capture_global_weights()
  on.exit(.mc_restore_global_weights(global_weights), add = TRUE)
  .mc_fit_models(
    exprObj = exprObj,
    formula = formula,
    data = data,
    REML = REML,
    useWeights = useWeights,
    weightsMatrix = weightsMatrix,
    showWarnings = showWarnings,
    fxn = fxn,
    control = control,
    quiet = quiet,
    BPPARAM = BPPARAM,
    colinearityCutoff = colinearityCutoff,
    allow_response_na = FALSE,
    dots = list(...)
  )
}

.fitExtractVarPartModel <- function(exprObj,
                                    formula,
                                    data,
                                    REML = FALSE,
                                    useWeights = TRUE,
                                    weightsMatrix = NULL,
                                    showWarnings = TRUE,
                                    colinearityCutoff = 0.999,
                                    control = NULL,
                                    quiet = FALSE,
                                    BPPARAM = BiocParallel::SerialParam(),
                                    ...) {
  models <- .fitVarPartModel(
    exprObj = exprObj,
    formula = formula,
    data = data,
    REML = REML,
    useWeights = useWeights,
    weightsMatrix = weightsMatrix,
    showWarnings = showWarnings,
    fxn = identity,
    colinearityCutoff = colinearityCutoff,
    control = control,
    quiet = quiet,
    BPPARAM = BPPARAM,
    ...
  )
  extractVarPart(models, showWarnings = FALSE)
}

.mc_tabix_functions <- function(functions) {
  .mc_methylkit_internals(functions)
}

.mc_bind_tabix_results <- function(results, data_table = FALSE) {
  results <- Filter(function(result) !is.null(result) && nrow(result) > 0L, results)
  if (!length(results)) {
    if (data_table) return(data.table::data.table())
    return(data.frame())
  }
  combined <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  if (data_table) combined else as.data.frame(combined)
}

applyTbxByOverlap_list <- function(tbxFile,
                                   ranges,
                                   chunk.size = 1e6,
                                   dir,
                                   filename,
                                   return.type = c("tabix", "data.frame", "data.table"),
                                   FUN,
                                   ...,
                                   tabixHead = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("This tabix helper requires data.table and Rsamtools", call. = FALSE)
  }
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  if (methods::is(ranges, "GRanges")) ranges <- list(ranges)
  if (!is.list(ranges) || !length(ranges)) {
    stop("ranges must be a non-empty list of GRanges objects", call. = FALSE)
  }
  functions <- .mc_tabix_functions(c(
    "getTabixByOverlap", ".write.table.noSci", "catsub2tabix"
  ))
  result_kind <- if (return.type == "data.table") "data.table" else "data.frame"
  results <- lapply(ranges, function(region) {
    value <- tryCatch(
      functions$getTabixByOverlap(
        tbxFile,
        granges = region,
        return.type = result_kind
      ),
      error = function(error) NULL
    )
    if (is.null(value) || !nrow(value)) return(NULL)
    FUN(value, ...)
  })
  if (return.type == "data.frame") {
    return(.mc_bind_tabix_results(results, data_table = FALSE))
  }
  if (return.type == "data.table") {
    return(.mc_bind_tabix_results(results, data_table = TRUE))
  }

  if (!dir.exists(dir)) stop("dir does not exist", call. = FALSE)
  nonempty <- Filter(function(result) !is.null(result) && nrow(result) > 0L, results)
  if (!length(nonempty)) stop("No tabix records overlapped ranges", call. = FALSE)
  prefix <- paste0(
    paste(sample(c(0:9, letters, LETTERS), 9L, replace = TRUE), collapse = ""),
    filename
  )
  for (index in seq_along(nonempty)) {
    functions$.write.table.noSci(
      nonempty[[index]],
      file.path(path.expand(dir), paste(index, prefix, sep = "_")),
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
  }
  path <- functions$catsub2tabix(
    dir,
    pattern = prefix,
    filename = filename,
    sort = TRUE,
    tabixHead = tabixHead
  )
  sub("[.]tbi$", "", path)
}

applyTbxByChunk_multiple <- function(
    tbxFile1,
    tbxFile2,
    chunk.size = 1e6,
    dir,
    filename,
    return.type = c("tabix", "data.frame", "data.table", "text"),
    FUN,
    ...,
    tabixHead = NULL,
    textHeader = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("This tabix helper requires data.table and Rsamtools", call. = FALSE)
  }
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  if (!is.numeric(chunk.size) || length(chunk.size) != 1L ||
      !is.finite(chunk.size) || chunk.size < 1 || chunk.size %% 1 != 0) {
    stop("chunk.size must be a positive integer", call. = FALSE)
  }
  functions <- .mc_tabix_functions(c(
    "getTabixByChunk", ".write.table.noSci", "catsub2tabix"
  ))
  file1 <- if (inherits(tbxFile1, "TabixFile")) tbxFile1 else {
    Rsamtools::TabixFile(tbxFile1, yieldSize = as.integer(chunk.size))
  }
  file2 <- if (inherits(tbxFile2, "TabixFile")) tbxFile2 else {
    Rsamtools::TabixFile(tbxFile2, yieldSize = as.integer(chunk.size))
  }
  if (Rsamtools::isOpen(file1)) close(file1)
  if (Rsamtools::isOpen(file2)) close(file2)
  Rsamtools::yieldSize(file1) <- as.integer(chunk.size)
  Rsamtools::yieldSize(file2) <- as.integer(chunk.size)
  records1 <- unname(Rsamtools::countTabix(file1)[[1L]])
  records2 <- unname(Rsamtools::countTabix(file2)[[1L]])
  if (!identical(as.numeric(records1), as.numeric(records2))) {
    stop("tabix files must contain the same number of records", call. = FALSE)
  }
  chunk_count <- ceiling(records1 / chunk.size)
  open(file1)
  open(file2)
  on.exit({
    if (Rsamtools::isOpen(file1)) close(file1)
    if (Rsamtools::isOpen(file2)) close(file2)
  }, add = TRUE)

  results <- lapply(seq_len(chunk_count), function(index) {
    first <- functions$getTabixByChunk(
      file1,
      chunk.size = NULL,
      return.type = if (return.type == "data.table") "data.table" else "data.frame"
    )
    second <- functions$getTabixByChunk(
      file2,
      chunk.size = NULL,
      return.type = if (return.type == "data.table") "data.table" else "data.frame"
    )
    FUN(first, second, ...)
  })
  close(file1)
  close(file2)

  if (return.type == "data.frame") {
    return(.mc_bind_tabix_results(results, data_table = FALSE))
  }
  if (return.type == "data.table") {
    return(.mc_bind_tabix_results(results, data_table = TRUE))
  }
  if (!dir.exists(dir)) stop("dir does not exist", call. = FALSE)
  prefix <- paste0(
    paste(sample(c(0:9, letters, LETTERS), 9L, replace = TRUE), collapse = ""),
    filename
  )
  for (index in seq_along(results)) {
    functions$.write.table.noSci(
      results[[index]],
      file.path(path.expand(dir), paste(index, prefix, sep = "_")),
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
  }
  if (return.type == "tabix") {
    path <- functions$catsub2tabix(
      dir,
      pattern = prefix,
      filename = filename,
      sort = TRUE,
      tabixHead = tabixHead
    )
    return(sub("[.]tbi$", "", path))
  }

  output <- file.path(path.expand(dir), filename)
  if (file.exists(output)) file.remove(output)
  if (!is.null(textHeader)) {
    connection <- file(output, open = "w", blocking = TRUE)
    tryCatch(
      write(textHeader, file = connection, ncolumns = length(textHeader), sep = "\t"),
      finally = close(connection)
    )
  } else {
    file.create(output)
  }
  temporary <- list.files(dir, pattern = prefix, full.names = TRUE)
  temporary <- temporary[order(as.integer(sub("_.*", "", basename(temporary))))]
  for (source in temporary) file.append(output, source)
  unlink(temporary)
  output
}

# Restore the original S4 generic surface while dispatching to the local,
# tested implementation above. An ANY method covers the original matrix,
# data.frame, EList, ExpressionSet, and sparseMatrix inputs.
methods::setGeneric(
  "checkModelStatus",
  signature = "fit",
  function(fit,
           showWarnings = TRUE,
           dream = FALSE,
           colinearityCutoff = 0.999,
           immediate = FALSE) {
    methods::standardGeneric("checkModelStatus")
  }
)

methods::setMethod(
  "checkModelStatus",
  "ANY",
  function(fit,
           showWarnings = TRUE,
           dream = FALSE,
           colinearityCutoff = 0.999,
           immediate = FALSE) {
    if (inherits(fit, "merMod")) {
      run_model_check_mixed(
        fit, showWarnings, dream, colinearityCutoff, immediate
      )
    } else {
      .mc_check_model(
        fit,
        showWarnings = isTRUE(showWarnings) || isTRUE(dream),
        colinearityCutoff = colinearityCutoff
      )
    }
    invisible(fit)
  }
)

methods::setGeneric(
  "fitVarPartModel",
  signature = "exprObj",
  function(exprObj,
           formula,
           data,
           REML = FALSE,
           useWeights = TRUE,
           weightsMatrix = NULL,
           showWarnings = TRUE,
           fxn = identity,
           control = NULL,
           quiet = FALSE,
           BPPARAM = BiocParallel::SerialParam(),
           ...) {
    methods::standardGeneric("fitVarPartModel")
  }
)

methods::setMethod(
  "fitVarPartModel",
  "ANY",
  function(exprObj,
           formula,
           data,
           REML = FALSE,
           useWeights = TRUE,
           weightsMatrix = NULL,
           showWarnings = TRUE,
           fxn = identity,
           control = NULL,
           quiet = FALSE,
           BPPARAM = BiocParallel::SerialParam(),
           ...) {
    .fitVarPartModel(
      exprObj, formula, data, REML, useWeights, weightsMatrix,
      showWarnings, fxn, 0.999, control, quiet, BPPARAM, ...
    )
  }
)

methods::setGeneric(
  "fitExtractVarPartModel",
  signature = "exprObj",
  function(exprObj,
           formula,
           data,
           REML = FALSE,
           useWeights = TRUE,
           weightsMatrix = NULL,
           showWarnings = TRUE,
           control = NULL,
           quiet = FALSE,
           BPPARAM = BiocParallel::SerialParam(),
           ...) {
    methods::standardGeneric("fitExtractVarPartModel")
  }
)

methods::setMethod(
  "fitExtractVarPartModel",
  "ANY",
  function(exprObj,
           formula,
           data,
           REML = FALSE,
           useWeights = TRUE,
           weightsMatrix = NULL,
           showWarnings = TRUE,
           control = NULL,
           quiet = FALSE,
           BPPARAM = BiocParallel::SerialParam(),
           ...) {
    .fitExtractVarPartModel(
      exprObj, formula, data, REML, useWeights, weightsMatrix,
      showWarnings, 0.999, control, quiet, BPPARAM, ...
    )
  }
)
