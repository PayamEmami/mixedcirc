library(mixedcirc)

factory <- getFromNamespace(".mc_boot_factory", "mixedcirc")
runner <- factory(
  function(i) {
    message("replicate message")
    warning("replicate warning")
    if (i == 2L) stop("replicate error")
    c(estimate = i)
  },
  errval = c(estimate = NA_real_)
)

successful <- runner(1L)
failed <- runner(2L)

stopifnot(
  identical(successful[[1L]], 1L),
  identical(attr(successful, "factory-message"), "replicate message\n"),
  identical(attr(successful, "factory-warning"), "replicate warning"),
  is.null(attr(successful, "factory-error")),
  is.na(failed[[1L]]),
  identical(attr(failed, "factory-error"), "replicate error")
)
