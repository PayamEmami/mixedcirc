search_before <- search()
warnings_seen <- character()
messages_seen <- character()

withCallingHandlers(
  library(mixedcirc),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    messages_seen <<- c(messages_seen, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)

new_search_entries <- setdiff(search(), search_before)
exports <- getNamespaceExports("mixedcirc")
missing_exports <- exports[!vapply(exports, function(name) {
  !inherits(
    try(getExportedValue("mixedcirc", name), silent = TRUE),
    "try-error"
  )
}, logical(1))]

stopifnot(
  length(warnings_seen) == 0L,
  length(messages_seen) == 0L,
  identical(new_search_entries, "package:mixedcirc"),
  length(missing_exports) == 0L
)
