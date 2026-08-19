library(mixedcirc)

partial <- list(
  transcript = matrix(
    1:6,
    nrow = 3,
    dimnames = list(paste0("sample", 1:3), c("comp1", "comp2"))
  ),
  protein = matrix(
    7:12,
    nrow = 3,
    dimnames = list(paste0("sample", 1:3), c("comp1", "comp2"))
  )
)
loadings <- list(
  transcript = matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("comp1", "comp2"))
  ),
  protein = matrix(
    5:10,
    nrow = 3,
    dimnames = list(c("prot1", "prot2", "prot3"), c("comp1", "comp2"))
  )
)

integration_result <- methods::new(
  "mixedcirc_integration",
  partial = partial,
  average = matrix(1:6, nrow = 3),
  loadings = loadings,
  model = list(ok = TRUE)
)

merged_partial <- mixedcirc_getscore(
  integration_result,
  type = "partial",
  merge = TRUE
)
merged_loadings <- mixedcirc_getscore(
  integration_result,
  type = "loadings",
  merge = TRUE
)

stopifnot(
  identical(dim(merged_partial), c(3L, 4L)),
  identical(
    colnames(merged_partial),
    c(
      "transcript_comp1", "transcript_comp2",
      "protein_comp1", "protein_comp2"
    )
  ),
  identical(dim(merged_loadings), c(5L, 2L)),
  identical(
    rownames(merged_loadings),
    c(
      "transcript_gene1", "transcript_gene2",
      "protein_prot1", "protein_prot2", "protein_prot3"
    )
  ),
  identical(
    mixedcirc_getscore(integration_result, "loadings", "protein"),
    loadings$protein
  )
)

bad_partial <- integration_result
bad_partial@partial$protein <- bad_partial@partial$protein[1:2, , drop = FALSE]
partial_error <- try(
  mixedcirc_getscore(bad_partial, type = "partial", merge = TRUE),
  silent = TRUE
)
stopifnot(
  inherits(partial_error, "try-error"),
  grepl("same number of samples", partial_error, fixed = TRUE)
)

bad_loadings <- integration_result
bad_loadings@loadings$protein <- cbind(
  bad_loadings@loadings$protein,
  comp3 = 11:13
)
loadings_error <- try(
  mixedcirc_getscore(bad_loadings, type = "loadings", merge = TRUE),
  silent = TRUE
)
stopifnot(
  inherits(loadings_error, "try-error"),
  grepl("same number of components", loadings_error, fixed = TRUE)
)
