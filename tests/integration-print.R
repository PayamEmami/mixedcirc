library(mixedcirc)

integration_result <- methods::new(
  "mixedcirc_integration",
  partial = list(
    transcriptome = matrix(1, nrow = 1, ncol = 1),
    metabolome = matrix(2, nrow = 1, ncol = 1)
  ),
  average = matrix(1.5, nrow = 1, ncol = 1),
  loadings = list(
    transcriptome = matrix(1, nrow = 1, ncol = 1),
    metabolome = matrix(1, nrow = 1, ncol = 1)
  ),
  model = NULL
)

print_output <- capture.output(print_result <- print(integration_result))
show_output <- capture.output(methods::show(integration_result))

stopifnot(
  identical(print_result, integration_result),
  any(grepl("Experimental mixedcirc integration result", print_output, fixed = TRUE)),
  any(grepl("Omics blocks: 2", print_output, fixed = TRUE)),
  any(grepl("independent data", print_output, fixed = TRUE)),
  identical(print_output, show_output)
)
