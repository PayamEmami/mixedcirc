library(mixedcirc)

model <- stats::lm(mpg ~ wt, data = mtcars)
fit_one <- methods::new(
  "mixedcirc_fit",
  results = data.frame(mesor = 1, amplitude = 2, phase = -0.5),
  fit = model,
  exp_design = data.frame(time = mtcars$wt, measure = mtcars$mpg),
  type = "expression",
  call = quote(mixedcirc_detect()),
  params = list(period = 24)
)
fit_two <- fit_one
fit_two@results$mesor <- 3

fit_list <- methods::new(
  "mixedcirc_fit_list",
  results = list(feature_a = fit_one, feature_b = fit_two)
)

subset_one <- fit_list["feature_b"]
subset_list <- fit_list[c("feature_a", "feature_b")]
extracted <- fit_list[["feature_a"]]
result_column <- fit_one[, "amplitude"]

print_one <- capture.output(returned_one <- print(fit_one))
print_list <- capture.output(returned_list <- print(fit_list))

stopifnot(
  length(fit_list) == 2L,
  methods::is(subset_one, "mixedcirc_fit"),
  identical(subset_one@results$mesor, 3),
  methods::is(subset_list, "mixedcirc_fit_list"),
  identical(names(subset_list@results), c("feature_a", "feature_b")),
  methods::is(extracted, "mixedcirc_fit"),
  is.data.frame(result_column),
  identical(names(result_column), "amplitude"),
  identical(returned_one, fit_one),
  identical(returned_list, fit_list),
  length(print_one) > 0L,
  any(grepl("Total number of variables: 2", print_list, fixed = TRUE))
)

missing_index <- try(fit_list[[]], silent = TRUE)
stopifnot(inherits(missing_index, "try-error"))
