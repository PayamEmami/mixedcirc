library(mixedcirc)

set.seed(101)
test_data <- transform(mtcars, cyl_factor = factor(cyl))
test_data$mpg[c(2, 7)] <- NA_real_

weighted_fit <- lm(
  mpg ~ wt + factor(cyl),
  data = test_data,
  weights = seq_len(nrow(test_data)),
  na.action = na.exclude,
  x = TRUE,
  y = TRUE
)

simulated_response <- simulate(weighted_fit, nsim = 1)[[1L]]
refitted <- lme4::refit(weighted_fit, simulated_response)

stopifnot(
  inherits(refitted, "lm"),
  identical(names(coef(refitted)), names(coef(weighted_fit))),
  nrow(model.frame(refitted)) == stats::nobs(weighted_fit),
  length(refitted$weights) == stats::nobs(weighted_fit),
  all(is.finite(coef(refitted)))
)

bootstrap <- mixedcirc:::boot_lm(
  weighted_fit,
  coef,
  nsim = 5,
  seed = 102
)

stopifnot(
  identical(attr(bootstrap, "bootFail"), 0L),
  all(is.finite(bootstrap$t))
)

evaluation <- 0L
failed_bootstrap <- tryCatch(
  {
    mixedcirc:::boot_lm(
      weighted_fit,
      function(fit) {
        evaluation <<- evaluation + 1L
        if (evaluation > 1L) stop("deliberate test failure")
        coef(fit)
      },
      nsim = 2,
      seed = 103
    )
    NULL
  },
  error = conditionMessage
)

stopifnot(
  is.character(failed_bootstrap),
  grepl("All bootstrap runs failed (2/2)", failed_bootstrap, fixed = TRUE),
  grepl("deliberate test failure", failed_bootstrap, fixed = TRUE)
)
