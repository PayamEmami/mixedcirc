library(mixedcirc)

wrap <- getFromNamespace(".mc_circular_wrap", "mixedcirc")
boot_ci <- getFromNamespace(".mc_boot_ci", "mixedcirc")
boot_pvalues <- getFromNamespace(".mc_boot_pvalues", "mixedcirc")

set.seed(816)

# A phase near zero should have an ordinary narrow interval.
t0_zero <- c(phase_x = 0.1)
tmat_zero <- matrix(
  rnorm(2000, mean = 0.1, sd = 0.05),
  ncol = 1,
  dimnames = list(NULL, "phase_x")
)
ci_zero <- boot_ci(t0_zero, tmat_zero, type = "perc")
stopifnot(
  ci_zero$boot.lower < t0_zero,
  ci_zero$boot.upper > t0_zero,
  (ci_zero$boot.upper - ci_zero$boot.lower) < 0.5
)

# Replicates near pi cross the branch cut but remain close on the circle.
t0_pos <- c(phase_x = 3.0)
raw_pos <- rnorm(2000, mean = 3.0, sd = 0.15)
tmat_pos <- matrix(wrap(raw_pos), ncol = 1, dimnames = list(NULL, "phase_x"))
stopifnot(any(tmat_pos[, 1] < 0))
ci_pos <- boot_ci(t0_pos, tmat_pos, type = "perc")
width_pos <- (pi - ci_pos$boot.lower) + (ci_pos$boot.upper - (-pi))
stopifnot(width_pos < 1.0)

# The equivalent case near -pi should have a comparable width.
t0_neg <- c(phase_x = -3.0)
raw_neg <- rnorm(2000, mean = -3.0, sd = 0.15)
tmat_neg <- matrix(wrap(raw_neg), ncol = 1, dimnames = list(NULL, "phase_x"))
ci_neg <- boot_ci(t0_neg, tmat_neg, type = "perc")
width_neg <- (pi - ci_neg$boot.lower) + (ci_neg$boot.upper - (-pi))
stopifnot(abs(width_neg - width_pos) < 0.3)

# A linear quantile interval would be much wider at the branch cut.
naive_width <- diff(as.numeric(stats::quantile(tmat_pos[, 1], c(0.025, 0.975))))
stopifnot(naive_width > 3 * width_pos)

# Adding a full turn must not change the interval.
t0_shifted <- c(phase_x = wrap(3.0 + 2 * pi))
tmat_shifted <- matrix(wrap(raw_pos + 2 * pi), ncol = 1, dimnames = list(NULL, "phase_x"))
ci_shifted <- boot_ci(t0_shifted, tmat_shifted, type = "perc")
stopifnot(
  isTRUE(all.equal(ci_shifted$boot.lower, ci_pos$boot.lower)),
  isTRUE(all.equal(ci_shifted$boot.upper, ci_pos$boot.upper))
)

# With comparable spread, circular distance from zero determines the p-value.
p_pos <- boot_pvalues(t0_pos, tmat_pos)
p_zero <- boot_pvalues(t0_zero, tmat_zero)
stopifnot(p_pos < p_zero)

# Linear parameters retain their ordinary percentile interval.
t0_lin <- c(mesor_x = 3.0)
tmat_lin <- matrix(
  rnorm(500, mean = 3.0, sd = 0.15),
  ncol = 1,
  dimnames = list(NULL, "mesor_x")
)
ci_lin <- boot_ci(t0_lin, tmat_lin, type = "perc")
qs_lin <- stats::quantile(tmat_lin[, 1], c(0.025, 0.975), names = FALSE)
stopifnot(isTRUE(all.equal(as.numeric(ci_lin$boot.lower), qs_lin[1])))
