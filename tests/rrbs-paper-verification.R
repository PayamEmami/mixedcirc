library(mixedcirc)

# Verify, without changing it, the paired methylated/unmethylated construction
# described by Chen et al. (2017): sample columns condition on total coverage,
# while methylation and harmonic effects are active on the methylated channel.
simulation <- mixedcirc_simulate_data(
  n_subjects_per_group = 4,
  n_timepoints = 8,
  n_features = 5,
  groups = c("control", "treated"),
  rrbs = TRUE,
  coverage_range = c(30, 100),
  seed = 1911
)

rownames(simulation$data_input) <- as.vector(rbind(
  paste0(simulation$meta_input$replicate_id, "_methylated"),
  paste0(simulation$meta_input$replicate_id, "_unmethylated")
))
qc <- mixedcirc_rrbs_qc(
  simulation$data_input,
  simulation$meta_input,
  replicate_id = "replicate_id",
  min_coverage = 8
)
stopifnot(
  inherits(qc, "mixedcirc_rrbs_qc"),
  qc$pairing$n_pairs == nrow(simulation$meta_input),
  qc$pairing$n_features == ncol(simulation$data_input),
  isTRUE(qc$pairing$row_names_verify_pairing),
  qc$pairing$counts_are_integer,
  all(qc$feature$mean_methylation >= 0 & qc$feature$mean_methylation <= 1)
)
for (type in c("feature_coverage", "sample_coverage", "methylation")) {
  plot <- mixedcirc_plot_rrbs_qc(qc, type = type)
  stopifnot(
    inherits(plot, "ggplot"),
    length(ggplot2::ggplot_build(plot)$data) > 0L
  )
}

fit <- suppressWarnings(mixedcirc_detect(
  simulation$data_input,
  simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = simulation$obs_weights,
  lm_method = "lm",
  random_effect_type = "none"
))
object <- fit[[1L]]
design <- object@exp_design
matrix_used <- stats::model.matrix(object@fit)

stopifnot(
  identical(as.numeric(design$scaler), rep(c(1, 0), nrow(simulation$meta_input))),
  nrow(matrix_used) == 2L * nrow(simulation$meta_input),
  sum(grepl("replicate_id", colnames(matrix_used), fixed = TRUE)) ==
    nrow(simulation$meta_input),
  any(grepl("scaler", colnames(matrix_used), fixed = TRUE)),
  any(grepl("scaler", colnames(matrix_used)) &
        grepl("inphase", colnames(matrix_used))),
  any(grepl("scaler", colnames(matrix_used)) &
        grepl("outphase", colnames(matrix_used))),
  all(c("beta_cos", "beta_sin", "p_differential_global") %in%
        names(mixedcirc_results(fit, "estimates")))
)

# The mixedcirc matrix and edgeR::modelMatrixMeth() have the same column space:
# they are equivalent parameterizations of the cited paper's model.
if (requireNamespace("edgeR", quietly = TRUE)) {
  sample_design <- simulation$meta_input
  sample_design$group <- factor(sample_design$group, levels = c("control", "treated"))
  sample_design$inphase <- cos(2 * pi * sample_design$time / 24)
  sample_design$outphase <- sin(2 * pi * sample_design$time / 24)
  treatment_design <- stats::model.matrix(
    ~ group + group:inphase + group:outphase,
    sample_design
  )
  paper_matrix <- edgeR::modelMatrixMeth(treatment_design)
  rank_paper <- qr(paper_matrix)$rank
  rank_mixedcirc <- qr(matrix_used)$rank
  rank_combined <- qr(cbind(paper_matrix, matrix_used))$rank
  stopifnot(
    nrow(paper_matrix) == nrow(matrix_used),
    ncol(paper_matrix) == ncol(matrix_used),
    rank_paper == rank_mixedcirc,
    rank_combined == rank_paper
  )
}

# Equal methylated and unmethylated values create no rhythmic methylation
# contrast even though paired coverage varies over time.
set.seed(1912)
n <- nrow(simulation$meta_input)
coverage_signal <- round(60 + 20 * cos(2 * pi * simulation$meta_input$time / 24))
null_pair <- as.vector(rbind(coverage_signal / 2, coverage_signal / 2))
null_data <- cbind(null_methylation = null_pair)
null_weights <- matrix(1, nrow = 2L * n, ncol = 1L)
null_fit <- suppressWarnings(mixedcirc_detect(
  null_data,
  simulation$meta_input,
  time = "time",
  group = "group",
  id = "id",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = null_weights,
  lm_method = "lm",
  random_effect_type = "none"
))
null_estimates <- mixedcirc_results(null_fit, "estimates", adjust_method = NULL)
stopifnot(
  all(abs(null_estimates$beta_cos) < 1e-8),
  all(abs(null_estimates$beta_sin) < 1e-8),
  all(null_estimates$amplitude < 1e-8)
)
