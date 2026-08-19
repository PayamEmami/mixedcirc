# mixedcirc

`mixedcirc` is a high-throughput harmonic-regression toolkit for circadian and
differential-circadian analysis. It supports ordinary linear models,
subject-level random effects for repeated measurements, observation weights,
paired bisulfite-sequencing counts, bootstrap uncertainty, direct joint rhythm
contrasts, model diagnostics, subject-specific rhythms, and
publication-oriented plots.

The core model represents time with cosine and sine terms. For a grouped
analysis, each group receives its own harmonic coefficients. Rhythmicity and
between-group changes are tested jointly over both coefficients, avoiding the
invalid shortcut of declaring a change because one group is significant and
another is not.



## Installation

Install the Bioconductor dependencies, then install `mixedcirc` from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Biobase", "BiocParallel", "limma"
))

install.packages("remotes")
remotes::install_github("mixedcirc/mixedcirc")
```

Optional packages such as `edgeR`, `lmerTest`, `performance`, and `rmcorr`
enable count input, additional mixed-model inference, model summaries, and
repeated-measures correlation, respectively.

## Quick start

```r
library(mixedcirc)

sim <- mixedcirc_simulate_data(
  n_subjects_per_group = 4,
  n_timepoints = 8,
  n_features = 5,
  groups = c("control", "treated"),
  seed = 1
)

fit <- mixedcirc_detect(
  data_input = sim$data_input,
  meta_input = sim$meta_input,
  time = "time",
  group = "group",
  id = "id",
  period = 24,
  lm_method = "lme",
  random_effect_type = "mesor"
)

fit[[1]]@results
plot(fit[[1]])

# Analysis-ready estimates, direct contrasts, diagnostics, and overview.
estimates <- mixedcirc_results(fit, "estimates")
contrasts <- mixedcirc_contrasts(fit)
diagnostics <- mixedcirc_results(fit, "diagnostics")
summary(fit)

# Dataset-level publication plots.
mixedcirc_plot_heatmap(fit)
mixedcirc_plot_differential(fit)
mixedcirc_plot_polar(fit)
```

For independent observations, use `lm_method = "lm"` and
`random_effect_type = "none"`. Additional nuisance variables can be included
with a one-sided formula such as `formula_extra = ~ age + batch`.

## Phase and peak time

The returned phase is the signed shift \(\phi\) in

\[
M + A\cos(2\pi t / period + \phi).
\]

It is represented on `[-2*pi, 0)`; a rhythm peaking at time 6 in a 24-hour
cycle therefore has phase `-pi/2`. Convert phases to peak times with:

```r
phase_columns <- grep("^phase_", names(fit[[1]]@results), value = TRUE)
phases <- unlist(fit[[1]]@results[phase_columns], use.names = FALSE)
mixedcirc_detectpeaks(phases, period = 24)
```

Plotting and interpolation use the period stored by `mixedcirc_detect()` when
`period` is omitted.

## Main functions

- `mixedcirc_detect()` fits feature-wise rhythmic models.
- `mixedcirc_results()`, `mixedcirc_contrasts()`, and `mixedcirc_summary()`
  produce analysis-ready result tables and overview statistics.
- `mixedcirc_classify()` combines within-group rhythmicity and direct
  between-group tests into conservative gain, loss, altered, unchanged,
  arrhythmic, or uncertain classes. Optional effect-size thresholds can require
  a minimum amplitude change, peak-time shift, or harmonic distance.
- `mixedcirc_boottest()` estimates parametric-bootstrap uncertainty.
- `mixedcirc_adjust()` adjusts selected p-value columns across features.
- `mixedcirc_detectpeaks()` converts signed phase to peak time.
- `mixedcirc_fit_plot()`, `mixedcirc_interpolate()`, and the
  `mixedcirc_plot_*()` family generate fitted curves and publication figures.
- `mixedcirc_subject_rhythms()` reconstructs individual rhythms from random
  effects.
- `mixedcirc_rrbs_qc()` verifies and summarizes paired bisulfite counts without
  modifying them.
- `mixedcirc_fitted()` extracts fitted values aligned to input rows.
- `mixedcirc_remove_trend()` creates a residualized input matrix;
  `mixedcirc_detrend()` removes a linear trend from stored fitted curves for
  descriptive downstream use.
- `mixedcirc_add_lag()` and `mixedcirc_compare()` support lagged association,
  Pearson correlation, and repeated-measures correlation.
- `mixedcirc_simulate_data()` simulates continuous or RRBS-like workflows.
- `mixedcirc_integrate()` performs experimental supervised SGCCA integration.

## High-throughput operation

Use `ncores` for feature-parallel fitting. For large studies, set
`store_models = FALSE` to retain estimates, tests, diagnostics, and model
quality statistics while omitting the fitted model objects. Prediction and
bootstrap functions require `store_models = TRUE` and report a direct refitting
instruction otherwise.

By default, a failed feature stops the analysis. Set `on_error = "keep"` to
retain a diagnostic record for the failure and continue fitting the remaining
features:

```r
fit_large <- mixedcirc_detect(
  data_input = large_matrix,
  meta_input = metadata,
  time = "time",
  group = "condition",
  id = "subject",
  ncores = 4,
  store_models = FALSE,
  on_error = "keep"
)

mixedcirc_results(fit_large, "diagnostics")
```

## RNA-seq counts

Supply raw RNA-seq counts as an `edgeR::DGEList`, with genes in rows and samples
in columns. Filter unexpressed and very low-count genes before fitting, for
example with `edgeR::filterByExpr()`. `mixedcirc_detect()` converts the counts
to log-CPM values and estimates voom-style mean-variance precision weights. For
repeated measurements (`lm_method = "lme"`), the weight estimation is
implemented locally in mixedcirc and does not require the `variancePartition`
package.

```r
keep <- edgeR::filterByExpr(dge, group = metadata$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]

rna_fit <- mixedcirc_detect(
  dge,
  metadata,
  time = "time",
  group = "condition",
  id = "subject",
  lm_method = "lme"
)
```

## Paired bisulfite-sequencing counts

The RRBS/BS-seq analysis follows the paired methylated/unmethylated design
described by Chen et al. (2017), in which sample-specific columns condition on
coverage and methylation effects act on the methylated-count indicator. This
pathway is part of the main package, not the experimental integration feature.

Input rows must alternate methylated then unmethylated counts for every
biological sample. `meta_input` has one row per biological sample and a unique
`replicate_id`. Supply positive observation weights or request internal weight
estimation. Run QC before fitting:

```r
qc <- mixedcirc_rrbs_qc(
  rrbs_counts,
  rrbs_metadata,
  replicate_id = "replicate_id",
  min_coverage = 8
)
print(qc)
mixedcirc_plot_rrbs_qc(qc, "feature_coverage")

rrbs_fit <- mixedcirc_detect(
  rrbs_counts,
  rrbs_metadata,
  time = "time",
  group = "condition",
  id = "subject",
  RRBS = TRUE,
  replicate_id = "replicate_id",
  obs_weights = rrbs_weights
)
```

For a very large data set already stored as a `methylKit::methylBaseDB`,
`mixedcirc_rrbs_voom_mixed()` estimates log-CPM values and precision weights
in tabix chunks. This compatibility workflow is experimental and requires the
optional Bioconductor package `methylKit`. The optional
`mixed_circ_tss_calculation()` and `construct_introns()` provide optional
BED/GFF annotation for assigning RRBS loci to promoters, TSSes, exons, and
introns.

Reference: Chen Y, Pal B, Visvader JE, Smyth GK. *Differential methylation
analysis of reduced representation bisulfite sequencing experiments using
edgeR*. F1000Research (2017), <https://pmc.ncbi.nlm.nih.gov/articles/PMC5747346/>.

## Important interpretation notes

- `data_input` matrices and data frames use observations in rows and features
  in columns. A `DGEList` follows the Bioconductor convention: features in rows
  and samples in columns.
- RRBS input uses two consecutive rows per biological observation: methylated,
  then unmethylated. `mixedcirc_rrbs_qc()` checks the representation but never
  reorders, filters, or transforms the input.
- `mixedcirc_compare()` is associative, not causal. Lag selection over the same
  data is exploratory unless validated on independent data or inside a complete
  resampling procedure. Lag values count ordered observations rather than
  elapsed time units.
- If a linear time trend is part of the scientific model, include it alongside
  the harmonic terms with `formula_extra = ~ time`. Two-stage detrending is
  intended for workflows that specifically require a residualized matrix.
- `mixedcirc_integrate()` is experimental and supervised by the supplied time
  and group design. Scores learned and tested on the same samples are not
  independent confirmatory evidence.

See `vignette("mixedcirc", package = "mixedcirc")` and the individual function
help pages for complete examples.

## License

[Apache License 2.0](LICENSE.md)

Maintainer: [Payam Emami](mailto:payam.emami@scilifelab.se)
