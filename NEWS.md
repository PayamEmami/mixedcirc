# mixedcirc 0.2.2

- Phase uncertainty from `mixedcirc_boottest()` now respects circular geometry,
  including intervals that cross the `-pi`/`pi` boundary. The corresponding
  plot draws these intervals locally around their point estimates.
- Snow-cluster bootstrap workers now load the package namespaces and S4 classes
  required for both ordinary and mixed-model refits.
- Documentation now uses the correct conversion from phase difference to
  peak-time difference, `peak = -phase / omega`.
- Lagged simulations now report harmonic truth values using observation-step
  lags and the package's signed phase convention.
- The polar plot legend uses the reader-facing labels `Rhythmic: Yes/No`.
- The main vignette now presents the continuous, RNA-seq, paired-count,
  inference, plotting, uncertainty, and high-throughput workflows together.
- Package and citation links now point to
  <https://github.com/mixedcirc/mixedcirc>.

# mixedcirc 0.2.1

- Fitted-curve plots use readable x-axis breaks in dense sampling designs.
- Polar plots label each feature once rather than once per group.
- `mixedcirc_plot_design()` can display one row per subject, making repeated
  measurements visible in the sampling design.

# mixedcirc 0.2.0

- Added analysis-ready estimate, contrast, diagnostic, classification, and
  analysis-summary tables with within-family multiple-testing correction.
- Added an omnibus differential-rhythm test for grouped analyses, distinct from
  the omnibus test that any group is rhythmic.
- Added heatmaps, differential-effect and polar plots, temporal-design plots,
  model diagnostics, subject-specific rhythm plots, bootstrap interval plots,
  and paired-count QC plots.
- Added subject-specific mesor, amplitude, phase, and peak-time estimates from
  mixed-effects fits.
- Added compact results (`store_models = FALSE`) and feature-level failure
  retention (`on_error = "keep"`) for large analyses.
- Added RRBS input QC and tests of the paired-count design against
  `edgeR::modelMatrixMeth()`.
- Added experimental chunked `methylBaseDB` weight estimation and optional
  genomic annotation helpers.

# mixedcirc 0.1.0

- Added feature-wise harmonic regression with ordinary and mixed-effects
  models, grouped rhythmicity tests, signed phase estimates, and peak-time
  conversion.
- Added parametric bootstrap inference for `lm` and `lmerMod` fits.
- Added weighted continuous, `DGEList`, and paired RRBS workflows.
- Added fitted-value extraction, detrending helpers, interpolation, plotting,
  p-value adjustment, lagged comparison, and data simulation.
- Added package-local feature-wise model fitting and mixed-model voom-style
  precision weights.
- Added experimental supervised SGCCA integration.
