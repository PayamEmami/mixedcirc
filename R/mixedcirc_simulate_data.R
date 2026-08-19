#' Simulate circadian data for mixedcirc
#'
#' Simulates continuous or RRBS-like circadian data with optional covariates and
#' optional copied lagged relationships between selected features.
#'
#' The function is designed for testing and demonstrating `mixedcirc_detect()`,
#' `mixedcirc_compare()`, and related package functionality.
#'
#' @param n_subjects_per_group Number of subjects per group. Default `8`.
#' @param n_timepoints Number of repeated time points per subject. Default `6`.
#' @param n_features Number of features to simulate. Default `20`.
#' @param groups Character vector of group labels. Default `c("A")`.
#' @param include_covariates Logical. If `TRUE`, subject-level covariates
#'   (`age`, `batch`) are simulated and added to the signal. Default `FALSE`.
#' @param rrbs Logical. If `TRUE`, simulate RRBS-like methylated/unmethylated
#'   counts. If `FALSE`, simulate continuous data. Default `FALSE`.
#' @param seed Random seed. Default `123`.
#'
#' @param age_mean Mean age for the first group when covariates are simulated.
#' @param age_group_shift Shift in mean age added for each successive group.
#' @param age_sd Standard deviation of age.
#' @param batch_levels Character vector of batch levels.
#'
#' @param age_beta Default age effect used when covariates are included.
#'   If `NULL`, a default value is chosen depending on `rrbs`.
#' @param batch_effects Named numeric vector of batch effects.
#'   If `NULL`, defaults are chosen depending on `rrbs`.
#' @param age_beta_sd Standard deviation of feature-specific jitter added to
#'   `age_beta`. Default `0`.
#' @param batch_effect_sd Standard deviation of feature-specific jitter added to
#'   each batch effect. Default `0`.
#'
#' @param subject_random_intercept_sd Standard deviation of subject-specific
#'   random intercepts. Default `0.6`.
#' @param residual_sd Residual standard deviation for continuous data.
#'   Default `0.8`.
#' @param coverage_range Numeric vector of length 2 giving the minimum and
#'   maximum sequencing depth for RRBS simulation. Default `c(15, 80)`.
#'
#' @param cross_lag_spec Optional list describing copied lagged relationships
#'   between features. Each element must itself be a list with:
#'   \describe{
#'     \item{from}{Source feature index or feature name.}
#'     \item{to}{Target feature index or feature name.}
#'     \item{lag}{Integer lag in observation steps, or a named integer vector
#'     with one lag per group.}
#'     \item{beta}{Numeric scaling factor applied to the copied lagged source.
#'     Default behavior is multiplication of the shifted source by `beta`.}
#'   }
#'
#'   Positive lag means that the target at one observation is copied from an
#'   earlier source observation within the same subject; a negative lag uses a
#'   later source observation. Lags refer to positions in the time-ordered data,
#'   not elapsed time units.
#'
#'   In this implementation, the target feature is constructed from the observed
#'   source feature:
#'   \itemize{
#'     \item for continuous data, after residual noise has been added
#'     \item for RRBS data, after methylated / unmethylated counts have been generated
#'   }
#'   Therefore, `lag = 0` and `beta = 1` produce an exact observed-data copy.
#' @param period Positive rhythm period used for the time grid and simulated
#'   harmonic signal. Default `24`.
#'
#' @return
#' A list with:
#' \describe{
#'   \item{data_input}{Simulated feature matrix or data.frame.}
#'   \item{meta_input}{Sample-level metadata.}
#'   \item{obs_weights}{Observation weights for RRBS data, otherwise `NULL`.}
#'   \item{truth}{A data frame describing the generating parameters. Phase
#'   columns use the same signed model-shift convention as `mixedcirc_detect()`.
#'   For copied features with a nonzero lag, the harmonic entries describe the
#'   intended shifted signal; boundary observations are copied without a shift
#'   because no source observation exists beyond the sampled series.}
#' }
#'
#' @details
#' For continuous data, each feature is generated from a group-specific circadian
#' pattern plus subject random intercepts and optional covariate effects.
#'
#' For RRBS-like data, each feature is generated on a latent logit scale and then
#' converted to methylated/unmethylated counts using a binomial model.
#'
#' If `cross_lag_spec` is provided, selected target features are replaced by
#' lagged copies of source features.
#'
#' @examples
#' sim <- mixedcirc_simulate_data(
#'   n_subjects_per_group = 10,
#'   n_timepoints = 8,
#'   n_features = 5,
#'   groups = c("A", "B"),
#'   cross_lag_spec = list(
#'     list(from = 1, to = 2, lag = 1, beta = 1.0),
#'     list(from = 3, to = 4, lag = c(A = 1, B = 2), beta = 0.7)
#'   ),
#'   seed = 1
#' )
#'
#' @export
mixedcirc_simulate_data <- function(
    n_subjects_per_group = 8,
    n_timepoints = 6,
    n_features = 20,
    groups = c("A"),
    include_covariates = FALSE,
    rrbs = FALSE,
    seed = 123,
    age_mean = 45,
    age_group_shift = 5,
    age_sd = 8,
    batch_levels = c("B1", "B2"),
    age_beta = NULL,
    batch_effects = NULL,
    age_beta_sd = 0,
    batch_effect_sd = 0,
    subject_random_intercept_sd = 0.6,
    residual_sd = 0.8,
    coverage_range = c(15, 80),
    cross_lag_spec = NULL,
    period = 24
) {
  groups <- as.character(groups)
  n_groups <- length(groups)

  if (length(groups) < 1 || anyNA(groups) || any(groups == "") || anyDuplicated(groups)) {
    stop("groups must contain unique, non-missing, non-empty group names")
  }
  if (!is.logical(include_covariates) || length(include_covariates) != 1 || is.na(include_covariates)) {
    stop("include_covariates must be TRUE or FALSE")
  }
  if (!is.logical(rrbs) || length(rrbs) != 1 || is.na(rrbs)) stop("rrbs must be TRUE or FALSE")
  if (!is.numeric(n_subjects_per_group) || length(n_subjects_per_group) != 1 ||
      !is.finite(n_subjects_per_group) || n_subjects_per_group < 1 ||
      n_subjects_per_group %% 1 != 0) {
    stop("n_subjects_per_group must be a positive integer")
  }
  if (!is.numeric(n_timepoints) || length(n_timepoints) != 1 ||
      !is.finite(n_timepoints) || n_timepoints < 2 || n_timepoints %% 1 != 0) {
    stop("n_timepoints must be an integer >= 2")
  }
  if (!is.numeric(n_features) || length(n_features) != 1 ||
      !is.finite(n_features) || n_features < 1 || n_features %% 1 != 0) {
    stop("n_features must be a positive integer")
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed %% 1 != 0) {
    stop("seed must be one finite integer")
  }
  if (seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be between 0 and .Machine$integer.max")
  }
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0) {
    stop("period must be one positive finite number")
  }
  scalar_finite <- list(
    age_mean = age_mean,
    age_group_shift = age_group_shift
  )
  bad_scalar <- !vapply(
    scalar_finite,
    function(z) is.numeric(z) && length(z) == 1L && is.finite(z),
    logical(1)
  )
  if (any(bad_scalar)) {
    stop(names(scalar_finite)[which(bad_scalar)[1L]], " must be one finite number")
  }
  scalar_nonnegative <- list(
    age_sd = age_sd,
    age_beta_sd = age_beta_sd,
    batch_effect_sd = batch_effect_sd,
    subject_random_intercept_sd = subject_random_intercept_sd,
    residual_sd = residual_sd
  )
  bad_sd <- !vapply(
    scalar_nonnegative,
    function(z) is.numeric(z) && length(z) == 1L && is.finite(z) && z >= 0,
    logical(1)
  )
  if (any(bad_sd)) {
    stop(names(scalar_nonnegative)[which(bad_sd)[1L]], " must be one non-negative finite number")
  }
  if (!is.numeric(coverage_range) || length(coverage_range) != 2 ||
      any(!is.finite(coverage_range)) || any(coverage_range %% 1 != 0) ||
      coverage_range[1] <= 0 || coverage_range[2] < coverage_range[1]) {
    stop("coverage_range must be a numeric vector of length 2 like c(15, 80)")
  }
  n_subjects_per_group <- as.integer(n_subjects_per_group)
  n_timepoints <- as.integer(n_timepoints)
  n_features <- as.integer(n_features)
  coverage_range <- as.integer(coverage_range)
  set.seed(as.integer(seed))

  if (include_covariates) {
    if (!is.character(batch_levels) || length(batch_levels) < 1L ||
        anyNA(batch_levels) || any(batch_levels == "") || anyDuplicated(batch_levels)) {
      stop("batch_levels must contain unique, non-missing, non-empty names")
    }
    if (is.null(age_beta)) age_beta <- if (rrbs) 0.01 else 0.03
    if (!is.numeric(age_beta) || length(age_beta) != 1L || !is.finite(age_beta)) {
      stop("age_beta must be NULL or one finite number")
    }
    if (is.null(batch_effects)) {
      batch_effects <- if (rrbs) c(B1 = 0, B2 = 0.4) else c(B1 = 0, B2 = 1.2)
    }
    if (!is.numeric(batch_effects) || length(batch_effects) < 1L ||
        any(!is.finite(batch_effects)) || is.null(names(batch_effects)) ||
        anyNA(names(batch_effects)) || any(names(batch_effects) == "") ||
        anyDuplicated(names(batch_effects))) {
      stop("batch_effects must be a named numeric vector, e.g. c(B1 = 0, B2 = 3)")
    }
    if (!all(batch_levels %in% names(batch_effects))) {
      stop("All batch_levels must appear in names(batch_effects)")
    }
  }

  bind_rows_fill <- function(xlist) {
    all_cols <- unique(unlist(lapply(xlist, names)))
    out <- lapply(xlist, function(x) {
      miss <- setdiff(all_cols, names(x))
      if (length(miss) > 0) {
        for (m in miss) x[[m]] <- NA
      }
      x[all_cols]
    })
    out <- do.call(rbind, lapply(out, function(x) as.data.frame(x, check.names = FALSE)))
    rownames(out) <- NULL
    out
  }

  safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

  calc_phase_mixedcirc <- function(inphase, outphase) {
    if (!is.finite(inphase) || !is.finite(outphase) ||
        (inphase == 0 && outphase == 0)) {
      return(NA_real_)
    }
    phi <- atan2(-outphase, inphase)
    if (phi >= 0) phi <- phi - 2 * pi
    phi
  }

  calc_amp_phase <- function(beta_cos, beta_sin) {
    amp <- sqrt(beta_cos^2 + beta_sin^2)
    phi <- calc_phase_mixedcirc(beta_cos, beta_sin)
    list(amplitude = amp, phase = phi)
  }

  shift_vector <- function(x, lag_k) {
    n <- length(x)
    out <- rep(NA_real_, n)

    if (lag_k > 0) {
      if (n > lag_k) out[(lag_k + 1):n] <- x[1:(n - lag_k)]
    } else if (lag_k < 0) {
      kk <- abs(lag_k)
      if (n > kk) out[1:(n - kk)] <- x[(kk + 1):n]
    } else {
      out <- x
    }

    out
  }

  resolve_feature_index <- function(z, feature_names) {
    if (length(z) != 1) stop("Feature selector must have length 1.")
    if (is.numeric(z)) {
      if (!is.finite(z) || z %% 1 != 0) {
        stop("Numeric feature index must be a finite integer.")
      }
      z <- as.integer(z)
      if (is.na(z) || z < 1 || z > length(feature_names)) stop("Numeric feature index out of range.")
      return(z)
    }
    if (is.character(z)) {
      idx <- match(z, feature_names)
      if (is.na(idx)) stop("Unknown feature name in cross_lag_spec: ", z)
      return(idx)
    }
    stop("Feature selector must be numeric index or character feature name.")
  }

  normalize_cross_lag_spec <- function(cross_lag_spec, feature_names, groups) {
    if (is.null(cross_lag_spec)) return(list())
    if (!is.list(cross_lag_spec)) stop("cross_lag_spec must be NULL or a list.")

    out <- vector("list", length(cross_lag_spec))
    to_seen <- integer(0)

    for (i in seq_along(cross_lag_spec)) {
      spec_i <- cross_lag_spec[[i]]
      if (!is.list(spec_i)) stop("Each element of cross_lag_spec must be a list.")

      required_names <- c("from", "to", "lag", "beta")
      if (!all(required_names %in% names(spec_i))) {
        stop("Each cross_lag_spec element must contain: from, to, lag, beta")
      }

      from_idx <- resolve_feature_index(spec_i$from, feature_names)
      to_idx   <- resolve_feature_index(spec_i$to, feature_names)

      if (from_idx == to_idx) stop("cross_lag_spec cannot have the same feature in 'from' and 'to'.")
      if (to_idx %in% to_seen) stop("Each target feature may appear only once in cross_lag_spec.")
      to_seen <- c(to_seen, to_idx)

      if (!is.numeric(spec_i$beta) || length(spec_i$beta) != 1 || !is.finite(spec_i$beta)) {
        stop("cross_lag_spec[[", i, "]]$beta must be a single finite numeric value.")
      }

      lag_i <- spec_i$lag
      if (length(lag_i) == 1) {
        if (!is.numeric(lag_i) || !is.finite(lag_i) || lag_i %% 1 != 0) {
          stop("Each lag must be an integer or a named integer vector by group.")
        }
        lag_by_group <- setNames(rep(as.integer(lag_i), length(groups)), groups)
      } else {
        if (is.null(names(lag_i))) {
          stop("If lag has length > 1, it must be a named vector with names matching groups.")
        }
        if (!all(groups %in% names(lag_i))) {
          stop("All groups must appear in names(lag) for grouped lag specification.")
        }
        lag_selected <- lag_i[groups]
        if (!is.numeric(lag_selected) || any(!is.finite(lag_selected)) ||
            any(lag_selected %% 1 != 0)) {
          stop("All group-specific lags must be finite integers.")
        }
        lag_by_group <- as.integer(lag_selected)
        names(lag_by_group) <- groups
      }

      out[[i]] <- list(
        from = from_idx,
        to = to_idx,
        beta = as.numeric(spec_i$beta),
        lag_by_group = lag_by_group
      )
    }

    out
  }

  shift_phase <- function(phi, lag_k) {
    phase_shift <- 2 * pi * lag_k / n_timepoints
    wrapped <- (phi - phase_shift) %% (2 * pi)
    if (wrapped == 0) -2 * pi else wrapped - 2 * pi
  }

  # copy observed continuous feature matrix
  apply_copy_lag_matrix <- function(mat, meta_input, groups, cross_specs) {
    if (length(cross_specs) == 0) return(mat)

    out <- mat
    ids <- unique(meta_input$id)
    id_group_map <- stats::setNames(meta_input$group[match(ids, meta_input$id)], ids)

    for (spec_i in cross_specs) {
      from_j <- spec_i$from
      to_j   <- spec_i$to
      beta_j <- spec_i$beta
      target_vec <- rep(NA_real_, nrow(out))

      for (g in groups) {
        lag_g <- spec_i$lag_by_group[g]
        ids_g <- ids[id_group_map[ids] == g]

        for (id_i in ids_g) {
          idx <- which(meta_input$id == id_i & meta_input$group == g)
          ord <- order(meta_input$time[idx])
          idx_ord <- idx[ord]

          source_vals <- out[idx_ord, from_j]
          shifted_source <- shift_vector(source_vals, lag_g)
          shifted_source[is.na(shifted_source)] <- source_vals[is.na(shifted_source)]

          target_vec[idx_ord] <- beta_j * shifted_source
        }
      }

      out[, to_j] <- target_vec
    }

    out
  }

  # copy observed RRBS counts feature-wise
  apply_copy_lag_rrbs_counts <- function(M_mat, U_mat, meta_input, groups, cross_specs) {
    if (length(cross_specs) == 0) {
      return(list(M = M_mat, U = U_mat))
    }

    M_out <- M_mat
    U_out <- U_mat
    ids <- unique(meta_input$id)
    id_group_map <- stats::setNames(meta_input$group[match(ids, meta_input$id)], ids)

    for (spec_i in cross_specs) {
      from_j <- spec_i$from
      to_j   <- spec_i$to
      beta_j <- spec_i$beta

      if (!isTRUE(all.equal(beta_j, 1))) {
        warning("For RRBS copying, beta != 1 changes methylated/unmethylated counts by rescaling and rounding. Exact copying only occurs for beta = 1.")
      }

      target_M <- rep(NA_real_, nrow(M_out))
      target_U <- rep(NA_real_, nrow(U_out))

      for (g in groups) {
        lag_g <- spec_i$lag_by_group[g]
        ids_g <- ids[id_group_map[ids] == g]

        for (id_i in ids_g) {
          idx <- which(meta_input$id == id_i & meta_input$group == g)
          ord <- order(meta_input$time[idx])
          idx_ord <- idx[ord]

          src_M <- M_out[idx_ord, from_j]
          src_U <- U_out[idx_ord, from_j]

          sh_M <- shift_vector(src_M, lag_g)
          sh_U <- shift_vector(src_U, lag_g)

          sh_M[is.na(sh_M)] <- src_M[is.na(sh_M)]
          sh_U[is.na(sh_U)] <- src_U[is.na(sh_U)]

          if (beta_j != 1) {
            sh_M <- pmax(0, round(beta_j * sh_M))
            sh_U <- pmax(0, round(beta_j * sh_U))
          }

          target_M[idx_ord] <- sh_M
          target_U[idx_ord] <- sh_U
        }
      }

      M_out[, to_j] <- target_M
      U_out[, to_j] <- target_U
    }

    list(M = M_out, U = U_out)
  }

  # Build metadata
  meta_list <- vector("list", n_groups)
  subj_counter <- 1

  for (g in seq_along(groups)) {
    ids <- paste0("id", subj_counter:(subj_counter + n_subjects_per_group - 1))
    subj_counter <- subj_counter + n_subjects_per_group
    time_vec <- seq(0, period - period / n_timepoints, length.out = n_timepoints)

    meta_g <- expand.grid(
      id = ids,
      time = time_vec,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    meta_g$group <- groups[g]
    meta_g$replicate_id <- paste0(meta_g$id, "_t", gsub("\\.", "_", as.character(meta_g$time)))

    if (include_covariates) {
      age_map <- setNames(
        round(rnorm(length(ids), mean = age_mean + age_group_shift * (g - 1), sd = age_sd)),
        ids
      )
      batch_map <- setNames(sample(batch_levels, length(ids), replace = TRUE), ids)

      meta_g$age <- age_map[meta_g$id]
      meta_g$batch <- batch_map[meta_g$id]
    }

    meta_list[[g]] <- meta_g
  }

  meta_input <- do.call(rbind, meta_list)
  rownames(meta_input) <- NULL

  inphase <- cos(2 * pi * meta_input$time / period)
  outphase <- sin(2 * pi * meta_input$time / period)
  n_samples <- nrow(meta_input)

  subject_ids <- unique(meta_input$id)
  rand_intercepts <- rnorm(length(subject_ids), 0, subject_random_intercept_sd)
  names(rand_intercepts) <- subject_ids

  feature_names <- if (!rrbs) paste0("feature_", seq_len(n_features)) else paste0("cpg_", seq_len(n_features))

  cross_specs <- normalize_cross_lag_spec(
    cross_lag_spec = cross_lag_spec,
    feature_names = feature_names,
    groups = groups
  )

  truth_list <- vector("list", n_features)

  if (!rrbs) {
    data_input <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
    colnames(data_input) <- feature_names

    mu_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
    colnames(mu_mat) <- feature_names

    feature_param_store <- vector("list", n_features)

    for (j in seq_len(n_features)) {
      base_mesor <- rnorm(1, mean = 10, sd = 2)
      group_mesor_shift <- setNames(rnorm(n_groups, 0, 1.0), groups)
      group_beta_cos <- setNames(rnorm(n_groups, 1.5, 0.5), groups)
      group_beta_sin <- setNames(rnorm(n_groups, 0.0, 0.5), groups)

      if (include_covariates) {
        age_beta_j <- rnorm(1, mean = age_beta, sd = age_beta_sd)
        batch_effects_j <- batch_effects
        for (nm in names(batch_effects_j)) {
          batch_effects_j[nm] <- rnorm(1, mean = batch_effects[nm], sd = batch_effect_sd)
        }
      } else {
        age_beta_j <- NA_real_
        batch_effects_j <- NULL
      }

      mu <- base_mesor +
        group_mesor_shift[meta_input$group] +
        group_beta_cos[meta_input$group] * inphase +
        group_beta_sin[meta_input$group] * outphase +
        rand_intercepts[meta_input$id]

      if (include_covariates) {
        mu <- mu +
          age_beta_j * meta_input$age +
          unname(batch_effects_j[meta_input$batch])
      }

      mu_mat[, j] <- mu

      feature_param_store[[j]] <- list(
        base_mesor = base_mesor,
        group_mesor_shift = group_mesor_shift,
        group_beta_cos = group_beta_cos,
        group_beta_sin = group_beta_sin,
        age_beta_j = age_beta_j,
        batch_effects_j = batch_effects_j
      )

      truth_row <- list(
        feature = feature_names[j],
        data_type = "continuous",
        copied_from = NA,
        copied_beta = NA,
        base_intercept = base_mesor,
        age_beta = age_beta_j,
        residual_sd = residual_sd,
        subject_random_intercept_sd = subject_random_intercept_sd
      )

      if (include_covariates) {
        for (b in batch_levels) {
          truth_row[[paste0("batch_effect_", safe_name(b))]] <- unname(batch_effects_j[b])
        }
      }

      for (g in groups) {
        gp <- safe_name(g)
        mesor_g <- base_mesor + group_mesor_shift[g]
        beta_cos_g <- group_beta_cos[g]
        beta_sin_g <- group_beta_sin[g]
        amp_phase <- calc_amp_phase(beta_cos_g, beta_sin_g)

        truth_row[[paste0("mesor_", gp)]] <- mesor_g
        truth_row[[paste0("group_shift_", gp)]] <- group_mesor_shift[g]
        truth_row[[paste0("beta_cos_", gp)]] <- beta_cos_g
        truth_row[[paste0("beta_sin_", gp)]] <- beta_sin_g
        truth_row[[paste0("amplitude_", gp)]] <- amp_phase$amplitude
        truth_row[[paste0("phase_", gp)]] <- amp_phase$phase
      }

      truth_list[[j]] <- truth_row
    }

    # First generate observed data with noise
    data_input <- mu_mat + matrix(rnorm(n_samples * n_features, 0, residual_sd), nrow = n_samples)

    # Then overwrite copied features using observed source signals
    data_input <- apply_copy_lag_matrix(
      mat = data_input,
      meta_input = meta_input,
      groups = groups,
      cross_specs = cross_specs
    )

    if (length(cross_specs) > 0) {
      for (spec_i in seq_along(cross_specs)) {
        from_j <- cross_specs[[spec_i]]$from
        to_j   <- cross_specs[[spec_i]]$to
        beta_i <- cross_specs[[spec_i]]$beta
        lag_i  <- cross_specs[[spec_i]]$lag_by_group
        from_nm <- feature_names[from_j]

        truth_list[[to_j]]$copied_from <- from_nm
        truth_list[[to_j]]$copied_beta <- beta_i

        source_params <- feature_param_store[[from_j]]

        for (g in groups) {
          gp <- safe_name(g)
          unshifted_cos <- beta_i * source_params$group_beta_cos[g]
          unshifted_sin <- beta_i * source_params$group_beta_sin[g]
          phi_new <- shift_phase(
            calc_phase_mixedcirc(unshifted_cos, unshifted_sin),
            lag_i[g]
          )
          amp_new <- abs(beta_i) * sqrt(source_params$group_beta_cos[g]^2 + source_params$group_beta_sin[g]^2)
          beta_cos_new <- amp_new * cos(phi_new)
          beta_sin_new <- -amp_new * sin(phi_new)

          truth_list[[to_j]][[paste0("copied_lag_", gp)]] <- lag_i[g]
          truth_list[[to_j]][[paste0("mesor_", gp)]] <- beta_i * (source_params$base_mesor + source_params$group_mesor_shift[g])
          truth_list[[to_j]][[paste0("group_shift_", gp)]] <- beta_i * source_params$group_mesor_shift[g]
          truth_list[[to_j]][[paste0("beta_cos_", gp)]] <- beta_cos_new
          truth_list[[to_j]][[paste0("beta_sin_", gp)]] <- beta_sin_new
          truth_list[[to_j]][[paste0("amplitude_", gp)]] <- amp_new
          truth_list[[to_j]][[paste0("phase_", gp)]] <- phi_new
        }

        if (include_covariates) {
          truth_list[[to_j]]$age_beta <- beta_i * source_params$age_beta_j
          for (b in batch_levels) {
            truth_list[[to_j]][[paste0("batch_effect_", safe_name(b))]] <-
              beta_i * unname(source_params$batch_effects_j[b])
          }
        } else {
          truth_list[[to_j]]$age_beta <- NA_real_
        }
      }
    }

    truth <- bind_rows_fill(truth_list)

    return(list(
      data_input = as.data.frame(data_input),
      meta_input = meta_input,
      obs_weights = NULL,
      truth = truth
    ))

  } else {
    data_input <- matrix(NA_real_, nrow = 2 * n_samples, ncol = n_features)
    obs_weights <- matrix(NA_real_, nrow = 2 * n_samples, ncol = n_features)
    colnames(data_input) <- feature_names
    colnames(obs_weights) <- feature_names

    eta_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
    colnames(eta_mat) <- feature_names

    feature_param_store <- vector("list", n_features)

    for (j in seq_len(n_features)) {
      base_logit <- rnorm(1, mean = -0.5, sd = 0.8)
      group_shift <- setNames(rnorm(n_groups, 0, 0.5), groups)
      group_beta_cos <- setNames(rnorm(n_groups, 0.8, 0.3), groups)
      group_beta_sin <- setNames(rnorm(n_groups, 0.0, 0.3), groups)

      if (include_covariates) {
        age_beta_j <- rnorm(1, mean = age_beta, sd = age_beta_sd)
        batch_effects_j <- batch_effects
        for (nm in names(batch_effects_j)) {
          batch_effects_j[nm] <- rnorm(1, mean = batch_effects[nm], sd = batch_effect_sd)
        }
      } else {
        age_beta_j <- NA_real_
        batch_effects_j <- NULL
      }

      eta <- base_logit +
        group_shift[meta_input$group] +
        group_beta_cos[meta_input$group] * inphase +
        group_beta_sin[meta_input$group] * outphase +
        rand_intercepts[meta_input$id]

      if (include_covariates) {
        eta <- eta +
          age_beta_j * meta_input$age +
          unname(batch_effects_j[meta_input$batch])
      }

      eta_mat[, j] <- eta

      feature_param_store[[j]] <- list(
        base_logit = base_logit,
        group_shift = group_shift,
        group_beta_cos = group_beta_cos,
        group_beta_sin = group_beta_sin,
        age_beta_j = age_beta_j,
        batch_effects_j = batch_effects_j
      )

      truth_row <- list(
        feature = feature_names[j],
        data_type = "rrbs",
        copied_from = NA,
        copied_beta = NA,
        base_logit_intercept = base_logit,
        age_beta = age_beta_j,
        subject_random_intercept_sd = subject_random_intercept_sd,
        coverage_min = coverage_range[1],
        coverage_max = coverage_range[2]
      )

      if (include_covariates) {
        for (b in batch_levels) {
          truth_row[[paste0("batch_effect_", safe_name(b))]] <- unname(batch_effects_j[b])
        }
      }

      for (g in groups) {
        gp <- safe_name(g)
        intercept_g <- base_logit + group_shift[g]
        beta_cos_g <- group_beta_cos[g]
        beta_sin_g <- group_beta_sin[g]
        amp_phase <- calc_amp_phase(beta_cos_g, beta_sin_g)

        truth_row[[paste0("logit_intercept_", gp)]] <- intercept_g
        truth_row[[paste0("group_shift_", gp)]] <- group_shift[g]
        truth_row[[paste0("beta_cos_", gp)]] <- beta_cos_g
        truth_row[[paste0("beta_sin_", gp)]] <- beta_sin_g
        truth_row[[paste0("amplitude_logit_", gp)]] <- amp_phase$amplitude
        truth_row[[paste0("phase_", gp)]] <- amp_phase$phase
        truth_row[[paste0("baseline_prob_", gp)]] <- plogis(intercept_g)
      }

      truth_list[[j]] <- truth_row
    }

    # First generate all observed counts
    M_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
    U_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
    colnames(M_mat) <- feature_names
    colnames(U_mat) <- feature_names

    for (j in seq_len(n_features)) {
      p <- plogis(eta_mat[, j])
      coverage <- sample(seq.int(coverage_range[1], coverage_range[2]), n_samples, replace = TRUE)
      M_mat[, j] <- rbinom(n_samples, size = coverage, prob = p)
      U_mat[, j] <- coverage - M_mat[, j]
    }

    # Then overwrite copied targets using observed counts
    rrbs_copied <- apply_copy_lag_rrbs_counts(
      M_mat = M_mat,
      U_mat = U_mat,
      meta_input = meta_input,
      groups = groups,
      cross_specs = cross_specs
    )
    M_mat <- rrbs_copied$M
    U_mat <- rrbs_copied$U

    for (j in seq_len(n_features)) {
      idx_M <- seq(1, 2 * n_samples, by = 2)
      idx_U <- seq(2, 2 * n_samples, by = 2)

      data_input[idx_M, j] <- M_mat[, j]
      data_input[idx_U, j] <- U_mat[, j]

      obs_weights[idx_M, j] <- M_mat[, j] + U_mat[, j]
      obs_weights[idx_U, j] <- M_mat[, j] + U_mat[, j]
    }

    if (length(cross_specs) > 0) {
      for (spec_i in seq_along(cross_specs)) {
        from_j <- cross_specs[[spec_i]]$from
        to_j   <- cross_specs[[spec_i]]$to
        beta_i <- cross_specs[[spec_i]]$beta
        lag_i  <- cross_specs[[spec_i]]$lag_by_group
        from_nm <- feature_names[from_j]

        truth_list[[to_j]]$copied_from <- from_nm
        truth_list[[to_j]]$copied_beta <- beta_i

        source_params <- feature_param_store[[from_j]]

        for (g in groups) {
          gp <- safe_name(g)
          unshifted_cos <- beta_i * source_params$group_beta_cos[g]
          unshifted_sin <- beta_i * source_params$group_beta_sin[g]
          phi_new <- shift_phase(
            calc_phase_mixedcirc(unshifted_cos, unshifted_sin),
            lag_i[g]
          )
          amp_new <- abs(beta_i) * sqrt(source_params$group_beta_cos[g]^2 + source_params$group_beta_sin[g]^2)
          beta_cos_new <- amp_new * cos(phi_new)
          beta_sin_new <- -amp_new * sin(phi_new)

          truth_list[[to_j]][[paste0("copied_lag_", gp)]] <- lag_i[g]
          truth_list[[to_j]][[paste0("logit_intercept_", gp)]] <- beta_i * (source_params$base_logit + source_params$group_shift[g])
          truth_list[[to_j]][[paste0("group_shift_", gp)]] <- beta_i * source_params$group_shift[g]
          truth_list[[to_j]][[paste0("beta_cos_", gp)]] <- beta_cos_new
          truth_list[[to_j]][[paste0("beta_sin_", gp)]] <- beta_sin_new
          truth_list[[to_j]][[paste0("amplitude_logit_", gp)]] <- amp_new
          truth_list[[to_j]][[paste0("phase_", gp)]] <- phi_new
          truth_list[[to_j]][[paste0("baseline_prob_", gp)]] <- plogis(
            beta_i * (source_params$base_logit + source_params$group_shift[g])
          )
        }

        if (include_covariates) {
          truth_list[[to_j]]$age_beta <- beta_i * source_params$age_beta_j
          for (b in batch_levels) {
            truth_list[[to_j]][[paste0("batch_effect_", safe_name(b))]] <-
              beta_i * unname(source_params$batch_effects_j[b])
          }
        } else {
          truth_list[[to_j]]$age_beta <- NA_real_
        }
      }
    }

    truth <- bind_rows_fill(truth_list)

    return(list(
      data_input = as.data.frame(data_input),
      meta_input = meta_input,
      obs_weights = as.data.frame(obs_weights),
      truth = truth
    ))
  }
}
