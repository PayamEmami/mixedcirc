#' Perform circadian and differential circadian analysis using linear or mixed-effects models
#'
#' Fits one model per variable to estimate mesor, amplitude, phase, within-group rhythmicity,
#' and, when group information is provided, between-group differences in rhythmicity.
#' The function supports ordinary linear regression (`lm`) and linear mixed-effects regression
#' (`lmer` via `lme4`). It can also analyze RRBS-style methylation data using a paired
#' methylated / unmethylated representation.
#'
#' @param data_input A numeric matrix, data.frame, or `DGEList` containing the input data.
#'   Rows correspond to observations and columns correspond to variables when `data_input`
#'   is a matrix or data.frame. If a `DGEList` is provided, rows correspond to variables and
#'   columns correspond to samples, following the standard `limma` / `edgeR` convention.
#'   When a `DGEList` is supplied, expression values and observational weights are estimated
#'   internally using `voom()` or `voomWithDreamWeights()`.
#' @param meta_input A data.frame or matrix containing the sample-level metadata used to build
#'   the design matrix. This must include the columns referenced by `time`, `group`, `id`,
#'   `replicate_id`, and any terms in `formula_extra`.
#' @param time A character string giving the name of the column in `meta_input` that contains
#'   circadian time.
#' @param group Optional character string giving the name of the grouping column in `meta_input`.
#'   If provided, group-specific mesor, amplitude, phase, and rhythmicity p-values are returned,
#'   together with between-group rhythm-difference tests and a global rhythmicity test across groups.
#'   If `NULL`, all samples are treated as belonging to a single group.
#' @param id Optional character string giving the name of the subject identifier column in
#'   `meta_input`. This is used for random effects when `lm_method = "lme"`.
#'   If `NULL`, the function switches to ordinary linear regression (`lm`).
#' @param period A single numeric value specifying the assumed period of the rhythm.
#'   Default is `24`.
#' @param lm_method Regression framework to use. Must be one of `"lm"` or `"lme"`.
#'   `"lm"` fits an ordinary linear model. `"lme"` fits a mixed-effects model using
#'   `lme4::lmer()`. Default is `"lme"`.
#' @param formula_extra Optional formula containing additional fixed-effect covariates to include
#'   in the model. The variables in this formula must be column names in `meta_input`.
#'   For RRBS data, these extra terms are automatically interacted with `scaler`.
#' @param random_effect_type Type of random-effect structure to use when `lm_method = "lme"`.
#'   Must be one of `"none"`, `"mesor"`, or `"mesor_phase"`.
#'   `"none"` adds no random effect,
#'   `"mesor"` adds a random intercept `(1 | id)`,
#'   and `"mesor_phase"` adds random intercept and rhythm terms.
#'   For RRBS this becomes `(1 + scaler:inphase + scaler:outphase | id)`;
#'   otherwise it becomes `(1 + inphase + outphase | id)`.
#' @param lmer.df Method used by `emmeans` / `emtrends` for denominator degrees of freedom when
#'   `lm_method = "lme"`. Must be one of `"Satterthwaite"` or `"Kenward-Roger"`.
#'   Default is `"Satterthwaite"`.
#' @param obs_weights Optional matrix of observational weights with the same dimensions as
#'   `data_input` when `data_input` is a matrix or data.frame. Each column corresponds to one
#'   variable and each row to one observation. For RRBS input, `obs_weights` is required unless
#'   weights are estimated by setting `force_weight_estimation = TRUE`.
#' @param RRBS Logical. If `TRUE`, the data are assumed to represent RRBS / BS-seq style
#'   methylation data, with each biological observation represented by two consecutive rows:
#'   methylated count followed by unmethylated count. In this mode, `obs_weights` must be supplied
#'   unless weights are estimated internally, and `replicate_id` must be provided.
#'   Default is `FALSE`.
#' @param replicate_id Character string giving the name of the column in `meta_input` that
#'   identifies unique technical / replicate observations for RRBS data. This is required when
#'   `RRBS = TRUE`.
#' @param force_weight_estimation Logical. If `TRUE`, weights are estimated internally using
#'   `voom()` or `voomWithDreamWeights()` even when `data_input` is not a `DGEList`.
#'   Default is `FALSE`.
#' @param ncores Number of cores to use for parallel processing. Default is `1`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default is `FALSE`.
#' @param store_models Logical. If `TRUE`, retain each fitted model for downstream
#'   prediction, bootstrap and diagnostic methods. Set to `FALSE` for a much
#'   smaller result object when only the tabular results are needed.
#' @param on_error How feature-specific fitting errors are handled. `"stop"`
#'   preserves fail-fast behavior. `"keep"` returns a failed result with the
#'   error recorded in [mixedcirc_results()] diagnostics and continues fitting
#'   the remaining features.
#' @param ... Additional arguments passed to the underlying model-fitting function.
#'
#' @return
#' An object of class `mixedcirc_fit_list`. Each element corresponds to one analyzed variable
#' and contains:
#' \describe{
#'   \item{results}{A one-row data.frame with estimated mesor, amplitude, phase, and p-values.}
#'   \item{fit}{The fitted model object (`lm` or `lmerMod`).}
#'   \item{exp_design}{The experimental design used for fitting, including generated `inphase`
#'   and `outphase` columns.}
#'   \item{type}{A character string indicating the analysis type, e.g. `"expression"` or `"RRBS"`.}
#' }
#'
#' @details
#' Circadian time is encoded using
#' \deqn{inphase = cos(2\pi t / period)}
#' and
#' \deqn{outphase = sin(2\pi t / period).}
#'
#' For non-RRBS data with group information, the default fixed-effects model is
#' \preformatted{
#' measure ~ 1 + group + group:inphase + group:outphase
#' }
#'
#' For non-RRBS data without group information, the fixed-effects model is
#' \preformatted{
#' measure ~ 1 + inphase + outphase
#' }
#'
#' For RRBS data with group information, the default fixed-effects model is
#' \preformatted{
#' measure ~ 0 + replicate_id + scaler + group:scaler + group:inphase:scaler + group:outphase:scaler
#' }
#'
#' For RRBS data without group information, the fixed-effects model is
#' \preformatted{
#' measure ~ 0 + replicate_id + scaler + inphase:scaler + outphase:scaler
#' }
#'
#' Additional fixed effects supplied through `formula_extra` are appended to the model.
#'
#' Rhythmicity within each group is tested jointly from the estimated `inphase` and `outphase`
#' coefficients using `emmeans` / `emtrends`. In grouped analyses, the function also computes
#' pairwise between-group tests for differences in rhythm coefficients, as well as a global
#' joint rhythmicity test across all groups.
#'
#' Amplitude is computed as
#' \deqn{\sqrt{\beta_{in}^2 + \beta_{out}^2}}
#' and phase is derived from the estimated sine and cosine coefficients using `calc_phase()`.
#' The returned phase is the signed shift \eqn{\phi} in
#' \deqn{M + A cos(2\pi t / period + \phi).}
#' Consequently, peak time is obtained from \eqn{-\phi}; use
#' `mixedcirc_detectpeaks()` to perform this conversion and wrap the result into
#' the observed cycle.
#'
#' For RRBS input, the rows of `data_input` are assumed to alternate between methylated and
#' unmethylated counts for each biological observation. For example, if the biological samples are
#' `A1_1`, `A1_2`, `A2_1`, and `B1_1`, then the rows of `data_input` must be ordered as
#' \preformatted{
#' A1_1_methylated, A1_1_unmethylated,
#' A1_2_methylated, A1_2_unmethylated,
#' A2_1_methylated, A2_1_unmethylated,
#' B1_1_methylated, B1_1_unmethylated
#' }
#' In this case, `meta_input` should contain one row per biological observation
#' (`A1_1`, `A1_2`, `A2_1`, `B1_1`), and the function internally expands it to match the
#' methylated / unmethylated row structure.
#'
#' If weights are not supplied directly, they may be estimated from count-like matrix data using
#' `limma::voom()` for ordinary linear models or `voomWithDreamWeights()`
#' for mixed-effects models.
#'
#' No multiple-testing adjustment is performed inside this function.
#'
#' @examples
#' data("circa_data")
#'
#' res <- mixedcirc_detect(
#'   data_input = circa_data$data_matrix,
#'   meta_input = data.frame(
#'     time = circa_data$time,
#'     group = circa_data$group,
#'     id = circa_data$id
#'   ),
#'   time = "time",
#'   group = "group",
#'   id = "id",
#'   period = 24,
#'   lm_method = "lm",
#'   verbose = TRUE
#' )
#'
#' @export

mixedcirc_detect <- function(
    data_input = NULL,
    meta_input = NULL,
    time = NULL,
    group = NULL,
    id = NULL,
    period = 24,
    lm_method = c("lm", "lme")[2],
    formula_extra = NULL,
    random_effect_type = c("none", "mesor", "mesor_phase")[2],
    lmer.df = c("Satterthwaite", "Kenward-Roger")[1],
    obs_weights = NULL,
    RRBS = FALSE,
    replicate_id = NULL,
    force_weight_estimation = FALSE,
    ncores = 1,
    verbose = FALSE,
    store_models = TRUE,
    on_error = c("stop", "keep"),
    ...
) {
  call_in<-match.call()
  fit_arguments <- list(...)
  protected_fit_arguments <- intersect(names(fit_arguments), c("formula", "data", "weights"))
  if (length(protected_fit_arguments)) {
    stop(
      "Do not supply ", paste(protected_fit_arguments, collapse = ", "),
      " through ...; mixedcirc constructs these arguments",
      call. = FALSE
    )
  }
  previous_plan <- future::plan("list")
  on.exit(future::plan(previous_plan), add = TRUE)

  if (isTRUE(verbose)) cat("Checking inputs ...\n")

  # -----------------------------
  # Input checks
  # -----------------------------
  if (is.null(data_input)) {
    stop("data_input must be a data frame, matrix, or DGEList")
  }

  if (is.null(meta_input)) {
    stop("meta_input must be a data frame or matrix")
  }

  if (!is.matrix(data_input) && !is.data.frame(data_input) && !is(data_input, "DGEList")) {
    stop("data_input must be a data frame, matrix, or DGEList")
  }

  if (!is.matrix(meta_input) && !is.data.frame(meta_input)) {
    stop("meta_input must be a data frame or matrix")
  }

  meta_input <- as.data.frame(meta_input)
  if (anyDuplicated(names(meta_input))) {
    stop("meta_input must have unique column names")
  }
  reserved_metadata <- intersect(
    names(meta_input),
    c("measure", "inphase", "outphase", "scaler", ".used_in_fit",
      ".fitted", ".residual", ".weight", ".mixedcirc_weight")
  )
  if (length(reserved_metadata) > 0L) {
    stop(
      "meta_input uses column names reserved by mixedcirc: ",
      paste(reserved_metadata, collapse = ", ")
    )
  }

  if (!is.character(time) || length(time) != 1L || is.na(time) || time == "") {
    stop("time must be a column name in meta_input")
  }

  if (!time %in% colnames(meta_input)) {
    stop("time must be a column name in meta_input")
  }
  if (make.names(time) != time) {
    stop("time must refer to a syntactically valid R column name")
  }

  if (!is.numeric(meta_input[[time]]) || all(!is.finite(meta_input[[time]])) ||
      any(is.infinite(meta_input[[time]]))) {
    stop("the time column must be numeric, contain a finite value, and not contain Inf or -Inf")
  }

  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores %% 1 != 0) {
    stop("ncores must be a positive integer")
  }
  ncores <- as.integer(ncores)

  logical_args <- list(
    RRBS = RRBS,
    force_weight_estimation = force_weight_estimation,
    verbose = verbose,
    store_models = store_models
  )
  bad_logical <- !vapply(
    logical_args,
    function(z) is.logical(z) && length(z) == 1L && !is.na(z),
    logical(1)
  )
  if (any(bad_logical)) {
    stop(names(logical_args)[which(bad_logical)[1L]], " must be TRUE or FALSE")
  }

  on_error <- match.arg(on_error)

  if (is.data.frame(data_input)) {
    numeric_columns <- vapply(data_input, is.numeric, logical(1))
    if (any(!numeric_columns)) {
      stop("all columns of data_input must be numeric")
    }
  }
  if (!is(data_input, "DGEList") &&
      (!is.numeric(as.matrix(data_input)) || any(is.infinite(as.matrix(data_input))))) {
    stop("data_input must contain numeric values without Inf or -Inf")
  }

  if (RRBS && is(data_input, "DGEList")) {
    stop("RRBS = TRUE is supported for paired matrix/data.frame input, not DGEList input")
  }

  if (!is(data_input, "DGEList")) {
    if(RRBS){
      if (nrow(meta_input)*2 != nrow(data_input)) {
        stop(
          "For RRBS input, meta_input must have one row per biological sample ",
          "and data_input must have two rows per sample"
        )
      }
    }else{
      if (nrow(meta_input) != nrow(data_input)) {
        stop("The number of rows of meta_input is not equal to the number of rows of data_input")
      }
    }
  } else {
    if (nrow(meta_input) != ncol(data_input)) {
      stop("The number of rows of meta_input is not equal to the number of columns of data_input")
    }
  }

  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || is.na(group) || group == "") {
      stop("group must be NULL or a single column name in meta_input")
    }
    if (!group %in% colnames(meta_input)) {
      stop("group must be the name of a column in meta_input")
    }
    if (make.names(group) != group) {
      stop("group must refer to a syntactically valid R column name")
    }
    observed_groups <- unique(meta_input[[group]][!is.na(meta_input[[group]])])
    if (length(observed_groups) == 0L) {
      stop("the supplied group column has no observed values")
    }
    multiple_groups <- length(observed_groups) > 1L
    if (!multiple_groups) {
      warning("The supplied group column has fewer than two observed groups; fitting a single-group model")
    }
  } else {
    warning("No group information was provided! Skipping differential analysis!")
    multiple_groups <- FALSE
    meta_input$dummy <- "all"
    group <- "dummy"
  }

  if (!is.null(id)) {
    if (!is.character(id) || length(id) != 1L || is.na(id) || id == "") {
      stop("id must be NULL or a single column name in meta_input")
    }
    if (!id %in% colnames(meta_input)) {
      stop("id must be the name of a column in meta_input")
    }
    if (make.names(id) != id) {
      stop("id must refer to a syntactically valid R column name")
    }
  } else {
    warning("No id information was provided! Switching to normal linear regression!")
    lm_method <- "lm"
    meta_input$dummy_id <- "all"
    id <- "dummy_id"
  }

  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0) {
    stop("period must be a single positive finite numeric value")
  }

  lm_method <- match.arg(lm_method, c("lm", "lme"))
  lmer.df <- match.arg(lmer.df, c("Satterthwaite", "Kenward-Roger"))
  random_effect_type <- match.arg(random_effect_type, c("none", "mesor", "mesor_phase"))

  if (lm_method == "lme" && random_effect_type == "none") {
    warning("random_effect_type = 'none' has no random effect; using lm_method = 'lm'")
    lm_method <- "lm"
  }

  if (lm_method == "lme" && length(unique(meta_input[[id]][!is.na(meta_input[[id]])])) < 2L) {
    stop("mixed-effects fitting requires at least two observed id levels")
  }

  if (RRBS) {
    if (is.null(obs_weights) && !force_weight_estimation) {
      stop("For RRBS data, obs_weights must be set unless force_weight_estimation = TRUE")
    }

    if (!is.character(replicate_id) || length(replicate_id) != 1L ||
        is.na(replicate_id) || replicate_id == "" ||
        !replicate_id %in% colnames(meta_input)) {
      stop("For RRBS data, replicate_id must be set and must be a column name in meta_input")
    }
    if (make.names(replicate_id) != replicate_id) {
      stop("replicate_id must refer to a syntactically valid R column name")
    }
    if (anyNA(meta_input[[replicate_id]]) || anyDuplicated(meta_input[[replicate_id]])) {
      stop("For RRBS data, replicate_id must be non-missing and unique for every biological observation")
    }
  }

  if (!is.null(obs_weights) && !is(data_input, "DGEList")) {
    if (!is.matrix(obs_weights) && !is.data.frame(obs_weights)) {
      stop("obs_weights must be a numeric matrix or data.frame")
    }
    obs_weights <- as.matrix(obs_weights)
    if (!identical(dim(obs_weights), dim(data_input))) {
      stop("obs_weights must be the same size as data_input")
    }
    if (!is.numeric(obs_weights) || any(!is.finite(obs_weights)) || any(obs_weights <= 0)) {
      stop("obs_weights must contain finite positive numeric values")
    }
  }


  if (!is.null(formula_extra)) {
    if (!inherits(formula_extra, "formula")) {
      stop("formula_extra must be a formula")
    }
    if (length(formula_extra) != 2L) {
      stop("formula_extra must be one-sided, for example ~ age + batch")
    }

    missing_formula_extra_terms <- setdiff(all.vars(formula_extra), colnames(meta_input))
    if (length(missing_formula_extra_terms) > 0) {
      stop(
        "formula_extra terms must be column names in meta_input: ",
        paste(missing_formula_extra_terms, collapse = ", ")
      )
    }
  }

  if (!is(data_input, "DGEList")) {
    eset <- as.matrix(data_input[, , drop = FALSE])
  } else {
    eset <- NA
  }

  feature_count <- if (is(data_input, "DGEList")) nrow(data_input) else ncol(data_input)
  if (feature_count < 1L) {
    stop("data_input must contain at least one feature")
  }

  if (ncores > 1L) {
    future::plan(future::multisession, workers = ncores)
  } else {
    future::plan(future::sequential)
  }

  # -----------------------------
  # Build experiment design
  # -----------------------------
  if (verbose) cat("Building experiment ...\n")

  exp_design <- meta_input

  if (RRBS) {
    exp_design <- meta_input[rep(seq_len(nrow(meta_input)), each = 2), , drop = FALSE]
    exp_design$scaler <- rep(c(1, 0), nrow(meta_input))
  }

  exp_design[, group] <- factor(exp_design[, group])

  # -----------------------------
  # Helper for random effects
  # -----------------------------
  build_random_term <- function(random_effect_type, id, RRBS = FALSE) {
    if (random_effect_type == "none") {
      return(NULL)
    }

    if (random_effect_type == "mesor") {
      return(paste0("(1|", id, ")"))
    }

    if (random_effect_type == "mesor_phase") {
      if (RRBS) {
        return(paste0("(1 + scaler:inphase + scaler:outphase|", id, ")"))
      } else {
        return(paste0("(1 + inphase + outphase|", id, ")"))
      }
    }

    stop("Unsupported random_effect_type")
  }

  # -----------------------------
  # Build fixed-effects formula
  # -----------------------------
  if (verbose) cat("Building formula ...\n")

  if (RRBS) {
    if (multiple_groups) {
      formula_input <- as.formula(
        paste0(
          "measure ~ 0 + ", replicate_id,
          " + scaler",
          " + ", group, ":scaler",
          " + ", group, ":inphase:scaler",
          " + ", group, ":outphase:scaler"
        )
      )
    } else {
      formula_input <- as.formula(
        paste0(
          "measure ~ 0 + ", replicate_id,
          " + scaler",
          " + inphase:scaler",
          " + outphase:scaler"
        )
      )
    }
  } else {
    if (multiple_groups) {
      formula_input <- as.formula(
        paste0(
          "measure ~ 1 + ", group,
          " + ", group, ":inphase",
          " + ", group, ":outphase"
        )
      )
    } else {
      formula_input <- as.formula("measure ~ 1 + inphase + outphase")
    }
  }

  formula_final <- formula_input
  # -----------------------------
  # Add extra fixed effects
  # -----------------------------
  if (!is.null(formula_extra)) {
    rhs1 <- paste(deparse(formula_input[[3]]), collapse = "")

    if (RRBS) {
      formula_extra_mod <- update(formula_extra, . ~ scaler:(.))
      rhs2 <- deparse(formula_extra_mod[[3]])
    } else {
      rhs2 <- deparse(formula_extra[[2]])
    }

    formula_final <- as.formula(
      paste(all.vars(formula_input)[1], "~", paste(rhs1, rhs2, sep = " + "))
    )


  } else {
    formula_final <- formula_input
  }


  # -----------------------------
  # Add random effects for lme
  # -----------------------------

  formula_without_random <- formula_final
  if (lm_method == "lme") {
    rand_term <- build_random_term(
      random_effect_type = random_effect_type,
      id = id,
      RRBS = RRBS
    )

    if (!is.null(rand_term)) {
      rhs <- paste(paste(deparse(formula_final[[3]]), collapse = ""), rand_term, sep = " + ")
      formula_final <- as.formula(
        paste(all.vars(formula_final)[1], "~", rhs)
      )
    }
  }

  if (lm_method == "lm" && random_effect_type != "none") {
    warning("random_effect_type was ignored because lm_method = 'lm'")
  }
  exp_design <- base::cbind(exp_design,
                            inphase = cos(2 * pi * exp_design[,time] / period),
                            outphase = sin(2 * pi * exp_design[,time] / period))


  model_matrix <- stats::model.matrix(
    stats::as.formula(paste(
      "~",
      paste(deparse(formula_without_random[[3L]]), collapse = "")
    )),
    exp_design
  )


  if (is(data_input, "DGEList") || isTRUE(force_weight_estimation)) {
    if (verbose) cat("Estimating weights for the input variables ...\n")

    voom_input <- if (is(data_input, "DGEList")) data_input else t(as.matrix(data_input))

    if (lm_method == "lm") {
      if (verbose) cat("Using voom to estimate precision weights ...\n")
      voom_res <- limma::voom(voom_input, model_matrix, plot = FALSE)
    } else {
      drop_lhs <- function(f) {
        if (length(f) == 3L) {
          f[[2]] <- NULL
        }
        f
      }
      formula_rhs <- drop_lhs(formula_final)
      voom_res <- voomWithDreamWeights(
        voom_input,
        formula = formula_rhs,
        data = exp_design,
        BPPARAM = BiocParallel::SnowParam(workers = ncores)
      )
    }
    obs_weights <- t(voom_res$weights)
    eset <- t(voom_res$E)
  }



  feature_names<-colnames(eset)
  if(length(feature_names)==0){
    feature_names<-paste0("feature_", seq_len(ncol(eset)))
  }
  if (anyDuplicated(feature_names)) {
    warning("Duplicated feature names were made unique")
    feature_names <- make.unique(as.character(feature_names))
  }
  chunks <- parallel::splitIndices(ncol(eset), min(ncol(eset),  ncores))
  if (verbose) cat("Splitting data into", length(chunks), "chunks ...\n")

  res<-foreach(chIndx=chunks)%dofuture%
    {

      outputs_fn<-foreach(i=chIndx)%do%{

        feature_name <- feature_names[i]
        fit_warnings <- character()

        tryCatch(withCallingHandlers({

        if (verbose) cat("Processing variable", i, "...\n")

        data_grouped <- cbind(
          measure = as.numeric(eset[, i]),
          exp_design
        )
        data_grouped <- as.data.frame(data_grouped)

        # ensure stable unique row names for mapping back after model fitting
        if (is.null(rownames(data_grouped)) ||
            anyNA(rownames(data_grouped)) ||
            any(rownames(data_grouped) == "") ||
            anyDuplicated(rownames(data_grouped))) {
          rownames(data_grouped) <- paste0("row_", seq_len(nrow(data_grouped)))
        }

        this_formula <- formula_final
        environment(this_formula) <- environment()
        if (verbose) cat("Fitting the model for variable", i, "...\n")

        weights_i <- if (!is.null(obs_weights)) obs_weights[, i] else NULL

        model_arguments <- c(
          list(formula = this_formula, data = data_grouped),
          if (is.null(weights_i)) list() else list(weights = weights_i),
          fit_arguments
        )
        model_in <- switch(
          lm_method,
          lm = do.call(stats::lm, model_arguments),
          lme = do.call(lme4::lmer, model_arguments)
        )

        mf_used <- stats::model.frame(model_in)
        used_rows <- rownames(mf_used)

        data_grouped$.used_in_fit <- FALSE
        data_grouped[used_rows, ".used_in_fit"] <- TRUE

        data_grouped$.fitted <- NA_real_
        data_grouped[used_rows, ".fitted"] <- stats::fitted(model_in)

        data_grouped$.residual <- NA_real_
        data_grouped[used_rows, ".residual"] <- stats::residuals(model_in)

        data_grouped$.weight <- NA_real_
        if (!is.null(weights_i)) {
          data_grouped$.weight <- weights_i
        }
        rh_information<-list()
        if(multiple_groups==FALSE && RRBS==FALSE){
          em_obj<-as.data.frame(emmeans(
            model_in,
            ~ 1,
            at = list(
              inphase = 0,
              outphase = 0
            ),lmer.df = lmer.df
          ))
          rh_information[["mesor"]]<-as.numeric(em_obj[2])
          em_obj_in<-emtrends(
            model_in,
            ~1,
            var="inphase",lmer.df = lmer.df
          )
          em_obj_out<-emtrends(
            model_in,
            ~1,
            var="outphase",lmer.df = lmer.df
          )
          joint_phase<-rbind(inphase=em_obj_in,outphase=em_obj_out)
          rhy_params<-(data.frame(joint_phase))[,2]
          rh_information[["beta_cos"]]<-rhy_params[1]
          rh_information[["beta_sin"]]<-rhy_params[2]
          rh_information[["amps"]]<-sqrt(sum(rhy_params^2))
          rh_information[["phase"]]<-calc_phase(rhy_params[1],rhy_params[2])
          test_res<-test(joint_phase, null = 0, joint = TRUE)
          rh_information[["p.value"]]<-as.numeric(test_res["p.value"])
        }else if(multiple_groups==FALSE && RRBS==TRUE){ ## For RRBS
          em_obj<-as.data.frame(emtrends(
            model_in,
            ~ 1,
            var = "scaler",
            at = list(
              inphase = 0,
              outphase = 0,scaler=1
            ),lmer.df = lmer.df
          ))
          rh_information[["mesor"]]<-as.numeric(em_obj[2])
          em_obj_in<-emtrends(
            model_in,
            ~scaler,
            var="inphase",
            at=list(scaler=1),lmer.df = lmer.df
          )
          em_obj_out<-emtrends(
            model_in,
            ~scaler,
            var="outphase",
            at=list(scaler=1),lmer.df = lmer.df
          )
          joint_phase<-rbind(inphase=em_obj_in,outphase=em_obj_out)
          rhy_params<-(data.frame(joint_phase))[,2]
          rh_information[["beta_cos"]]<-rhy_params[1]
          rh_information[["beta_sin"]]<-rhy_params[2]
          rh_information[["amps"]]<-sqrt(sum(rhy_params^2))
          rh_information[["phase"]]<-calc_phase(rhy_params[1],rhy_params[2])
          test_res<-test(joint_phase, null = 0, joint = TRUE)
          rh_information[["p.value"]]<-as.numeric(test_res["p.value"])
        }else if(multiple_groups==TRUE && RRBS==FALSE){ ## For Multigroup
          formula_for_group<- as.formula(paste0( "~ 1 + ", group))
          em_obj<-data.frame(emmeans(
            model_in,
            formula_for_group,
            at = list(
              inphase = 0,
              outphase = 0
            ),lmer.df = lmer.df
          ))

          em_obj_in<-emtrends(
            model_in,
            formula_for_group,
            var="inphase",lmer.df = lmer.df
          )
          em_obj_out<-emtrends(
            model_in,
            formula_for_group,
            var="outphase",lmer.df = lmer.df
          )
          for(gr in 1:nrow(em_obj))
          {
            group_inf<-paste0("mesor_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-as.numeric(em_obj[gr,2])

            joint_phase<-rbind(inphase=em_obj_in[gr],outphase=em_obj_out[gr])
            rhy_params<-(data.frame(joint_phase))[,2]
            group_inf<-paste0("beta_cos_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-rhy_params[1]
            group_inf<-paste0("beta_sin_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-rhy_params[2]
            group_inf<-paste0("amps_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-sqrt(sum(rhy_params^2))
            group_inf<-paste0("phase_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-calc_phase(rhy_params[1],rhy_params[2])
            test_res<-test(joint_phase, null = 0, joint = TRUE)
            group_inf<-paste0("p.value_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-as.numeric(test_res["p.value"])
          }

          # calculate contrasts
          pairs_in<-pairs(em_obj_in,lmer.df = lmer.df)
          pairs_out<-pairs(em_obj_out,lmer.df = lmer.df)

          pairs_in_df<-data.frame(pairs_in)
          for(gr in 1:nrow(pairs_in_df))
          {
            group_inf<-paste0("p.value_",as.character(pairs_in_df[gr,1]))
            joint_phase<-rbind(inphase=pairs_in[gr],outphase=pairs_out[gr])
            test_res<-test(joint_phase, null = 0, joint = TRUE)
            rh_information[[group_inf]]<-as.numeric(test_res["p.value"])
          }

          joint_phase<-rbind(inphase=em_obj_in,outphase=em_obj_out)
          test_res<-test(joint_phase, null = 0, joint = TRUE)
          rh_information[["p.value_global"]]<-as.numeric(test_res["p.value"])

          differential_in <- emmeans::contrast(
            em_obj_in, method = "trt.vs.ctrl", ref = 1, adjust = "none"
          )
          differential_out <- emmeans::contrast(
            em_obj_out, method = "trt.vs.ctrl", ref = 1, adjust = "none"
          )
          differential_joint <- rbind(
            inphase = differential_in,
            outphase = differential_out
          )
          differential_test <- test(differential_joint, null = 0, joint = TRUE)
          rh_information[["p.value_differential_global"]] <-
            as.numeric(differential_test["p.value"])

        }else if(multiple_groups==TRUE && RRBS==TRUE){ ## For Multigroup and RRBS
          formula_for_group<- as.formula(paste0( "~ 1 + ", group))
          em_obj<-data.frame(emtrends(
            model_in,
            formula_for_group,var="scaler",
            at = list(
              inphase = 0,
              outphase = 0,scaler=1
            ),lmer.df = lmer.df
          ))

          em_obj_in<-emtrends(
            model_in,
            formula_for_group,
            var="inphase",at=list(scaler=1),lmer.df = lmer.df
          )
          em_obj_out<-emtrends(
            model_in,
            formula_for_group,
            var="outphase",at=list(scaler=1),lmer.df = lmer.df
          )
          for(gr in 1:nrow(em_obj))
          {
            group_inf<-paste0("mesor_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-as.numeric(em_obj[gr,2])

            joint_phase<-rbind(inphase=em_obj_in[gr],outphase=em_obj_out[gr])
            rhy_params<-(data.frame(joint_phase))[,2]
            group_inf<-paste0("beta_cos_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-rhy_params[1]
            group_inf<-paste0("beta_sin_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-rhy_params[2]
            group_inf<-paste0("amps_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-sqrt(sum(rhy_params^2))
            group_inf<-paste0("phase_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-calc_phase(rhy_params[1],rhy_params[2])
            test_res<-test(joint_phase, null = 0, joint = TRUE)
            group_inf<-paste0("p.value_",as.character(em_obj[gr,1]))
            rh_information[[group_inf]]<-as.numeric(test_res["p.value"])
          }

          # calculate contrasts
          pairs_in<-pairs(em_obj_in,lmer.df = lmer.df)
          pairs_out<-pairs(em_obj_out,lmer.df = lmer.df)

          pairs_in_df<-data.frame(pairs_in)
          for(gr in 1:nrow(pairs_in_df))
          {
            group_inf<-paste0("p.value_",as.character(pairs_in_df[gr,1]))
            joint_phase<-rbind(inphase=pairs_in[gr],outphase=pairs_out[gr])
            test_res<-test(joint_phase, null = 0, joint = TRUE)
            rh_information[[group_inf]]<-as.numeric(test_res["p.value"])
          }

          joint_phase<-rbind(inphase=em_obj_in,outphase=em_obj_out)
          test_res<-test(joint_phase, null = 0, joint = TRUE)
          rh_information[["p.value_global"]]<-as.numeric(test_res["p.value"])

          differential_in <- emmeans::contrast(
            em_obj_in, method = "trt.vs.ctrl", ref = 1, adjust = "none"
          )
          differential_out <- emmeans::contrast(
            em_obj_out, method = "trt.vs.ctrl", ref = 1, adjust = "none"
          )
          differential_joint <- rbind(
            inphase = differential_in,
            outphase = differential_out
          )
          differential_test <- test(differential_joint, null = 0, joint = TRUE)
          rh_information[["p.value_differential_global"]] <-
            as.numeric(differential_test["p.value"])

        }
        dt_out<-data.frame(rh_information,check.names = FALSE)
        rownames(dt_out)<-feature_name

        fit_converged <- if (inherits(model_in, "merMod")) {
          convergence_messages <- model_in@optinfo$conv$lme4$messages
          optimizer_code <- model_in@optinfo$conv$opt
          is.null(convergence_messages) &&
            (is.null(optimizer_code) || all(as.numeric(optimizer_code) == 0))
        } else if (inherits(model_in, "glm")) {
          isTRUE(model_in$converged)
        } else {
          TRUE
        }
        fit_singular <- if (inherits(model_in, "merMod")) {
          tryCatch(lme4::isSingular(model_in), error = function(e) NA)
        } else {
          NA
        }
        fit_rank_deficient <- if (inherits(model_in, "merMod")) {
          anyNA(lme4::fixef(model_in, add.dropped = TRUE))
        } else if (inherits(model_in, "lm")) {
          model_in$rank < length(stats::coef(model_in))
        } else {
          NA
        }
        safe_fit_statistic <- function(value) {
          value <- tryCatch(suppressWarnings(as.numeric(value)[1L]),
                            error = function(e) NA_real_)
          if (!length(value) || !is.finite(value)) NA_real_ else value
        }

        new("mixedcirc_fit",results = dt_out,
            fit=if (isTRUE(store_models)) model_in else NULL,
            exp_design=data_grouped,
            type=if (RRBS) "RRBS" else "expression",
            call=call_in,params=list(
              multiple_groups=multiple_groups,
              RRBS=RRBS,
              period=period,
              time=time,
              group=if (multiple_groups) group else NULL,
              id=if (id == "dummy_id") NULL else id,
              replicate_id=replicate_id,
              lm_method=lm_method,
              fit_status="ok",
              fit_message=NA_character_,
              fit_warnings=unique(fit_warnings),
              model_class=class(model_in)[1L],
              n_obs=stats::nobs(model_in),
              converged=fit_converged,
              singular=fit_singular,
              rank_deficient=fit_rank_deficient,
              AIC=safe_fit_statistic(stats::AIC(model_in)),
              BIC=safe_fit_statistic(stats::BIC(model_in)),
              logLik=safe_fit_statistic(stats::logLik(model_in)),
              store_models=store_models
            ))

        }, message = function(condition) {
          if (isTRUE(verbose)) {
            return(invisible(NULL))
          }
          invokeRestart("muffleMessage")
        }, warning = function(condition) {
          fit_warnings <<- c(fit_warnings, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }), error = function(condition) {
          if (identical(on_error, "stop")) stop(condition)

          failed_result <- data.frame(p.value = NA_real_, check.names = FALSE)
          rownames(failed_result) <- feature_name
          failed_design <- cbind(
            measure = as.numeric(eset[, i]),
            exp_design
          )
          failed_design <- as.data.frame(failed_design)
          new(
            "mixedcirc_fit",
            results = failed_result,
            fit = NULL,
            exp_design = failed_design,
            type = if (RRBS) "RRBS" else "expression",
            call = call_in,
            params = list(
              multiple_groups = multiple_groups,
              RRBS = RRBS,
              period = period,
              time = time,
              group = if (multiple_groups) group else NULL,
              id = if (id == "dummy_id") NULL else id,
              replicate_id = replicate_id,
              lm_method = lm_method,
              fit_status = "failed",
              fit_message = conditionMessage(condition),
              fit_warnings = unique(fit_warnings),
              model_class = NA_character_,
              n_obs = NA_integer_,
              converged = FALSE,
              singular = NA,
              rank_deficient = NA,
              AIC = NA_real_,
              BIC = NA_real_,
              logLik = NA_real_,
              store_models = store_models
            )
          )
        })

      }


      outputs_fn
    }

  result_list <- unlist(res, recursive = FALSE, use.names = FALSE)
  names(result_list) <- feature_names
  res_tmp<-new("mixedcirc_fit_list",results = result_list)
  return((res_tmp))

}

calc_phase<-function(inphase,outphase){
  if (length(inphase) != 1L || length(outphase) != 1L ||
      !is.finite(inphase) || !is.finite(outphase) ||
      (inphase == 0 && outphase == 0)) {
    return(NA_real_)
  }

  # atan2(-sin coefficient, cos coefficient) is the signed shift in
  # M + A*cos(omega*t + phi). Phases are represented on [-2*pi, 0).
  phi <- atan2(-outphase, inphase)
  if (phi >= 0) phi <- phi - 2 * pi
  phi
}
