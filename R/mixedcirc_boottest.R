#' Bootstrap Testing for Circadian Data
#'
#' This function performs bootstrap testing on fitted (mixed) circadian models.
#'
#' @param object A fitted model object of class \code{mixedcirc_fit}.
#' @param nsim Number of bootstrap simulations.
#' @param seed seed for reproducibitly
#' @param parallel The type of parallel processing to use (\code{"no"}, \code{"multicore"}, \code{"snow"}).
#' @param ncpus Number of CPU cores to use for parallel processing.
#' @param boot.type Type of bootstrap confidence interval (\code{"perc"}, \code{"basic"}, \code{"norm"}).
#' @param conf_level Confidence level for the bootstrap intervals.
#' @param abs_phase Logical; if \code{TRUE}, absolute phases are returned.
#' @param ... Additional arguments passed to \code{bootMer}.
#'
#' @return A \code{mixedcirc_boot} object containing:
#'   \itemize{
#'     \item \code{results}: A data frame with bootstrap statistics.
#'     \item \code{metadata}: Metadata about the bootstrap procedure.
#'     \item \code{original_model}: The original fitted model object.
#'   }
#' @examples
#' \dontrun{
#'   # Example usage
#'   result <- mixedcirc_boottest(fitted_model, nsim = 1000)
#' }
#' @export
mixedcirc_boottest<-function(object,nsim,seed,parallel = c("no", "multicore", "snow"),
                             ncpus = getOption("boot.ncpus", 1L),
                             boot.type=c("perc","basic","norm")[1],conf_level=0.95,abs_phase=TRUE,...){

  if (!"mixedcirc_fit" %in% class(object)) {
    stop("Object must be mixedcirc_fit!")
  }

  fit <- object@fit

  model_type<-NA
  if (class(fit) == "lm") {
    model_type <- "lm"
  }
  else {
    model_type <- "lme4"
  }

  grouped<-length(unique(na.omit(object@exp_design$group)))>1
  group_id <- base::levels(object@exp_design$group)
  bk_transform<-function(coeffs,grouped){
    if(grouped){
      if (any(base::grepl(paste0(group_id[1], ":"),
                          colnames(coeffs)))) {
        rhy_params <- coeffs[, base::paste("group",group_id[1],":",
                                           c(ifelse(object@type=="RRBS", yes = "scaler:inphase", no = "inphase"),
                                             ifelse(object@type=="RRBS", yes = "scaler:outphase", no = "outphase")),
                                           sep = ""),drop=F]
        mesor_A <- coeffs[, 1]

        amps_A <- sqrt(base::rowSums(rhy_params^2))
        sb <- sign(rhy_params[1])
        sg <- sign(rhy_params[2])
        theta <- atan(abs(rhy_params[, 2]/rhy_params[,
                                                     1]))
        if ((sb == 1 | sb == 0) & sg == 1) {
          phi <- -theta
        }else if (sb == -1 & (sg == 1 | sg == 0)) {
          phi <- theta - pi
        }else if ((sb == -1 | sb == 0) & sg == -1) {
          phi <- -theta - pi
        }else if (sb == 1 & (sg == -1 | sg == 0)) {
          phi <- theta - (2 * pi)
        }
        phases_A <- phi

      }
      else {
        amps_A <- 0
        phases_A <- 0
      }
      if (any(base::grepl(paste0(group_id[2], ":"),
                          colnames(coeffs)))) {
        rhy_params <- coeffs[, base::paste("group",group_id[2],":",
                                           c(ifelse(object@type=="RRBS", yes = "scaler:inphase", no = "inphase"),
                                             ifelse(object@type=="RRBS", yes = "scaler:outphase", no = "outphase")),
                                           sep = ""),drop=F]
        mesor_B <- coeffs[, 2]

        amps_B <- sqrt(base::rowSums(rhy_params^2))
        sb <- sign(rhy_params[1])
        sg <- sign(rhy_params[2])
        theta <- atan(abs(rhy_params[, 2]/rhy_params[,
                                                     1]))
        if ((sb == 1 | sb == 0) & sg == 1) {
          phi <- -theta
        }
        else if (sb == -1 & (sg == 1 | sg == 0)) {
          phi <- theta - pi
        }
        else if ((sb == -1 | sb == 0) & sg == -1) {
          phi <- -theta - pi
        }
        else if (sb == 1 & (sg == -1 | sg == 0)) {
          phi <- theta - (2 * pi)
        }
        phases_B <- phi
      }
      else {
        amps_B <- 0
        phases_B <- 0
      }
      if (abs_phase) {
        dt_out <- data.frame(
          mesor_A = abs(mesor_A),
          mesor_B = abs(mesor_B),
          mesor_Avsmesor_B = mesor_A - mesor_B,
          amps_A = abs(amps_A),
          amps_B = abs(amps_B),
          amps_Avsamps_B = amps_A - amps_B,
          phases_A = abs(phases_A),
          phases_B = abs(phases_B),
          phases_Avsphases_B = phases_A - phases_B
        )
      } else {
        dt_out <- data.frame(
          mesor_A = mesor_A,
          mesor_B = mesor_B,
          mesor_Avsmesor_B = mesor_A - mesor_B,
          amps_A = amps_A,
          amps_B = amps_B,
          amps_Avsamps_B = amps_A - amps_B,
          phases_A = phases_A,
          phases_B = phases_B,
          phases_Avsphases_B = phases_A - phases_B
        )
      }

      colnames(dt_out) <- gsub(pattern = "_A", replacement = paste("_",
                                                                   group_id[1], sep = ""), x = colnames(dt_out),
                               fixed = T)
      colnames(dt_out) <- gsub(pattern = "_B", replacement = paste("_",
                                                                   group_id[2], sep = ""), x = colnames(dt_out),
                               fixed = T)
      dt_out
    }else{
      rhy_params <- coeffs[, c(ifelse(object@type=="RRBS", yes = "scaler:inphase",
                                      no = "inphase"),
                               ifelse(object@type=="RRBS", yes = "scaler:outphase",
                                      no = "outphase")),drop=F]



      mesor <- coeffs[1]
      amps_A <- sqrt(base::rowSums(rhy_params^2))
      sb <- sign(rhy_params[1])
      sg <- sign(rhy_params[2])
      theta <- atan(abs(rhy_params[, 2]/rhy_params[,
                                                   1]))
      if ((sb == 1 | sb == 0) & sg == 1) {
        phi <- -theta
      }
      else if (sb == -1 & (sg == 1 | sg == 0)) {
        phi <- theta - pi
      }
      else if ((sb == -1 | sb == 0) & sg == -1) {
        phi <- -theta - pi
      }
      else if (sb == 1 & (sg == -1 | sg == 0)) {
        phi <- theta - (2 * pi)
      }
      phases_A <- phi
      dt_out <- data.frame(amps = amps_A, phases = phases_A)
      if (abs_phase) {
        dt_out <- abs(dt_out)
      }
      dt_out
    }

  }


  coeff_fun<-NULL

  if(model_type=="lm"){
    coeff_fun<-function(x){
      unlist(bk_transform(data.frame(t(coef(x)),check.names = F),grouped = grouped))
    }
  }else if(model_type=="lme4"){
    coeff_fun<-function(x){
      unlist(bk_transform(data.frame(t(fixef(x)),check.names = F),grouped = grouped))
    }
  }


  boot_results<-boot_lm(fit,coeff_fun,nsim = nsim,seed = seed,parallel = parallel, ncpus = ncpus,...)

  p_values<-sapply(names(boot_results$t0),function(nm){
    mean(abs(boot_results$t[,nm] - mean(boot_results$t[,nm]))>abs(boot_results$t0[nm]))
  })

  names(p_values)<-paste(names(p_values),".p_value",sep="")

  results<-cbind(confint(boot_results, type = boot.type, level = conf_level),
        boot.estimate = boot_results$t0,
        boot.SE = apply(boot_results$t, 2, sd),
        pvalue = p_values)

  return(mixedcirc_boot(
    results = data.frame(results,check.names = F),
    metadata = list(
      nsim = nsim,
      seed=seed,
      parallel = parallel,
      conf_level = conf_level,
      boot.type = boot.type
    ),
    original_model = object
  ))
}

refit.lm <- function(object, new_response) {
  # Extract the original model frame
  mf <- model.frame(object)

  # Replace the response with the new response
  mf[[1]] <- new_response

  # Extract the weights from the original model
  wts <- weights(object)

  # Refit the model with the updated response and stored weights
  new_object <- lm(formula = formula(object), data = mf, weights = wts)

  return(new_object)
}


initialize.parallel <- expression({
  have_mc <- have_snow <- FALSE
  if (length(parallel)>1) parallel <- match.arg(parallel)
  do_parallel <- (parallel != "no" && ncpus > 1L)
  if (do_parallel) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!(have_mc || have_snow))
      do_parallel <- FALSE # (only for "windows")
  }
})

boot_lm<-
  function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
            type = c("parametric", "semiparametric"), verbose = FALSE,
            .progress = "none", PBargs = list(), parallel = c("no", "multicore",
                                                              "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL)
  {
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    if (.progress != "none") {
      pbfun <- get(paste0(.progress, "ProgressBar"))
      setpbfun <- get(paste0("set", .simpleCap(.progress),
                             "ProgressBar"))
      pb <- do.call(pbfun, PBargs)
    }
    do_parallel <- have_mc <- have_snow <- NULL
    eval(initialize.parallel)
    if (do_parallel && .progress != "none")
      message("progress bar disabled for parallel operations")
    FUN <- match.fun(FUN)
    type <- match.arg(type)
    if (!is.null(seed))
      set.seed(seed)
    else if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
      stop("This function currently only handles functions that return numeric vectors")

    if(is(x,"lm"))
    {
      if (type == "parametric") {
        argList <- list(x, nsim = nsim, na.action = na.exclude)
        if (!missing(re.form)) {
          argList <- c(argList)
        }
        else {
          argList <- c(argList)
        }
        ss <- do.call(simulate, argList)
      }
      else {
        stop("semiparametric bootstrapping with lm not yet implemented")
      }
      control <- if (!is(x, "lm"))
        NULL
      else eval.parent(x$call$control)
      ffun <- local({
        FUN
        refit
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- lme4:::factory(function(i) FUN(refit(x, ss[[i]])),
                             errval = rep(NA, length.t0))
        function(i) {
          ret <- f1(i)
          if (verbose) {
            cat(sprintf("%5d :", i))
            str(ret)
          }
          if (!do_parallel && .progress != "none") {
            setpbfun(pb, i/nsim)
          }
          ret
        }
      })
      simvec <- seq_len(nsim)
      res <- if (do_parallel) {
        if (have_mc) {
          parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        }
        else if (have_snow) {
          if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(rep("localhost",
                                                 ncpus))
            parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
              parallel::clusterSetRNGStream(cl)
            res <- parallel::parLapply(cl, simvec, ffun)
            parallel::stopCluster(cl)
            res
          }
          else parallel::parLapply(cl, simvec, ffun)
        }
      }
      else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning",
                                         "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      if ((numFail <- sum(msgs[["factory-error"]])) > 0) {
        warning("some bootstrap runs failed (", numFail, "/",
                nsim, ")")
      }
      fail.msgs <- if (numFail == 0)
        NULL
      else msgs[["factory-error"]]
      s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = model.frame(x),
                          seed = .Random.seed, statistic = FUN, sim = "parametric",
                          call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = FALSE),
                     class = c("bootMer", "boot"))
      attr(s, "bootFail") <- numFail
      attr(s, "boot.fail.msgs") <- fail.msgs
      attr(s, "boot.all.msgs") <- msgs
      attr(s, "boot_type") <- "boot"
      s
    }else{
      mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
      if (isLMM(x))
        mle <- c(mle, list(sigma = sigma(x)))
      if (type == "parametric") {
        argList <- list(x, nsim = nsim, na.action = na.exclude)
        if (!missing(re.form)) {
          argList <- c(argList, list(re.form = re.form))
        }
        else {
          argList <- c(argList, list(use.u = use.u))
        }
        ss <- do.call(simulate, argList)
      }
      else {
        if (!missing(re.form))
          stop(paste(sQuote("re.form")), "cannot be used with semiparametric bootstrapping")
        if (use.u) {
          if (isGLMM(x))
            warning("semiparametric bootstrapping is questionable for GLMMs")
          ss <- replicate(nsim, fitted(x) + sample(residuals(x,
                                                             "response"), replace = TRUE), simplify = FALSE)
        }
        else {
          stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
        }
      }
      control <- if (!is(x, "merMod"))
        NULL
      else eval.parent(x@call$control)
      ffun <- local({
        FUN
        refit
        x
        ss
        verbose
        do_parallel
        control
        length.t0 <- length(t0)
        f1 <- lme4:::factory(function(i) FUN(refit(x, ss[[i]], control = control)),
                             errval = rep(NA, length.t0))
        function(i) {
          ret <- f1(i)
          if (verbose) {
            cat(sprintf("%5d :", i))
            str(ret)
          }
          if (!do_parallel && .progress != "none") {
            setpbfun(pb, i/nsim)
          }
          ret
        }
      })
      simvec <- seq_len(nsim)
      res <- if (do_parallel) {
        if (have_mc) {
          parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        }
        else if (have_snow) {
          if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(rep("localhost",
                                                 ncpus))
            parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
              parallel::clusterSetRNGStream(cl)
            res <- parallel::parLapply(cl, simvec, ffun)
            parallel::stopCluster(cl)
            res
          }
          else parallel::parLapply(cl, simvec, ffun)
        }
      }
      else lapply(simvec, ffun)
      t.star <- do.call(cbind, res)
      rownames(t.star) <- names(t0)
      msgs <- list()
      for (mtype in paste0("factory-", c("message", "warning",
                                         "error"))) {
        msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
        msgs[[mtype]] <- table(msgs[[mtype]])
      }
      if ((numFail <- sum(msgs[["factory-error"]])) > 0) {
        warning("some bootstrap runs failed (", numFail, "/",
                nsim, ")")
      }
      fail.msgs <- if (numFail == 0)
        NULL
      else msgs[["factory-error"]]
      s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = model.frame(x),
                          seed = .Random.seed, statistic = FUN, sim = "parametric",
                          call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                     class = c("bootMer", "boot"))
      attr(s, "bootFail") <- numFail
      attr(s, "boot.fail.msgs") <- fail.msgs
      attr(s, "boot.all.msgs") <- msgs
      attr(s, "boot_type") <- "boot"
      s
    }

  }
