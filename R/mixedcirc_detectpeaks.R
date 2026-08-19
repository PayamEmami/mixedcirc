#' Convert a fitted phase to peak time
#'
#' Converts the signed phase returned by `mixedcirc_detect()` to the time at
#' which the fitted rhythm reaches its maximum.
#'
#' @param phi Numeric vector of signed phases returned by `mixedcirc_detect()`.
#' @param period Positive rhythm period. Default is `24`.
#' @param min_time Minimum time considered when `correct_min = TRUE`.
#' @param correct_period Logical; if `TRUE`, peak times are wrapped into
#'   `[0, period)`.
#' @param correct_min Logical; if `TRUE`, wrapped peak times below `min_time`
#'   are advanced by complete periods until they are at least `min_time`.
#' @export
#' @examples
#' data("circa_data")
#'
#' results <- mixedcirc_detect(
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
#' phase_columns <- grep("^phase_", names(results[[1]]@results), value = TRUE)
#' phases <- unlist(results[[1]]@results[phase_columns], use.names = FALSE)
#' mixedcirc_detectpeaks(phases, min_time = 0)
#'
#' @return A numeric vector of peak times.
#'
#' @details
#' `mixedcirc_detect()` reports the signed phase \eqn{\phi} from
#' \deqn{M + A cos(2\pi t / period + \phi).}
#' The corresponding peak time is
#' \deqn{(-\phi mod 2\pi) period / (2\pi).}
#'
#'
#'


mixedcirc_detectpeaks <- function(phi,
                                  period = 24,
                                  min_time = 0,
                                  correct_period = TRUE,
                                  correct_min = FALSE) {
  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(period) || length(period) != 1L ||
      !is.finite(period) || period <= 0) {
    stop("`period` must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(min_time) || length(min_time) != 1L || !is.finite(min_time)) {
    stop("`min_time` must be a single finite number.", call. = FALSE)
  }

  if (!is.logical(correct_period) || length(correct_period) != 1L ||
      is.na(correct_period)) {
    stop("`correct_period` must be either TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(correct_min) || length(correct_min) != 1L ||
      is.na(correct_min)) {
    stop("`correct_min` must be either TRUE or FALSE.", call. = FALSE)
  }

  tau <- period
  peak_time <- -phi * tau / (2 * pi)

  if (correct_period) {
    peak_time <- peak_time %% tau

    if (correct_min) {
      below_min <- !is.na(peak_time) & peak_time < min_time
      peak_time[below_min] <- peak_time[below_min] +
        ceiling((min_time - peak_time[below_min]) / tau) * tau
    }
  }

  peak_time
}
