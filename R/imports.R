# Selective namespace imports used throughout mixedcirc.
#
# Keeping these imports explicit prevents dependencies from being attached to
# the user's search path and avoids ambiguous symbols from broad imports.
#'
#' @importFrom emmeans emmeans emtrends test
#' @importFrom doFuture %dofuture%
#' @importFrom foreach %do% foreach
#' @importFrom graphics abline axis hist image layout lines mtext pairs par plot plot.new rect text title
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom lme4 getME isGLMM isLMM refit
#' @importFrom methods is new show
#' @importClassesFrom Biobase ExpressionSet
#' @importClassesFrom limma EList
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom rlang .data
#' @importFrom stats anova approxfun as.dendrogram as.dist as.formula ave cancor coef cov2cor density dist fitted fitted.values hclust lm median model.frame model.matrix na.exclude na.omit order.dendrogram p.adjust plogis predict pt qt quantile rbinom reorder residuals rnorm runif sd setNames sigma simulate terms update var vcov
#' @importFrom utils combn object.size str
NULL

# Non-standard evaluation symbols used by data.table and foreach.
utils::globalVariables(c(
  ".", ".N", ".SD", ".mixedcirc_weight", "chIndx", "end", "i",
  "start", "strand"
))
