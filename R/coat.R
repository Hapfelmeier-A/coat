#' Conditional Method Agreement Trees (COAT)
#'
#' Tree models capturing the dependence of method agreement on covariates.
#' The classic Bland-Altman analysis is used for modeling method agreement
#' while the covariate dependency can be learned either nonparametrically
#' via conditional inference trees (CTree) or using model-based recursive
#' partitioning (MOB).
#'
#' @param formula symbolic description of the model of type \code{y1 + y2 ~ x1 + ... + xk}.
#' The left-hand side should specify a pair of measurements (\code{y1} and \code{y2}) for the Bland-Altman analysis.
#' The right-hand side can specify any number of potential split variables for the tree.
#' @param data,subset,na.action arguments controlling the formula processing
#' via \code{\link[stats]{model.frame}}.
#' @param weights optional numeric vector of weights (case/frequency weights, by default).
#' @param means logical. Should the intra-individual mean values of measurements
#' be included as potential split variable?
#' @param type character string specifying the type of tree to be fit. Either \code{"ctree"} (default) or \code{"mob"}.
#' @param minsize,minbucket integer. The minimum number of observations in a subgroup.
#' Only one of the two arguments should be used (see also below).
#' @param minsplit integer. The minimum number of observations to consider splitting.
#' Must be at least twice the minimal subgroup size (\code{minsplit} or \code{minbucket}).
#' If set to \code{NULL} (the default) it is set to be at least 2.5 times the minimal
#' subgroup size.
#' @param ... further control arguments, either passed to \code{\link[partykit]{ctree_control}}
#' or \code{\link[partykit]{mob_control}}, respectively.
#'
#' @details Conditional method agreement trees (COAT) employ unbiased
#' recursive partitioning in order to detect and model dependency on covariates
#' in the classic Bland-Altman analysis. One of two recursive partitioning techniques
#' can be used to find subgroups defined by splits in covariates to a pair
#' of measurements, either nonparametric conditional inference trees (CTree)
#' or parametric model-based trees (MOB). In both cases, each subgroup is associated
#' with two parameter estimates: the mean of the measurement difference (\dQuote{Bias})
#' and the corresponding sample standard deviation (\dQuote{SD}) which can be
#' used to construct the limits of agreement (i.e., the corresponding confidence intervals).
#'
#' The minimum number of observations in a subgroup defaults to 10,
#' so that the mean and variance of the measurement differences can be estimated
#' reasonably for the Bland-Altman analysis. The default can be changed with
#' with the argument \code{minsize} or, equivalently, \code{minbucket}.
#' (The different names stem from slightly different conventions in the underlying
#' tree functions.) Consequently, the minimum number of observations to consider
#' splitting (\code{minsplit}) must be, at the very least, twice the minimum number
#' of observations per subgroup (which would allow only one possible split, though).
#' By default, \code{minsplit} is 2.5 times \code{minsize}.
#' Users are encouraged to consider whether for their application it is sensible
#' to increase or decrease these defaults. Finally, further control parameters
#' can be specified through the \code{...} argument, see
#' \code{\link[partykit]{ctree_control}} and \code{\link[partykit]{mob_control}},
#' respectively, for details.
#'
#' In addition to the standard specification of the two response measurements in the
#' formula via \code{y1 + y2 ~ ...}, it is also possible to use \code{y1 - y2 ~ ...}.
#' The latter may be more intuitive for users that think of it as a model for the
#' difference of two measurements. Finally \code{cbind(y1, y2) ~ ...} also works.
#' Internally, all of these are processed in the same way, namely as a bivariate
#' dependent variable that can then be modeled and plotted appropriately.
#'
#' To add the means of the measurement pair as a potential splitting variable,
#' there are also different equivalent strategies. The standard specification would
#' be via the \code{means} argument: \code{y1 + y2 ~ x1 + ..., means = TRUE}.
#' Alternatively, the user can also extend the formula argument via
#' \code{y1 + y2 ~ x1 + ... + means(y1, y2)}.
#'
#' The SD is estimated by the usual sample standard deviation in each subgroup,
#' i.e., divided by the sample size \eqn{n - 1}. Note that the inference in the
#' MOB algorithm internally uses the maximum likelihood estimate (divided by \eqn{n})
#' instead so the the fluctuation tests for parameter instability can be applied.
#'
#' @references Karapetyan S, Zeileis A, Henriksen A, Hapfelmeier A (2025).
#' \dQuote{Tree models for assessing covariate-dependent method agreement with an application to physical activity measurements.}
#' Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 74, Issue 3, June 2025, Pages 775–799.
#' \doi{10.1093/jrsssc/qlae077}
#'
#' @examples
#' \dontshow{ if(!requireNamespace("MethComp")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("the MethComp package is required for this example but is not installed")
#'   } else q() }
#' }
#' ## package and data (reshaped to wide format)
#' library("coat")
#' data("scint", package = "MethComp")
#' scint_wide <- reshape(scint, v.names = "y", timevar = "meth", idvar = "item", direction = "wide")
#'
#' ## coat based on ctree() without and with mean values of paired measurements as predictor
#' tr1 <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide)
#' tr2 <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide, means = TRUE)
#'
#' ## display
#' print(tr1)
#' plot(tr1)
#'
#' print(tr2)
#' plot(tr2)
#'
#' ## tweak various graphical arguments of the panel function (just for illustration):
#' ## different colors, nonparametric bootstrap percentile confidence intervals, ...
#' plot(tr1, tp_args = list(
#'   xscale = c(0, 150), linecol = "deeppink",
#'   confint = TRUE, B = 250, cilevel = 0.5, cicol = "gold"
#' ))
#' @return Object of class \code{coat}, inheriting either from \code{constparty} (if \code{\link[partykit]{ctree}}
#' is used) or \code{modelparty} (if \code{\link[partykit]{mob}} is used).
#'
#' @importFrom stats model.weights na.omit update weighted.mean
#' @importFrom partykit ctree_control mob_control
#'
#' @export
coat <- function(formula, data, subset, na.action, weights, means = FALSE, type = c("ctree", "mob"),
  minsize = 10L, minbucket = minsize, minsplit = NULL, ...)
{
  ## type of tree
  type <- match.arg(tolower(type), c("ctree", "mob"))

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## preparation of ctree/mob call
  m <- match.call(expand.dots = FALSE)
  if(missing(na.action)) m$na.action <- na.omit

  ## if desired "means(y1, y2)" is added as split variable
  if(means) {
    formula <- update(formula, . ~ . + `_means_`)
    formula[[3L]][[3L]] <- formula[[2]]
    formula[[3L]][[3L]][[1L]] <- as.name("means")
    m$formula <- formula
  }

  ## if measurements are specified as `y1 - y2`, switch to `y1 + y2` internally
  if(formula[[2L]][[1L]] == as.name("-")) {
    formula[[2L]][[1L]] <- as.name("+")
    m$formula <- formula
  }

  ## update/remove processed arguments
  m$means <- NULL
  m$type <- NULL

  ## process hyperparameters
  if(!missing(minsize) && !missing(minbucket)) {
    warning("the minimal subgroup size should either be specified by 'minsize' or 'minbucket' but not both, using 'minsize'")
    minbucket <- minsize
  }
  minsize <- minbucket
  if(is.null(minsplit)) minsplit <- ceiling(2.5 * minsize)
  if(minsplit < 2L * minsize) {
    warning("the minimal sample size to consider splitting ('minsplit') must be at least twice the minimal subgroup size ('minsize'), increased accordingly")
    minsplit <- 2L * minsize
  }

  ## add fit/trafo function
  if(type == "mob") {
    m[[1L]] <- as.call(quote(partykit::mob))
    m$fit <- bafit
    m$control <- partykit::mob_control(minsize = minsize, minsplit = minsplit, ...)
    m$control$ytype <- "matrix"
  } else {
    m[[1L]] <- as.call(quote(partykit::ctree))
    m$ytrafo <- batrafo
    m$control <- partykit::ctree_control(minbucket = minsize, minsplit = minsplit, ...)
  }

  ## fit tree
  rval <- eval(m, parent.frame())

  ## informative warning if tree considered splitting at all
  if(is.null(rval$node$split) && (is.null(rval$node$info) || is.null(rval$node$info$test))) {
    message("Info: The tree has no splits due to the hyperparameters ('minsize', 'minsplit', ...), no test were carried out, possibly consider adjusting the hyperparameters.")
  }

  ## unify output
  rval$info$call <- cl
  class(rval) <- c("coat", class(rval))
  if(type == "mob") {
    rval$fitted[["(weights)"]] <- model.weights(rval$data)
    if(is.null(rval$fitted[["(weights)"]])) rval$fitted[["(weights)"]] <- 1
    rval$fitted[["(response)"]] <- rval$data[, attr(rval$info$terms$response, "term.labels"), drop = FALSE]
  }
  return(rval)
}
