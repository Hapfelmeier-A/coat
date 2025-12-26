#' Conditional Method Agreement Trees (COAT)
#'
#' Tree models capturing the dependence of method agreement on covariates.
#' The Bland-Altman analysis for single or replicate measurements per subject (or item)
#' is used to model method agreement, while the covariate dependency can
#' be learned either nonparametrically via conditional inference trees (CTree)
#' or using model-based recursive partitioning (MOB).
#'
#' @param formula symbolic description of the model of type \code{y ~ x1 + ... + xk}.
#' The left-hand side should specify the measurements (\code{y}) for the Bland-Altman analysis.
#' The right-hand side can specify any number of potential split variables for the tree.
#' @param data,subset,na.action arguments controlling the formula processing
#' via \code{\link[stats]{model.frame}}. \code{data} must be provided in long format,
#' i.e. with two or more rows per subject (or item) if one or more measurements
#' are available per method.
#' @param id character referring to a column in \code{data} which indicates the data rows belonging to the same subject (or item).
#' @param meth character referring to a column in \code{data} which indicates the data rows belonging to the same method.
#' @param weights optional numeric vector of weights (case/frequency weights, by default).
#' @param means logical. Should the intra-individual mean values of measurements
#' be included as potential split variable?
#' @param replicates Does \code{data} contain replicate measurements?
#' @param paired Are replicate measurements paired (TRUE) or unpaired (FALSE)?
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
#' The SD is estimated by the usual sample standard deviation in each subgroup,
#' i.e., divided by the sample size \eqn{n - 1}. Note that the inference in the
#' MOB algorithm internally uses the maximum likelihood estimate (divided by \eqn{n})
#' instead so the the fluctuation tests for parameter instability can be applied.
#'
#' An extension of the Bland-Altman analysis addresses the frequent
#' cases of paired and unpaired replicate measurements per subject or item.
#' Both cases are analysed in COAT using CTree with appropriate data aggregation
#' to model Bias and SD at the subject/item level, with the latter being
#' decomposed into within-subject variance and between-subject variance.
#' Currently, only covariates that are constant between replicate measurements
#' can be used, and the first entry in \code{data} is used for model fitting.
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
#' Mean values of measurements can be added as a potential splitting variable
#' via the \code{means} argument, using the standard specification
#' \code{y ~ x1 + ..., means = TRUE}. It may not be appropriate to calculate
#' mean values across paired replicate measurements, as in such cases it is
#' often assumed that the underlying true values vary.
#'
#' @references Karapetyan S, Zeileis A, Henriksen A, Hapfelmeier A (2025).
#' \dQuote{Tree models for assessing covariate-dependent method agreement with an application to physical activity measurements.}
#' Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 74, Issue 3, June 2025, Pages 775–799.
#' \doi{10.1093/jrsssc/qlae077}
#' @references Karapetyan S, Zeileis A, Flick M, Saugel B, Hapfelmeier A (2026).
#' \dQuote{Tree models for covariate-dependent method agreement with repeated measurements.}
#' TBD, Volume XXX, Issue XXX, XXX 2026, Pages XXX–XXX.
#' \doi{XXX}
#'
#'
#' @examples
#' \dontshow{ if(!requireNamespace("MethComp")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("the MethComp package is required for this example but is not installed")
#'   } else q() }
#' }
#' ### single measurements per subject (or item)
#'
#' ## package and data
#' library("coat")
#' data("scint", package = "MethComp")
#' scint$DMSA <- factor(scint$meth == "DMSA")
#'
#' ## coat based on ctree()
#' tr1 <- coat(y ~ age + sex, data = scint, id = "item", meth = "DMSA")
#'
#' ## coat based on mob() including mean values of paired measurements as predictor
#' tr2 <- coat(y ~ age + sex, data = scint, id = "item", meth = "DMSA", means = TRUE, type = "mob")
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
#'
#'
#' ### replicate measurements per subject (or item)
#'
#' ## data
#' data("ox", package = "MethComp")
#'
#' ## coat with paired measurements and mean values as only predictor
#' tr3 <- coat(y ~ 1, data = ox, id = "item", meth = "meth", replicates = TRUE, paired = TRUE, means = TRUE)
#'
#' ## same coat, but treating measurements as unpaired (exchangeable)
#' tr4 <- coat(y ~ 1, data = ox, id = "item", meth = "meth", replicates = TRUE, paired = FALSE, means = TRUE)
#'
#' ## display
#' print(tr3)
#' plot(tr3)
#'
#' print(tr4)
#' plot(tr4)
#'
#' @return Object of class \code{coat}, inheriting either from \code{constparty} (if \code{\link[partykit]{ctree}}
#' is used) or \code{modelparty} (if \code{\link[partykit]{mob}} is used).
#'
#' @importFrom stats model.weights na.omit update weighted.mean
#' @importFrom partykit ctree_control mob_control
#'
#' @export
coat <- function(formula, data, id = NULL, meth = NULL, subset, na.action, weights, means = FALSE,
                 replicates = FALSE, paired = FALSE, type = c("ctree", "mob"), minsize = 10L,
                 alpha = 0.05, minbucket = minsize, minsplit = NULL, ...){
  cl <- match.call()

  if(is.null(id) | is.null(meth)) stop("arguments 'id' and 'meth' must be specified")

  if(paired & means) message("Info: Use of mean values across paired replicate measurements may not be meaningful")

  # if(!missing(subset)) data <- data[eval(substitute(subset), envir = data), , drop = FALSE]

  data <- coat.reshape(formula, data, id = id, meth = meth, replicates = replicates, paired = paired)
  fit <- coat.fit(formula, data, replicates = replicates, paired = paired, type = type, means = means, na.action = na.action,
                  minsize = minsize, minsplit = minsplit, alpha = alpha, ...)

  fit$info$call <- cl
  return(fit)
}


#' @describeIn coat fitter function for tree models.
#' @export
coat.fit <- function(formula, data, subset, na.action, weights, means = FALSE, type = c("ctree", "mob"),
                     replicates = FALSE, paired = FALSE, minsize = 10L, alpha = 0.05, minbucket = minsize, minsplit = NULL, ...){

  ## set replicates with paired measurements
  if(paired) replicates <- TRUE

  ## type of tree
  type <- match.arg(tolower(type), c("ctree", "mob"))
  if (replicates & type != "ctree") {
    type <- "ctree"
    warning("'type' was set to 'ctree' with replicate measurements")
  }

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## preparation of ctree/mob call
  m <- match.call(expand.dots = FALSE)
  if(missing(na.action)) m$na.action <- na.omit

  ## update formula
  if (replicates) {
    if (!paired) {
      formula <- as.formula(paste("`_mean_diff_` ~", paste(c("id", "nr1", "nr2", "rss1", "rss2", "`_means_`"), collapse = "+"), "|", formula[3], ifelse(means, "+ `_means_`", "")))
    } else {
      formula <- as.formula(paste("`_mean_diff_` ~", paste(c("id", "nr", "rss", "`_means_`"), collapse = "+"), "|", formula[3], ifelse(means, "+ `_means_`", "")))
    }
  } else {
    formula <- as.formula(paste("`_diff_` ~ id + `_means_`|", formula[3], ifelse(means, "+ `_means_`", "")))
  }
  m$formula <- formula

  ## update/remove processed arguments
  m$means <- NULL
  m$type <- NULL

  ## process hyperparameters
  if(!missing(minsize) && !missing(minbucket)) {
    warning("the minimum subgroup size should be specified using either 'minsize' or 'minbucket'. Try using only 'minsize'.")
    minbucket <- minsize
  }
  minsize <- minbucket
  if(is.null(minsplit)) minsplit <- ceiling(2.5 * minsize)
  if(minsplit < 2L * minsize) {
    warning("The minimum sample size for considering further splitting ('minsplit') must be at least twice the minimum subgroup size ('minsize') and has been increased accordingly.")
    minsplit <- 2L * minsize
  }

  ## add trafo function for ctree
  m[[1L]] <- as.call(quote(partykit::ctree))
  if(replicates) {
    if(paired) {m$ytrafo <- batrafo.repl.pair
    } else {m$ytrafo <- batrafo.repl.unpair}
  } else {m$ytrafo <- batrafo}
  m$control <- partykit::ctree_control(minbucket = minsize, minsplit = minsplit, alpha = alpha, ...)

  ## add fit function for mob
  if(type == "mob") {
    m[[1L]] <- as.call(quote(partykit::mob))
    m$fit <- gaussfit
    m$control <- partykit::mob_control(minsize = minsize, minsplit = minsplit, alpha = alpha, ...)
    m$control$ytype <- "matrix"
  }

  ## fit tree
  rval <- eval(m, parent.frame())

  ## informative warning if tree considered splitting at all
  if(is.null(rval$node$split) && (is.null(rval$node$info) || is.null(rval$node$info$test))) {
    message("Info: No splits were examined and performed due to the selected hyperparameters (‘minsize’, ‘minsplit’, ...). Consider other settings.")
  }

  ## unify output
  rval$info$call <- cl
  class(rval) <- c("coat", class(rval))

  # for the further functions we add an additional class paired/unpaired
  if(replicates) {
    if(paired) {class(rval) <- c(class(rval), "paired")
    } else class(rval) <- c(class(rval), "unpaired")
  }

  return(rval)
}

