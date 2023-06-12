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
#' @param ... further control arguments, either passed to \code{\link[partykit]{ctree_control}}
#' or \code{\link[partykit]{mob_control}}, respectively.
#'
#' @details The minimum number of observations required to model conditional
#' agreement defaults to 20. Users may choose to modify this value as needed. See
#' \code{\link[partykit]{ctree_control}} and \code{\link[partykit]{mob_control}} for details.
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
#' @return Object of class \code{coat}, inheriting either from \code{constparty} (if \code{\link[partykit]{ctree}}
#' is used) or \code{modelparty} (if \code{\link[partykit]{mob}} is used).
#'
#' @importFrom stats model.weights na.omit update weighted.mean
#' @importFrom partykit ctree_control mob_control
#'
#' @export
coat <- function(formula, data, subset, na.action, weights, means = FALSE, type = c("ctree", "mob"), ...)
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

  ## update/remove processed arguments
  m$means <- NULL
  m$type <- NULL

  ## add fit/trafo function
  if(type == "mob") {
    m[[1L]] <- as.call(quote(partykit::mob))
    m$fit <- bafit
    m$control <- partykit::mob_control(...)
    m$control$ytype <- "matrix"
    if(is.null(m$control$minsize) && is.null(m$control$minbucket)) m$control$minsize <- 20L
  } else {
    m[[1L]] <- as.call(quote(partykit::ctree))
    m$ytrafo <- batrafo
    m$control <- partykit::ctree_control(...)
  }

  ## fit tree
  rval <- eval(m, parent.frame())
  
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
