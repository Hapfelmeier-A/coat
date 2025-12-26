#' Bland-Altman Test of Method Agreement
#'
#' Function to perform a Bland-Altman test of differences in method agreement
#' with single measurements per subject (or item). Additional functions are
#' given for printing and plotting. The function \code{\link[coat]{coat}} must
#' be used in case of replicate measurements, setting the arguments
#' \code{alpha = 1} and \code{maxdepth = 1L}.
#'
#' @param formula symbolic description of the model used to perform the Bland-Altman test of type \code{y ~ x}.
#' The left-hand side should specify the measurements (\code{y}) to assess the agreement.
#' The right-hand side should specify a factor with two levels indicating two independent groups or samples to be compared. Alternatively, multilevel factors or continuously scaled variables can be specified to perform a Bland-Altman test of association, followed by binary splitting into two subgroups.
#' @param data,subset,na.action arguments controlling the formula processing
#' via \code{\link[stats]{model.frame}}. \code{data} must be provided in long format,
#' i.e. with two rows per subject (or item), one for each of the measurements made with the compared methods.
#' @param weights optional numeric vector of weights (case/frequency weights, by default).
#' @param x an object as returned by \code{\link[coat]{batest}}.
#' @param digits a numeric specifying the number of digits to display.
#' @param type character string specifying whether \code{"test"} statistics (default), the \code{"model"} or \code{"both"} should be printed.
#' @param ... further control arguments, passed to \code{\link[partykit]{ctree_control}}
#'
#' @details The minimum number of subjects (or items) in each of the two groups
#' or samples that are passed to \code{batest} or defined by it must be three.
#'
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
#' ## package and data
#' library("coat")
#' data("scint", package = "MethComp")
#' scint$DMSA <- factor(scint$meth == "DMSA")
#'
#' ## two-sample BA-test
#' testresult <- batest(y ~ sex, data = scint, id = "item", meth = "DMSA")
#'
#' ## display
#' testresult
#' print(testresult, digits = 1, type = "both")
#' plot(testresult)
#'
#' @return Object of class \code{batest} with elements
#' \item{\code{test}}{result of the Bland-Altman test.}
#' \item{\code{model}}{tree model used to perform the Bland-Altman test.}
#'
#' @importFrom partykit character_split sctest.constparty
#'
#' @export
batest <- function(formula, data, id = NULL, meth = NULL, subset, na.action, weights, ...)
{
  ## check whether a single covariate is used.
  if (length(formula[[3]]) > 1) {
    stop("Please provide a single variable on the right-hand side of the formula. Otherwise, consider using the coat() function.")
  }

  ## preparation of ctree call
  base <- match.call(expand.dots = TRUE)
  base[[1L]] <- quote(partykit::ctree)
  base$data <- coat.reshape(formula = formula, data = data, id = id, meth = meth, replicates = FALSE, paired = FALSE)
  base$formula <- update.formula(formula, as.formula("`_diff_` ~ ."))
  if (missing(na.action)) base$na.action <- na.omit

  ## add fit/trafo function
  m <- m.bias <- m.var <- base
  m$ytrafo <- batrafo
  m.bias$ytrafo <- batrafo.mean
  m.var$ytrafo <- batrafo.var

  ctrl <- ctree_control(alpha = 1, minsplit = 6L, minbucket = 3L, maxdepth = 1L, ...)
  m$control <- m.bias$control <- m.var$control <- ctrl

  ## fit tree
  rval <- eval(m, parent.frame())
  rval.bias <- eval(m.bias, parent.frame())
  rval.var <- eval(m.var, parent.frame())

  ## extract test statistics
  test <- matrix(NA, nrow = 3, ncol = 5)
  colnames(test) <- c(character_split(rval$node$split, data = rval$data)$levels, "Chisq", "df", "p-value")
  rownames(test) <- c("Bias", "SD", "Total")

  test[1, c(3, 5)] <- sctest.constparty(rval.bias, node = 1L)
  test[2, c(3, 5)] <- sctest.constparty(rval.var, node = 1L)
  test[3, c(3, 5)] <- sctest.constparty(rval, node = 1L)
  test[1:2, 1:2] <- t(coef.coat(rval))
  test[, 4] <- c(1, 1, 2)

  ## unify output
  rval$info$call <- match.call(expand.dots = FALSE)
  class(rval) <- c("coat", class(rval))
  trlist <- list("test" = test, "model" = rval)
  class(trlist) <- "batest"
  return(trlist)
}

#' @describeIn batest function to print the result of the Bland-Altman test.
#' @export
print.batest <- function(x, digits = 2, type = c("test", "model", "both"), ...) {

  ## type of output
  type <- match.arg(tolower(type), c("test", "model", "both"))

  if(type != "model") {
    out <- x$test
    out[, 1] <- round(out[, 1], digits)
    out[, 2] <- round(out[, 2], digits)
    out[, 3] <- round(out[, 3], 3)
    out[, 5] <- round(out[, 5], 3)

    if (type == "test") {
      print(out, na.print="", ...)
    } else {
      print(x$model)
      cat("\n")
      print(out, na.print="", ...)
    }
  } else print(x$model)
}

#' @describeIn batest function to plot the result of the Bland-Altman test.
#' @export
plot.batest <- function(x, ...) {
  plot(x$model, ...)
}

