#' Bland-Altman Test of Method Agreement
#'
#' Function to perform a Bland-Altman test for hypothesis testing of differences in method agreement. Additional functions are given for printing and plotting.
#'
#' @param formula symbolic description of the model used to perform the Bland-Altman test of type \code{y1 + y2 ~ x}.
#' The left-hand side should specify a pair of measurements (\code{y1} and \code{y2}) to assess the agreement.
#' The right-hand side should specify a factor with two levels indicating two independent groups or samples to be compared. Alternatively, multilevel factors or continuously scaled variables can be specified to perform a Bland-Altman test of association, followed by binary splitting into two subgroups.
#' @param data,subset,na.action arguments controlling the formula processing
#' via \code{\link[stats]{model.frame}}.
#' @param weights optional numeric vector of weights (case/frequency weights, by default).
#' @param x an object as returned by \code{\link[coat]{BAtest}}.
#' @param digits a numeric specifying the number of digits to display.
#' @param type character string specifying whether \code{"test"} statistics (default), the \code{"model"} or \code{"both"} should be printed.
#' @param ... further control arguments, passed to \code{\link[partykit]{ctree_control}}
#'
#' @examples
#' \dontshow{ if(!requireNamespace("MethComp")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("the MethComp package is required for this example but is not installed")
#'   } else q() }
#' }
#' ## package and data
#' library("coat")
#' data("VitCap", package = "MethComp")
#' ## transform data to required 'wide' format
#' VitCap_wide <- reshape(VitCap, v.names = "y", timevar = "instrument",
#'                        idvar = c("item", "user"), drop = "meth", direction = "wide")
#'
#' ## two-sample BA-test
#' testresult <- BAtest(y.St + y.Exp ~ user, data = VitCap_wide)
#'
#' ## display
#' testresult
#' print(testresult, digits = 1, type = "both")
#' plot(testresult)
#'
#' @return Object of class \code{BAtest} with elements
#' \item{\code{test}}{result of the Bland-Altman test.}
#' \item{\code{model}}{tree model used to perform the Bland-Altman test.}
#'
#' @export
BAtest <- function(formula, data, subset, na.action, weights, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)
  
  # check whether a single covariate is used.
  if (length(formula[[3]]) > 1) {
    stop("Please provide a single variable on the right-hand side of the formula. Otherwise, consider using the coat() function.")
  }
  
  ## preparation of ctree call
  m <- m.bias <- m.var <- match.call(expand.dots = FALSE)
  if(missing(na.action)) m$na.action <- m.bias$na.action <- m.var$na.action <- na.omit
  
  ## add fit/trafo function
  m[[1L]] <- m.bias[[1L]] <- m.var[[1L]] <- as.call(quote(partykit::ctree))
  m$ytrafo <- batrafo
  m.bias$ytrafo <- batrafo.diff
  m.var$ytrafo <- batrafo.var
  m$control <- m.bias$control <- m.var$control <- partykit::ctree_control(alpha = 1, minsplit = 6L, minbucket = 3L, maxdepth = 1L, ...)
  
  ## fit tree
  rval <- eval(m, parent.frame())
  rval.bias <- eval(m.bias, parent.frame())
  rval.var <- eval(m.var, parent.frame())
  
  # extract test statistics
  test <- matrix(NA, nrow = 3, ncol = 5)
  colnames(test) <- c(character_split(rval$node$split, data = rval$data)$levels, "Chisq", "df", "p-value")
  rownames(test) <- c("Bias", "SD", "Total")
  
  test[1, c(3, 5)] <- sctest.constparty(rval.bias, node = 1L)
  test[2, c(3, 5)] <- sctest.constparty(rval.var, node = 1L)
  test[3, c(3, 5)] <- sctest.constparty(rval, node = 1L)
  test[1:2, 1:2] <- t(coat:::coef.coat(rval))
  test[2, 1:2] <- sqrt(test[2, 1:2])
  test[, 4] <- c(1, 1, 2)
  
  ## unify output
  rval$info$call <- cl
  class(rval) <- c("coat", class(rval))
  trlist <- list("test" = test, "model" = rval)
  class(trlist) <- "BAtest"
  return(trlist)
}

#' @describeIn BAtest function to print the result of the Bland-Altman test.
#' @export
print.BAtest <- function(x, digits = 2, type = c("test", "both", "model"), ...) {
  
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

#' @describeIn BAtest function to plot the result of the Bland-Altman test.
#' @export
plot.BAtest <- function(x, ...) {
  plot(x$model, ...)
}
