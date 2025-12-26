#' Convenience Function for Bland-Altman Analysis
#'
#' Auxiliary function for preparing data in the long format for Bland-Altman analysis.
#'
#' @param formula symbolic description of the model of type \code{y ~ x1 + ... + xk}.
#' The left-hand side should specify the measurements (\code{y}) for the Bland-Altman analysis.
#' The right-hand side can specify any number of potential split variables for the tree.
#' @param data in long format, i.e. with two or more rows per subject (or item)
#' if one or more measurements are available per method.
#' @param id character referring to a column in \code{data} which indicates the data rows belonging to the same subject (or item).
#' @param meth character referring to a column in \code{data} which indicates the data rows belonging to the same method.
#' @param replicates Does \code{data} contain replicate measurements?
#' @param paired Are replicate measurements paired (TRUE) or unpaired (FALSE)?
#'
#' @return Data frame including the subject-specific components of a
#' Bland-Altman analysis and the covariates required for modelling in \code{\link[coat]{coat.fit}}.
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
#' d1 <- coat.reshape(y ~ age + sex, data = scint, id = "item", meth = "DMSA", replicates = FALSE, paired = FALSE)
#' head(d1)
#'
#'
#' ### replicate measurements per subject (or item)
#'
#' ## data
#' data("ox", package = "MethComp")
#'
#' ## coat with paired measurements and mean values as only predictor
#' d2 <- coat.reshape(y ~ 1, data = ox, id = "item", meth = "meth", replicates = TRUE, paired = TRUE)
#' head(d2)
#'

#' @importFrom stats lm.fit dnorm lm var sd

#' @rdname coat.reshape
#' @export
coat.reshape <- function(formula, data, id = NULL, meth = NULL, replicates = FALSE, paired = FALSE){
  y <- all.vars(formula)[1]
  x <- all.vars(formula)[-1]
  data[, y] <- as.numeric(data[, y])

  if (replicates) {
    if (!paired) {
      data <- data.frame(cbind(unique(data[, id]),
                               c(by(data[, c(y, meth)], data[, id], function(z) diff(by(z[, y], z[, meth], mean, na.rm = TRUE)))),
                               c(by(data[, y], data[, id], mean, na.rm = TRUE)),
                               data[!duplicated(data[, id]), x, drop = FALSE],
                               do.call(cbind, by(data[, c(y, id)], data[, meth], function(z1) c(by(z1[, y], z1[, id], function(z2) length(na.omit(z2)))))),
                               do.call(cbind, by(data[, c(y, id)], data[, meth], function(z1) c(by(z1[, y], z1[, id], function(z2) sum(aov(z2 ~ 1)$residuals^2)))))
      ))
      names(data) <- c("id", "_mean_diff_", "_means_", x, "nr1", "nr2", "rss1", "rss2")
    } else {
      data <- data.frame(cbind(unique(data[, id]),
                               c(by(data[, c(y, meth)], data[, id], function(z) mean(apply(do.call(cbind, split(z[, y], z[, meth])), 1, diff)))),
                               c(by(data[, y], data[, id], mean)),
                               data[!duplicated(data[, id]), x, drop = FALSE],
                               as.numeric(table(data[, id]) / 2),
                               c(by(data[, c(y, meth)], data[, id], function(z) sum(aov(apply(do.call(cbind, split(z[, y], z[, meth])), 1, diff) ~ 1)$residuals^2)))
      ))
      names(data) <- c("id", "_mean_diff_", "_means_", x, "nr", "rss")
    }
  }
  else {
    data <- data.frame(cbind(unique(data[, id]),
                             c(by(data[, y], data[, id], diff, na.rm = TRUE)),
                             c(by(data[, y], data[, id], mean, na.rm = TRUE)),
                             data[!duplicated(data[, id]), x, drop = FALSE]
    ))
    names(data) <- c("id", "_diff_", "_means_", x)
  }
  return(data)
}

#' Transformation function used in \code{\link[partykit]{ctree}} to model the mean and variance as bivariate outcome. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
meanvar <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = cbind(mean = y, var = (y - mean(y))^2), unweighted = TRUE)
  }
}

#' Transformation function used in \code{\link[partykit]{ctree}} to model the variance as outcome. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
.var <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = (y - mean(y))^2, unweighted = TRUE)
  }
}

#' \code{fit} function used in \code{\link[partykit]{mob}} for an intercept-only linear regression model with residual variance as additional parameter. See \code{fit} in \code{\link[partykit]{mob}} for details.
#' @noRd
gaussfit <- function(y, x = NULL, start = NULL, weights = NULL,
                     offset = NULL, ..., estfun = FALSE, object = FALSE)
{
  x <- cbind("(Intercept)" = rep.int(1, length(y)))
  m <- lm.fit(x, y)
  b <- m$coefficients
  s2 <- mean(m$residuals^2)
  s2ols <- var(m$residuals) ## n * s2 = (n-1) * s2ols
  list(
    coefficients = c("Bias" = as.numeric(b), "SD" = sqrt(s2ols)), ## parameter estimates employed for Bland-Altman analysis
    mobcoefficients = c(b, "(Variance)" = s2), ## parameter estimates underlying the mob() inference
    objfun = -sum(dnorm(y, mean = m$fitted.values, sd = sqrt(s2), log = TRUE)),
    estfun = if(estfun) cbind(m$residuals/s2 * x, "(Variance)" = (m$residuals^2 - s2)/(2 * s2^2)) else NULL,
    object = if(object) lm(y ~ 0 + x) else NULL
  )
}

#' Transformation function used in \code{\link[partykit]{ctree}} to transform a bivariate measurement pair to mean and variance of the differences. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
batrafo <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  n <- nrow(data$data)
  variance <- (y - mean(y))^2 * (n/(n - 1))
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = cbind(mean = y, var = variance), unweighted = TRUE)
  }
}

#' Transformation function used in \code{\link[partykit]{ctree}} to transform a bivariate measurement pair to mean of the differences. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
batrafo.mean <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = y, unweighted = TRUE)
  }
}

#' Transformation function used in \code{\link[partykit]{ctree}} to transform a bivariate measurement pair to variance of the differences. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
batrafo.var <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  n <- nrow(data$data)
  variance <- (y - mean(y))^2 * (n/(n - 1))
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = variance, unweighted = TRUE)
  }
}

#'Transformation function for unpaired measurements
batrafo.repl.unpair <- function(data, weights, control, ...){
  y <- data$data[, data$variables$y, drop = TRUE]
  n <- nrow(data$data)
  nr1 <- data$data$nr1
  nr2 <- data$data$nr2
  tau <- (y - mean(y))^2 * (n/(n - 1)) # between subject variance
  s1 <- data$data$rss1 / (mean(nr1) - 1) # within subject variance
  s2 <- data$data$rss2 / (mean(nr2) - 1)
  overallvar <- tau + (1 - (1/n) * sum(1/nr1)) * s1 + (1 - (1/n) * sum(1/nr2)) * s2
  function(subset, weights, info, estfun, object, ...){
    list(estfun = cbind(Bias = y, Variance = overallvar), unweighted = TRUE)
  }
}

#' Transformation function for paired measurements
batrafo.repl.pair <- function(data, weights, control, ...){
  y <- data$data[, data$variables$y, drop = TRUE]
  n <- nrow(data$data)
  nr <- data$data$nr
  bias <- y * nr * (n / sum(nr))
  tau <- (y - weighted.mean(y, nr))^2 * nr * (n/(n - 1)) # between subject variance
  sigma <- data$data$rss / (mean(nr) - 1) # within subject variance
  overallvar <- (tau - sigma) / ((sum(nr)^2 - sum(nr^2)) / ((n - 1) * sum(nr))) + sigma
  function(subset, weights, info, estfun, object, ...){
    list(estfun = cbind(Bias = bias, Variance = overallvar), unweighted = TRUE)
  }
}


