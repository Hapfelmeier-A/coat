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
  if(is.null(x)) x <- cbind("(Intercept)" = rep.int(1, length(y)))
  m <- lm.fit(x, y)
  b <- m$coefficients
  s2 <- mean(m$residuals^2)
  list(
    coefficients = c(b, "(Variance)" = s2),
    objfun = -sum(dnorm(y, mean = m$fitted.values, sd = sqrt(s2), log = TRUE)),
    estfun = if(estfun) cbind(m$residuals/s2 * x, "(Variance)" = (m$residuals^2 - s2)/(2 * s2^2)) else NULL,
    object = if(object) lm(y ~ 0 + x) else NULL
  )
}


#' Transformation function used in \code{\link[partykit]{ctree}} to transform a bivariate measurement pair to mean and variance of the differences. See \code{ytrafo} in \code{\link[partykit]{ctree}} for details.
#' @noRd
batrafo <- function(data, weights, control, ...) {
  y <- data$data[, data$variables$y, drop = TRUE]
  if(NCOL(y) != 2L) stop("'y' data must provide exactly 2 columns")
  y <- y[, 1L] - y[, 2L]
  function(subset, weights, info, estfun, object, ...) {
    list(estfun = cbind(mean = y, var = (y - mean(y))^2), unweighted = TRUE)
  }
}

#' \code{fit} function used in \code{\link[partykit]{mob}} for modeling the differences of a bivariate measurement pair by an intercept-only linear regression model with residual variance as additional parameter. See \code{fit} in \code{\link[partykit]{mob}} for details.
#' @noRd
bafit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE)
{
  if(NCOL(y) != 2L) stop("'y' must provide exactly 2 columns")
  y <- y[, 1L] - y[, 2L]
  if(is.null(x)) x <- cbind("(Intercept)" = rep.int(1, length(y)))
  m <- lm.fit(x, y)
  b <- m$coefficients
  s2 <- mean(m$residuals^2)
  list(
    coefficients = c(b, "(Variance)" = s2),
    objfun = -sum(dnorm(y, mean = m$fitted.values, sd = sqrt(s2), log = TRUE)),
    estfun = if(estfun) cbind(m$residuals/s2 * x, "(Variance)" = (m$residuals^2 - s2)/(2 * s2^2)) else NULL,
    object = if(object) lm(y ~ 0 + x) else NULL
  )
}

diffs <- function(y1, y2) {
  if(NCOL(y1) != 1L || NCOL(y2) != 1L || NROW(y1) != NROW(y2)) {
    stop("'y1' and 'y2' must be a pair of univariate measurements of the same length")
  }
  if(!is.numeric(y1)) y1 <- as.numeric(y1)
  if(!is.numeric(y2)) y2 <- as.numeric(y2)
  return(y1 - y2)
}

means <- function(y1, y2) {
  if(NCOL(y1) != 1L || NCOL(y2) != 1L || NROW(y1) != NROW(y2)) {
    stop("'y1' and 'y2' must be a pair of univariate measurements of the same length")
  }
  if(!is.numeric(y1)) y1 <- as.numeric(y1)
  if(!is.numeric(y2)) y2 <- as.numeric(y2)
  return((y1 + y2)/2)
}
