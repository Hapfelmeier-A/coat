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
