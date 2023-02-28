#' Conditional Method Agreement Forests
#'
#' Function to model conditional method agreement by \code{\link[disttree]{distforest}}. See \code{\link[COAT]{coat}} for a respective tree implementation.
#'
#' @param y1 a character string specifying the variable in data containing the measurements by one method.
#' @param y2 a character string specifying the variable in data containing the measurements by another method.
#' @param covars a character string or vector of a single or multiple covariates.
#' @param data a data frame containing \code{y1}, \code{y2} and \code{covars}.
#' @param means a logical indicating whether intraindividual mean values of measurements shall be included as covariate.
#' @param ... further arguments passed to \code{\link[disttree]{distforest}}.
#'
#' @details See \code{\link[disttree]{distforest}} for implementation details.
#'
#' @examples
#' ### load package ###
#' library("COAT")
#'
#' ### data ###
#' require(MethComp)
#' data(scint)
#' ## transform data to required 'wide' format
#' scint_wide <- reshape(scint, v.names = "y", timevar = "meth", idvar = "item", direction = "wide")
#'
#' ### model fit ###
#' myfor <- coarf("y.DTPA", "y.DMSA", c("age", "sex"), data = scint_wide)
#'
#' ### prediction of fitted bias and standard deviation of differences ###
#' predict(myfor)
#'
#' ### prediction and visualization of age dependent method agreement by sex ###
#' preds_females <- predict(myfor, newdata = data.frame(sex = "F", age = sort(scint_wide$age)))
#' preds_males <- predict(myfor, newdata = data.frame(sex = "M", age = sort(scint_wide$age)))
#' matplot(sort(scint_wide$age), cbind(preds_females$mu,
#'                                     preds_females$mu - qnorm(0.975) * preds_females$sigma,
#'                                     preds_females$mu + qnorm(0.975) * preds_females$sigma),
#'         type = "l", lty = c(1, 2, 2), col = "black", xlab = "age", ylab = "Difference", las = 1)
#' matlines(sort(scint_wide$age), cbind(preds_males$mu,
#'                                      preds_males$mu - qnorm(0.975) * preds_females$sigma,
#'                                      preds_males$mu + qnorm(0.975) * preds_females$sigma),
#'          lty = c(1, 2, 2), col = "red", xlab = "age", ylab = "Difference", las = 1)
#' rug(scint_wide$age[scint_wide$sex == "F"])
#' rug(scint_wide$age[scint_wide$sex == "M"], col = "red")
#' legend("topleft", c("F", "M"), lty = 1, col = c("black", "red"), bty = "n")
#'
#' @return Object of class \code{distforest}.
#'
#' @importFrom stats as.formula dnorm lm lm.fit qnorm sd
#'
#' @export
coarf <- function(y1, y2, covars, data, means = FALSE, ...) {

  # create model data
  moddat <- data.frame(data, "means." = rowMeans(data[, c(y1, y2)]))

  # create outcome
  moddat$diffs. <- moddat[, y1] - moddat[, y2]

  # remove missing values
  if (any(is.na(moddat[, c("diffs.", covars)]))) {
    moddat <- moddat[!apply(moddat[, c("diffs.", covars)], 1, function(x) any(is.na(x))), ]
    warning("Observations with missing values have been removed.")
  }

  # create model formula
  if (means) {
    model.formula <- as.formula(paste("diffs. ~ means. + ", paste(covars, collapse = " + ")))
  } else {
    model.formula <- as.formula(paste("diffs. ~ ", paste(covars, collapse = " + ")))
  }

  # fit the forest model
  model <- disttree::distforest(model.formula, data = moddat, ...)

  # Update data with means
  if (means) model$data$means. <- moddat$means.

  return(model)
}
