#' Conditional method agreement trees (coat)
#'
#' Functions to fit and plot \code{coat} models.
#'
#' @param y1 a character string specifying the variable in data containing the measurements by one method.
#' @param y2 a character string specifying the variable in data containing the measurements by another method.
#' @param covars a character string or vector of a single or multiple covariates.
#' @param data a data frame containing \code{y1}, \code{y2} and \code{covars}.
#' @param means a logical indicating whether intraindividual mean values of measurements shall be included as covariate.
#' @param type a character string specifying the type of coat model to be fit. Either \code{"ctree"} (default), \code{"disttree"} (equals \code{"ctree"}) or \code{"mob"}.
#' @param x an object as returned by \code{\link[coat]{coat}}.
#' @param digits a numeric specifying the number of digits to display.
#' @param xlim.max an optional numeric value to define the upper limit of the x-axis.
#' @param level a numeric specifying the desired coverage of the prediction interval.
#' @param ... further arguments passed to \code{\link[partykit]{ctree_control}} or the fit function of \code{\link[partykit]{mob}}.
#'
#' @details The minimum number of observations required to model conditional agreement defaults to 20. Users may choose to modify this value as needed. See \code{\link[partykit]{ctree_control}} and \code{\link[partykit]{mob_control}} for details.
#'
#' @examples
#' \dontshow{ if(!requireNamespace("MethComp")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("the MethComp package is required for this example but is not installed")
#'   } else q() }
#' }
#' ### load package ###
#' library("coat")
#'
#' ### data ###
#' data("scint", package = "MethComp")
#' ## transform data to required 'wide' format
#' scint_wide <- reshape(scint, v.names = "y", timevar = "meth", idvar = "item", direction = "wide")
#'
#' ### fit coat model using ctree() ###
#' mytree1 <- coat("y.DTPA", "y.DMSA", c("age", "sex"), data = scint_wide)
#' ## including mean values as predictor
#' mytree2 <- coat("y.DTPA", "y.DMSA", c("age", "sex"), data = scint_wide, means = TRUE)
#'
#' ### plot ###
#' plot(mytree1)
#' plot(mytree1, digits = 2, xlim.max = 120)
#' plot(mytree2, digits = 2, xlim.max = 120)
#'
#' @return Object of class \code{coat}.
#'
#' @importFrom stats as.formula dnorm lm lm.fit qnorm sd
#' @importFrom ggplot2 aes geom_hline geom_label geom_point theme_bw xlab ylab xlim
#' @importFrom gridExtra arrangeGrob grid.arrange tableGrob ttheme_minimal
#' @importFrom ggparty ggparty geom_edge geom_edge_label geom_node_label geom_node_plot geom_node_splitvar
#' @importFrom grid textGrob gpar
#'
#' @export
coat <- function(y1, y2, covars, data, means = FALSE, type = c("ctree", "disttree", "mob")[1], ...) {

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

  # fit the tree model
  if (type == "mob") {
    model <- partykit::mob(model.formula, data = moddat, fit = gaussfit, ...)
  } else {
    model <- partykit::ctree(model.formula, data = moddat, ytrafo = meanvar, ...)
  }

  # Update data with means
  model$data$means. <- moddat$means.

  # define class of returned object
  class(model) <- c("coat", class(model))

  return(model)
}

#' @describeIn coat function to print a coat model.
#' @export
print.coat <- function(x, digits = 2, ...) {
  if (inherits(x, "modelparty")) {
    x <- partykit::as.constparty(x)
  }
  
  partykit:::print.constparty(x, FUN = function(y1, w, digits) paste(c("Bias =", "SD ="), round(c(mean(y1), sd(y1)), digits), collapse = ", "), ...)
}
                              
#' @describeIn coat function to plot a coat model.
#' @export
plot.coat <- function(x, digits = 2, xlim.max = NULL, level = 0.95, ...) {
  diffs. <- id <- means. <- nodesize <-  p.value <- splitvar <- NULL # due to NSE notes in R CMD check

  if (is.null(xlim.max)) xlim.max <- max(x$data$means.)
  level <- 1 - (1 - level)/2

  p1 <- ggparty(x, terminal_space = 0.5)
  mean_diff <- sapply(p1$data$nodedata_diffs., mean)
  sd_diff <- sapply(p1$data$nodedata_diffs., sd)

  p1 + geom_edge() +
    geom_edge_label() +
    geom_node_splitvar() +
    geom_node_plot(gglist = list(aes(x = means., y = diffs.),
                                 geom_point(alpha = 0.8),
                                 geom_hline(aes(yintercept = mean_diff[id]), col = "blue"),
                                 geom_hline(aes(yintercept = mean_diff[id] + qnorm(level)*sd_diff[id]), col = "blue", linetype = "dashed"),
                                 geom_hline(aes(yintercept = mean_diff[id] - qnorm(level)*sd_diff[id]), col = "blue", linetype = "dashed"),
                                 geom_label(aes(x = xlim.max * 0.95, y = mean_diff[id], label = round(mean_diff[id], digits)), col = "blue"),
                                 geom_label(aes(x = xlim.max * 0.95, y = mean_diff[id] + qnorm(level)*sd_diff[id], label = round(mean_diff[id] + qnorm(level)*sd_diff[id], digits)), col = "blue"),
                                 geom_label(aes(x = xlim.max * 0.95, y = mean_diff[id] - qnorm(level)*sd_diff[id], label = round(mean_diff[id] - qnorm(level)*sd_diff[id], digits)), col = "blue"),
                                 theme_bw(base_size = 10), xlab("Mean values"), ylab("Differences"),
                                 xlim(NA, xlim.max))) +

      geom_node_label(line_list = list(aes(label = splitvar),
                                       aes(label = paste("p = ", round(p.value, 3)))),
                      line_gpar = list(list(size = 12, col = "black"),
                                       list(size = 11)),
                      ids = "inner") +
      geom_node_label(aes(label = paste0("Node ", id, ", N = ", nodesize)),
                      size = 4, nudge_x = 0.02, nudge_y = 0.01,
                      ids = "terminal")
}


