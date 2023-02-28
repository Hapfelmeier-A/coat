#' Two-sample Bland-Altman test of method agreement
#'
#' Function to perform a two-sample Bland-Altman test for hypothesis testing of differences in method agreement. Additional functions are given for printing and plotting.
#'
#' @param y1 a numeric vector of measurements by one method.
#' @param y2 a numeric vector of measurements by another method.
#' @param groups a factor with two levels indicating the two independent groups or samples to be compared.
#' @param x an object as returned by \code{\link[COAT]{BAtest}}.
#' @param statistics a logical indicating whether plotting of test statistics is required.
#' @param digits a numeric specifying the number of digits to display.
#' @param xlim.max an optional numeric value to define the upper limit of the x-axis.
#' @param level a numeric specifying the desired coverage of the prediction interval.
#' @param ... further arguments passed to \code{\link[partykit]{ctree_control}}.
#'
#' @examples
#' ### load package ###
#' library("COAT")
#'
#' ### data ###
#' require(MethComp)
#' data(VitCap)
#' ## transform data to required 'wide' format
#' VitCap_wide <- reshape(VitCap, v.names = "y", timevar = "instrument",
#'                        idvar = c("item", "user"), drop = "meth", direction = "wide")
#'
#' ### perform a two-sample BA-test ###
#' testresult <- BAtest(VitCap_wide$y.St, VitCap_wide$y.Exp, VitCap_wide$user)
#' testresult
#'
#' ### plot ###
#' plot(testresult)
#' ## expand xlim to avoid overlapping labeling
#' plot(testresult, digits = 2, xlim.max = 4800)
#' ## plot without table of statistics
#' plot(testresult, digits = 2, statistics = FALSE, xlim.max = 4800)
#'
#' @return Object of class \code{BAtest} with elements
#' \item{\code{test}}{result of the Bland-Altman test.}
#' \item{\code{model}}{tree model used to perform the Bland-Altman test.}
#'
#' @importFrom strucchange sctest
#' @importFrom ggplot2 aes geom_hline geom_label geom_point theme_bw xlab ylab xlim
#' @importFrom gridExtra arrangeGrob grid.arrange tableGrob ttheme_minimal
#' @importFrom ggparty ggparty geom_edge geom_edge_label geom_node_label geom_node_plot geom_node_splitvar
#' @importFrom grid textGrob gpar
#'
#' @export
BAtest <- function(y1, y2, groups, ...) {

  # coerce predictor to factor variable
  if (!is.factor(groups)) {
    groups <- as.factor(groups)
    warning(paste(groups, "has been coerced to a factor."))
  }

  # levels of predictor
  lvs <- levels(groups)

  # check whether there are only two groups
  if (length(lvs) != 2) {
    stop(paste(groups, "is required to have exactly two levels."))
  }

  # check whether samples are at least of size 3.
  if (any(table(groups) < 3)) {
    stop("Each sample is required to contain a minimum of three observations.")
  }

  # check whether y1 and y2 are numeric
  if (!(inherits(y1, c("numeric", "integer")) & inherits(y2, c("numeric", "integer")))) {
    stop(paste(y1, "and", y2, "need to be of class 'numeric' or 'integer'."))
  }

  # calculate differences
  diffs. <- y1 - y2

  # create model data
  moddat <- data.frame(y1, y2, groups, diffs., "means." = rowMeans(cbind(y1, y2)))

  # calculate descriptive statistics
  descriptives <- by(diffs., groups, function(x) c(mean(x), sd(x)))

  # fit the tree model
  tree_mean <- partykit::ctree(diffs. ~ groups, data = moddat, alpha = 1, minsplit = 6L, minbucket = 3L, ...)
  tree_var <- partykit::ctree(diffs. ~ groups, data = moddat, alpha = 1, ytrafo = .var, minsplit = 6L, minbucket = 3L, ...)
  tree_meanvar <- partykit::ctree(diffs. ~ groups, data = moddat, alpha = 1, ytrafo = meanvar, minsplit = 6L, minbucket = 3L, ...)

  # matrix of test statistics
  test <- matrix(NA, nrow = 3, ncol = 5)
  colnames(test) <- c(lvs, "Chisq", "df", "p-value")
  rownames(test) <- c("Bias", "SD", "Total")

  test[1, c(3, 5)] <- sctest(tree_mean)[[1]]
  test[2, c(3, 5)] <- sctest(tree_var)[[1]]
  test[3, c(3, 5)] <- sctest(tree_meanvar)[[1]]

  test[1:2, 1] <- descriptives[[1]]
  test[1:2, 2] <- descriptives[[2]]

  test[, 4] <- c(1, 1, 2)

  # Update data with means
  tree_meanvar$data$means. <- moddat$means.

  # list of output
  output <- list("test" = test, "model" = tree_meanvar)

  # define class of returned object
  class(output) <- "BAtest"

  return(output)
}

#' @describeIn BAtest function to print the result of the Bland-Altman test.
#' @export
print.BAtest <- function(x, digits = 2, ...) {

  out <- x$test
  out[, 1] <- round(out[, 1], digits)
  out[, 2] <- round(out[, 2], digits)
  out[, 3] <- round(out[, 3], 3)
  out[, 5] <- round(out[, 5], 3)

  print(out, na.print="", ...)
}

#' @describeIn BAtest function to plot the result of the Bland-Altman test.
#' @export
plot.BAtest <- function(x, statistics = TRUE, digits = 2, xlim.max = NULL, level = 0.95, ...) {
  diffs. <- id <- means. <- nodesize <-  p.value <- NULL # due to NSE notes in R CMD check

  if (is.null(xlim.max)) xlim.max <- max(x$model$data$means.)
  level <- 1 - (1 - level)/2

  p1 <- ggparty(x$model, terminal_space = 0.5)
  mean_diff <- sapply(p1$data$nodedata_diffs., mean)
  sd_diff <- sapply(p1$data$nodedata_diffs., sd)

  p1 <- p1 + geom_edge() +
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
                                 xlim(NA, xlim.max)),
                   shared_axis_labels = TRUE, legend_separator = TRUE) +

    geom_node_label(line_list = list(aes(label = "p-value:"),
                                     aes(label = round(p.value, 3))),
                    line_gpar = list(list(size = 12, col = "black"),
                                     list(size = 11)),
                    ids = "inner") +
    geom_node_label(aes(label = paste0("N = ", nodesize)),
                    size = 4, nudge_x = 0.02, nudge_y = 0.01,
                    ids = "terminal")

  if (statistics) {
    out <- x$test
    out[, 1] <- round(out[, 1], digits)
    out[, 2] <- round(out[, 2], digits)
    out[, 3] <- round(out[, 3], 3)
    out[, 5] <- round(out[, 5], 3)
    out[3, 1:2] <- ""
    p2 <- tableGrob(out, theme = ttheme_minimal())
    grid.arrange(p1, arrangeGrob(p2, top = textGrob("Bland-Altman test", gp = gpar(fontsize = 14, font = 3))), ncol = 1, nrow = 2, heights = c(3, 1))
  } else p1
}
