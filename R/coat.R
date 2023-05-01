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
#' @param label.align a numeric between 0 and 1 specifying the alignment of labels relative to the plot width or \code{xlim.max}.
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
#' mytree1 <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide)
#' ## including mean values as predictor
#' mytree2 <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide, means = TRUE)
#'
#' ### plot ###
#' plot(mytree1)
#' plot(mytree1, digits = 2, xlim.max = 120)
#' plot(mytree2, digits = 2, xlim.max = 120)
#'
#' @return Object of class \code{coat}.
#'
#' @importFrom stats as.formula dnorm lm lm.fit qnorm sd
#' @importFrom ggplot2 ggplot aes geom_hline geom_label geom_point theme_bw xlab ylab xlim theme ggtitle margin
#' @importFrom gridExtra arrangeGrob grid.arrange tableGrob ttheme_minimal
#' @importFrom ggtext element_markdown
#' @importFrom ggparty ggparty geom_edge geom_edge_label geom_node_label geom_node_plot geom_node_splitvar
#' @importFrom grid textGrob gpar
#'
#' @export
coat <- function(formula, data, subset, na.action, weights, means = FALSE, type = c("ctree", "disttree", "mob"), ...)
{
  ## type of tree
  type <- match.arg(tolower(type), c("ctree", "disttree", "mob"))

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## preparation of ctree/mob call
  m <- match.call(expand.dots = FALSE)
  if(missing(na.action)) m$na.action <- na.omit

  ## if desired "means(y1, y2)" is added as split variable
  if(means) {
    formula <- update(formula, . ~ `_means_` + .)
    formula[[3L]][[2L]] <- formula[[2]]
    formula[[3L]][[2L]][[1L]] <- as.name("means")
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

#' @describeIn coat function to print a coat model.
#' @export
print.coat <- function(x, FUN = NULL, digits = 2L,
  header = TRUE, footer = TRUE, title = "Conditional method agreement tree (COAT)", ...)
{
  header_panel <- if(header) function(party) {      
    c(title, "", "Model formula:", deparse(party$info$call$formula), "", "Fitted party:", "")
  } else function(party) ""
  
  footer_panel <- if(footer) function(party) {
    n <- width(party)
    n <- format(c(length(party) - n, n))
    c("", paste("Number of inner nodes:   ", n[1L]),
      paste("Number of terminal nodes:", n[2L]), "")
  } else function (party) ""

  node_labs <- nodeapply(x, nodeids(x), function(node) {
    y <- node$fitted[["(response)"]]
    y <- y[, 1L] - y[, 2L]
    w <- node$fitted[["(weights)"]]
    if (is.null(w)) w <- rep.int(1, NROW(y))
    m <- weighted.mean(y, w)
    s <- sqrt(weighted.mean((y - m)^2, w))
    paste(c("Bias =", "SD ="), format(round(c(m, s), digits = digits), nsmall = digits), collapse = ", ")
  }, by_node = FALSE)

  terminal_panel <- function(node) paste(":", node_labs[[id_node(node)]])

  print.party(x, terminal_panel = terminal_panel, header_panel = header_panel, footer_panel = footer_panel, ...)
  invisible(x)
}
                              
#' @describeIn coat function to plot a coat model.
#' @export
plot.coat <- function(x, terminal_panel = node_baplot, tnex = 2, drop_terminal = TRUE, ...) {
  partykit::plot.party(x, terminal_panel = terminal_panel, tnex = tnex, drop_terminal = drop_terminal, ...)
}

autoplot.coat <- function(x, digits = 2, xlim.max = NULL, level = 0.95, label.align = 0.95, ...) {
  diffs. <- id <- means. <- nodesize <-  p.value <- splitvar <- NULL # due to NSE notes in R CMD check
  
  ## augment data
  y <- x$fitted[["(response)"]]
  x$data$means. <- (y[, 1L] + y[, 2L])/2
  x$data$diffs. <- y[, 1L] - y[, 2L]
  
  if (is.null(xlim.max)) xlim.max <- max(x$data$means.)
  level <- 1 - (1 - level)/2
  
  if (length(x) == 1) {
    p1 <- ggplot(data_party(x), aes(x = means., y = diffs.))
    mean_diff <- mean(p1$data$diffs.)
    sd_diff <- sd(p1$data$diffs.)
    
    p2 <- p1 + geom_point(alpha = 0.8) + 
      geom_hline(aes(yintercept = mean_diff), col = 4) + 
      geom_hline(aes(yintercept = mean_diff + qnorm(level)*sd_diff), col = 4, linetype = "dashed") + 
      geom_hline(aes(yintercept = mean_diff - qnorm(level)*sd_diff), col = 4, linetype = "dashed") + 
      geom_label(aes(x = xlim.max * label.align, y = mean_diff, label = round(mean_diff, digits)), col = 4) +
      geom_label(aes(x = xlim.max * label.align, y = mean_diff + qnorm(level)*sd_diff, label = round(mean_diff + qnorm(level)*sd_diff, digits)), col = 4) + 
      geom_label(aes(x = xlim.max * label.align, y = mean_diff - qnorm(level)*sd_diff, label = round(mean_diff - qnorm(level)*sd_diff, digits)), col = 4) + 
      theme_bw(base_size = 10) + xlab("Mean values") + ylab("Differences") + xlim(NA, xlim.max)
    p2 + ggtitle(paste0("Node ", 1, ", N = ", length(p1$data$diffs.))) + 
      theme(plot.title = element_markdown(hjust = 0.5, linetype = 1, padding = unit(0.25, "lines"), r = grid::unit(0.15, "lines"), margin = margin(0, 0, 0, 0)))
  } else {
    p1 <- ggparty(x, terminal_space = 0.5)
    mean_diff <- sapply(p1$data$nodedata_diffs., mean)
    sd_diff <- sapply(p1$data$nodedata_diffs., sd)
    
    p1 + geom_edge() +
      geom_edge_label() +
      geom_node_splitvar() +
      geom_node_plot(gglist = list(aes(x = means., y = diffs.),
                                   geom_point(alpha = 0.8),
                                   geom_hline(aes(yintercept = mean_diff[id]), col = 4),
                                   geom_hline(aes(yintercept = mean_diff[id] + qnorm(level)*sd_diff[id]), col = 4, linetype = "dashed"),
                                   geom_hline(aes(yintercept = mean_diff[id] - qnorm(level)*sd_diff[id]), col = 4, linetype = "dashed"),
                                   geom_label(aes(x = xlim.max * label.align, y = mean_diff[id], label = round(mean_diff[id], digits)), col = 4),
                                   geom_label(aes(x = xlim.max * label.align, y = mean_diff[id] + qnorm(level)*sd_diff[id], label = round(mean_diff[id] + qnorm(level)*sd_diff[id], digits)), col = 4),
                                   geom_label(aes(x = xlim.max * label.align, y = mean_diff[id] - qnorm(level)*sd_diff[id], label = round(mean_diff[id] - qnorm(level)*sd_diff[id], digits)), col = 4),
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
}
