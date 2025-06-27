#' Methods for Conditional Method Agreement Trees (COAT)
#'
#' Extracting information from or visualization of conditional method agreement trees.
#' Visualizations use trees with Bland-Altman plots in terminal nodes, drawn
#' via grid graphics.
#'
#' Various methods are provided for trees fitted by \code{\link[coat]{coat}},
#' in particular \code{print}, \code{plot} (via \pkg{grid}/\pkg{partykit}) and
#' \code{coef}. The \code{plot} method draws Bland-Altman plots in the terminal panels by default,
#' using the function \code{node_baplot}.
#'
#' In addition to these dedicated \code{coat} methods, further methods are inherited
#' from \code{\link[partykit]{ctree}} or \code{\link[partykit]{mob}}, respectively,
#' depending on which \code{type} of \code{coat} was fitted.
#'
#' @param x,object,obj a \code{coat} object as returned by \code{\link[coat]{coat}}.
#' @param digits numeric. Number of digits used for rounding the displayed coefficients
#' or limits of agreement.
#' @param header,footer logical. Should a header/footer be printed for the tree?
#' @param title character with the title for the tree.
#' @param node integer. ID of the node for which the Bland-Altman parameters
#' (coefficients) should be extracted.
#' @param drop logical. Should the matrix attribute be dropped if the parameters
#' from only a single node are extracted?
#' @param terminal_panel a panel function or panel-generating function passed to
#' \code{\link[partykit]{plot.party}}. By default, \code{node_baplot} is used to
#' generate a suitable panel function for drawing Bland-Altman plots based on the
#' the provided \code{coat} object. It can be customized using the \code{tp_args} argument
#' (passed through \code{...}).
#' @param tnex numeric specification of the terminal node extension
#' relative to the inner nodes (default is twice the size).
#' @param drop_terminal logical. Should all terminal nodes be "dropped" to
#' the bottom row?
#' @param level numeric level for the limits of agreement.
#' @param pch,cex,col,linecol,lty,bg graphical parameters for the scatter plot and limits
#' of agreement in the Bland-Altman plot (scatter plot character, character extension, plot color,
#' line color, line types, and background color).
#' @param confint logical. Should nonparametric bootstrap percentile confidence intervals be plotted?
#' @param B numeric. Number of bootstrap samples to be used if \code{confint = TRUE}.
#' @param cilevel numeric. Level of the confidence intervals if \code{confint = TRUE}.
#' @param cicol color specification for the confidence intervals if \code{confint = TRUE}.
#' @param xscale,yscale numeric specification of scale of x-axis and y-axis, respectively.
#' By default the range of all scatter plots and limits of agreement across all nodes
#' are used.
#' @param ylines numeric. Number of lines for spaces in y-direction.
#' @param id logical. Should node IDs be plotted?
#' @param mainlab character or function. An optional title for the plots. Either
#' a character or a \code{function(id, nobs)}.
#' @param gp grid graphical parameters.
#' @param ... further arguments passed to methods.
#'
#' @return The \code{print()} method returns the printed object invisibly.
#' The \code{coef()} method returns the vector (for a single node) or matrix (for multiple nodes) of estimated parameters (bias and standard deviation).
#' The \code{plot()} method returns \code{NULL}.
#' The \code{node_baplot()} panel-generating function returns a function that can be plugged into the \code{plot()} method.
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
#' ## conditional method agreement tree
#' tr <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide)
#'
#' ## illustration of methods (including some customization)
#'
#' ## printing
#' print(tr)
#' print(tr, header = FALSE, footer = FALSE)
#'
#' ## extracting Bland-Altman parameters
#' coef(tr)
#' coef(tr, node = 1)
#'
#' ## visualization (via grid with node_baplot)
#' plot(tr)
#' plot(tr, ip_args = list(id = FALSE),
#'   tp_args = list(col = "slategray", id = FALSE, digits = 3, pch = 19))
#'

#' @rdname coat-methods
#' @method print coat
#' @export
#' @importFrom partykit nodeapply nodeids print.party width
#' @importFrom stats weighted.mean
print.coat <- function(x, digits = 2L,
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
    s <- sqrt(weighted.mean((y - m)^2, w) * sum(w)/(sum(w) - 1))
    paste(c("Bias =", "SD ="), format(round(c(m, s), digits = digits), nsmall = digits), collapse = ", ")
  }, by_node = FALSE)

  terminal_panel <- function(node) paste(":", node_labs[[id_node(node)]])

  print.party(x, terminal_panel = terminal_panel, header_panel = header_panel, footer_panel = footer_panel, ...)
  invisible(x)
}


#' @rdname coat-methods
#' @method coef coat
#' @export
#' @importFrom stats coef weighted.mean
#' @importFrom partykit data_party nodeapply nodeids
coef.coat <- function(object, node = NULL, drop = TRUE, ...) {
  if (is.null(node)) node <- nodeids(object, terminal = TRUE)
  cf <- if (inherits(object, "modelparty")) {
    nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  } else {
    lapply(node, function(n) {
      dat <- data_party(object, n)
      yn <- dat[["(response)"]]
      yn <- yn[, 1L] - yn[, 2L]
      wn <- dat[["(weights)"]]
      if(is.null(wn)) wn <- rep.int(1, length(yn))
      mv <- c("Bias" = weighted.mean(yn, wn))
      mv <- c(mv, "SD" = sqrt(weighted.mean((yn - mv)^2, wn) * sum(wn)/(sum(wn) - 1)))
    })
  }
  names(cf) <- node
  cf <- do.call(rbind, cf)
  if (drop) drop(cf) else cf
}


#' @rdname coat-methods
#' @method plot coat
#' @export
#' @importFrom partykit plot.party
plot.coat <- function(x, terminal_panel = node_baplot, tnex = 2, drop_terminal = TRUE, ...) {
  partykit::plot.party(x, terminal_panel = terminal_panel, tnex = tnex, drop_terminal = drop_terminal, ...)
}


#' @rdname coat-methods
#' @export
#' @importFrom stats coef qnorm weighted.mean quantile
#' @importFrom partykit id_node data_party info_node
#' @importFrom grid viewport gpar grid.clip grid.layout grid.lines grid.points grid.polygon grid.rect grid.text grid.xaxis grid.yaxis pushViewport popViewport upViewport unit
node_baplot <- function(obj,
                        level = 0.95,
                        digits = 2,
                        pch = 1,
			cex = 0.5,
                        col = 1,
                        linecol = 4,
		        lty = c(1, 2),
			bg = "white",
			confint = FALSE,
			B = 500,
			cilevel = 0.95,
                        cicol = "lightgray",
		        xscale = NULL,
		        yscale = NULL,
		        ylines = 3,
		        id = TRUE,
                        mainlab = NULL,
			gp = gpar())
{
    ## means and differences
    y <- obj$fitted[["(response)"]]
    x <- (y[, 1L] + y[, 2L])/2
    y <- y[, 1L] - y[, 2L]
    stopifnot(is.numeric(x), is.numeric(y))

    ## ## limits of agreement
    cf <- coef(obj, drop = FALSE)
    loa <- cbind(cf[, 1L] - qnorm((1 - level)/2) * cf[, 2L], cf[, 1L] + qnorm((1 - level)/2) * cf[, 2L])

    if (is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    if (is.null(yscale)) yscale <- range(c(y, loa)) + c(-0.1, 0.1) * diff(range(c(y, loa)))

    ### panel function for Bland-Altman plots in nodes
    rval <- function(node) {

        ## extract data
	nid <- id_node(node)
	dat <- data_party(obj, nid)
        yn <- dat[["(response)"]]
        xn <- (yn[, 1L] + yn[, 2L])/2
	yn <- yn[, 1L] - yn[, 2L]
	wn <- dat[["(weights)"]]
	if(is.null(wn)) wn <- rep.int(1, length(yn))

        ## extract mean and variance
        cf <- info_node(node)$coefficients
        if(is.null(cf)) {
          cf <- c("Bias" = weighted.mean(yn, wn))
          cf <- c(cf, "SD" = sqrt(weighted.mean((yn - cf)^2, wn) * sum(wn)/(sum(wn) - 1)))
        }

        ## grid
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1),
                                         c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"),
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_baplot", nid, sep = ""),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], sum(wn))
	}
        grid.text(mainlab)
        popViewport()

        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = xscale, yscale = yscale,
			 name = paste0("node_baplot", nid, "plot"),
			 clip = FALSE)

        pushViewport(plot)

	## confidence intervals
	if (confint) {
	  loa_boot <- sapply(1:B, function(z) {
	    boot_index <- sample(1:length(yn), length(yn), replace = TRUE)
	    wn_boot <- weighted.mean(yn[boot_index], wn[boot_index])
	    sd_boot <- sqrt(weighted.mean((yn[boot_index] - wn_boot)^2, wn[boot_index]) * sum(wn[boot_index])/(sum(wn[boot_index]) - 1))

	    wn_boot + c(1, 0, -1) * qnorm((1 - level)/2) * sd_boot
	  })

	  stats_boot <- apply(loa_boot, 1, function(z) quantile(z, probs = 0:1 + c(1, -1) * (1-cilevel)/2))

	  grid.polygon(unit(c(0, 1, 1, 0), "npc"), unit(rep(stats_boot[, 1L], each = 2L), "native"), gp = gpar(col = cicol, fill = cicol))
	  grid.polygon(unit(c(0, 1, 1, 0), "npc"), unit(rep(stats_boot[, 2L], each = 2L), "native"), gp = gpar(col = cicol, fill = cicol))
          grid.polygon(unit(c(0, 1, 1, 0), "npc"), unit(rep(stats_boot[, 3L], each = 2L), "native"), gp = gpar(col = cicol, fill = cicol))
	}

        ## box and axes
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
	grid.clip()

	## scatterplot
        grid.points(unit(xn, "native"), unit(yn, "native"), size = unit(cex, "char"), pch = pch, gp = gpar(col = col))

        ## limits of agreement
        loa <- cf[1L] + c(1, 0, -1) * qnorm((1 - level)/2) * cf[2L]
        grid.lines(unit(c(0, 1), "npc"), unit(loa[2L], "native"), gp = gpar(col = linecol, lty = lty[1L]))
        grid.lines(unit(c(0, 1), "npc"), unit(loa[1L], "native"), gp = gpar(col = linecol, lty = lty[2L]))
        grid.lines(unit(c(0, 1), "npc"), unit(loa[3L], "native"), gp = gpar(col = linecol, lty = lty[2L]))

        ## annotation
        if (isTRUE(digits)) digits <- 2L
        if (is.numeric(digits)) {
          loalab <- format(round(loa, digits = digits), nsmall = digits)
          for (i in 1L:3L) {
            grid.rect(
              x = unit(1, "npc") - unit(1, "lines") - max(unit(0.5, "strwidth", loalab)),
              y = unit(loa[i], "native"),
              width = unit(1, "lines") + max(unit(1, "strwidth", loalab)),
              height = unit(1, "lines") + max(unit(1, "strheight", loalab)),
              gp = gpar(col = linecol, fill = bg))
            grid.text(loalab[i],
              x = unit(1, "npc") - unit(1, "lines") - max(unit(0.5, "strwidth", loalab)),
              y = unit(loa[i], "native"),
              gp = gpar(col = linecol))
          }
        }

        upViewport(2)
    }

    return(rval)
}
class(node_baplot) <- "grapcon_generator"
