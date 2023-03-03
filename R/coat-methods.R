#' @rdname coat
#' @method coef coat
#' @export
#' @usage NULL
#' @importFrom stats coef qnorm weighted.mean
#' @importFrom partykit data_party nodeapply nodeids
coef.coat <- function(object, node = NULL, drop = TRUE, ...) {
  if (is.null(node)) node <- nodeids(object, terminal = TRUE)
  cf <- if (inherits(object, "modelparty")) {
    nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  } else {
    lapply(node, function(n) {
      dat <- data_party(object, n)
      yn <- dat[["(response)"]]
      wn <- dat[["(weights)"]]
      if(is.null(wn)) wn <- rep.int(1, length(yn))
      mv <- c("(Intercept)" = weighted.mean(yn, wn))
      mv <- c(mv, "(Variance)" = weighted.mean((yn - mv)^2, wn))
    })
  }
  names(cf) <- node
  cf <- do.call(rbind, cf)
  if (drop) drop(cf) else cf
}

#' @rdname coat
#' @export
#' @usage NULL
#' @importFrom stats coef qnorm weighted.mean
#' @importFrom partykit id_node data_party info_node
#' @importFrom grid viewport gpar grid.clip grid.layout grid.lines grid.points grid.rect grid.text grid.xaxis grid.yaxis pushViewport popViewport upViewport unit
node_baplot <- function(obj,
                        col = 1,
                        linecol = 4,
		        lty = c(1, 2),
                        level = 0.95,
			bg = "white",
		        xscale = NULL,
		        yscale = NULL,
                        digits = 2,
		        ylines = 3,
			cex = 0.5,
		        id = TRUE,
                        mainlab = NULL, 
			gp = gpar())
{
    ## means and differences
    x <- obj$data[["means."]]
    y <- obj$data[["diffs."]]
    stopifnot(is.numeric(x), is.numeric(y))

    ## limits of agreement
    cf <- coef(obj, drop = FALSE)
    loa <- cbind(cf[, 1L] - qnorm((1 - level)/2) * sqrt(cf[, 2L]), cf[, 1L] + qnorm((1 - level)/2) * sqrt(cf[, 2L]))

    if (is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    if (is.null(yscale)) yscale <- range(c(y, loa)) + c(-0.1, 0.1) * diff(range(c(y, loa)))

    ### panel function for Bland-Altman plots in nodes
    rval <- function(node) {

        ## extract data
	nid <- id_node(node)
	dat <- data_party(obj, nid)
        xn <- dat[["means."]]
	yn <- dat[["diffs."]]
	wn <- dat[["(weights)"]]
	if(is.null(wn)) wn <- rep.int(1, length(yn))

        ## extract mean and variance
        cf <- info_node(node)$coefficients
        if(is.null(cf)) {
          cf <- c("(Intercept)" = weighted.mean(yn, wn))
          cf <- c(cf, "(Variance)" = weighted.mean((yn - cf)^2, wn))
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
	
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
	grid.clip()

        ## scatterplot
        grid.points(unit(xn, "native"), unit(yn, "native"), size = unit(cex, "char"), gp = gpar(col = col))

        ## limits of agreement
        loa <- cf[1L] + c(1, 0, -1) * qnorm((1 - level)/2) * sqrt(cf[2L])
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
