
#------------------------------------------------
#' @title Plot coverage gaps in data
#'
#' @description Plot coverage gaps in data
#'
#' @details TODO
#'
#' @param dat TODO
#' @param type TODO
#' @param pch TODO
#' @param ylim TODO
#' @param xlab TODO
#' @param ylab TODO
#' @param log TODO
#' @param ... other arguments to \code{plot()}
#'
#' @export
#' @examples
#' # TODO

plot_coverage <- function(dat, type = "o", pch = 20, ylim = c(0,100), xlab = "Coverage (unique reads)", ylab = "data completeness (%)", log="x", ...) {

  # get dropout at different levels of coverage filtering
  x <- c(1:10, seq(20,200,10))
  y <- mapply(function(x) {mean(dat$Coverage >= x)}, x) * 100

  plot(x, y, type = type, log = log, pch = pch, ylim = ylim, xlab = xlab, ylab = ylab, ...)

  # return invisibly
  invisible(dat)
}
