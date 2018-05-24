
#------------------------------------------------
#' @title Plot coverage gaps in data
#'
#' @description Plot coverage gaps in data
#'
#' @details Plot a line graph showing the percentage of sites that have more 
#'   than a given number of unique reads (coverage). 
#'
#' @param dat MIP data. The data must have the following variables:
#'   \itemize{
#'       \item Coverage : The total read coverage as numerics
#'       }
#' @param ... other arguments to \code{plot}
#' @param pch Either an integer or a single character specifying a symbol to be 
#'   used as the default in plotting points. See \code{\link[graphics]{points}} 
#'   for possible values and their interpretation. Note that only integers and 
#'   single-character strings can be set as a graphics parameter 
#'   (and not NA nor NULL).
#' @importFrom graphics plot par
#' @inheritParams graphics::plot.default
#' @inheritParams graphics::points
#' @export
#' @examples
#' dat <- dummy_data()
#' dat2 <- filter_misc(dat)
#' plot_coverage(dat2)

plot_coverage <- function(dat, type = "o", pch = 20, ylim = c(0, 100),
                          xlab = "Coverage (unique reads)",
                          ylab = "data completeness (%)",
                          log = "x",
                          ...) {

  # get dropout at different levels of coverage filtering
  x <- c(1:10, seq(20, 200, 10))
  cover_func <- function(x) {
    mean(dat$Coverage >= x)
  }
  y <- mapply(cover_func, x) * 100

  plot(x, y, type = type, log = log, pch = pch, ylim = ylim,
       xlab = xlab, ylab = ylab, ...)

  # return invisibly
  invisible(dat)
}
