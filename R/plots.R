#  plot_coverage ---------------------------------------------------------------
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
#'
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

#  plot_pca_variance -----------------------------------------------------------
#' @title Plot variance explained by PCA components
#'
#' @description Plot the variance explained by each component. The number of 
#'   components shown is controlled by `num_components`, with up to the first
#'   10 componenets shown by default. If less than the requested number of 
#'   components exist, then all the components will be shown.
#'
#' @param pca output of \code{pca_mip_data}
#' @param num_components numeric for number of maximum components to be shown
#' 
#' @importFrom plotly plot_ly layout
#' @export
#' @examples
#' dat <- dummy_data()
#' dat <- filter_misc(dat = dat)
#' dat <- filter_coverage(dat = dat, min_coverage = 2)
#' dat <- melt_mip_data(dat = dat)
#' dat <- impute_mip_data(dat = dat)
#' pca <- pca_mip_data(dat = dat)
#' plot_pca_variance(pca, num_components = 3)

plot_pca_variance <- function(pca, num_components = 10) {

  # potential components
  nc <- min(length(colnames(pca$x)), num_components)

  # x has to be factored in order to preserve PC order
  x <- factor(colnames(pca$x)[seq_len(nc)],
              levels = colnames(pca$x)[seq_len(nc)])

  out <- plot_ly(x = x, y = pca$var[seq_len(nc)],
                 type = "bar", width = 500, height = 500) %>%
    plotly::layout(title = "Percentage of total variance explaned by PCA",
           xaxis = list(title = "Principal Component"),
           yaxis = list(title = "Percentage of variance explained"))

  # return invisibly
  print(out)
  invisible(out)
}

#  plot_pca -----------------------------------------------------------
#' @title Plot PCA
#'
#' @description Plots either the first 2 or 3 principal components, with the
#'   data points labelled accordng to chosen meta data in the mip data set.
#'   
#' @details Using the output of \code{pca_mip_data} and a specified variable
#'   within the mip data set, e.g the country of sample collection, a 
#'   scatterplot of the data is produced. Either the first 2 or 3 components
#'   can be used, as specified with `num_components`. The chosen `meta_var` 
#'   must be a variable found within the mip data set (this is included as 
#'   the last element in the object returned by \code{pca_mip_data}).
#'
#' @param pca output of \code{pca_mip_data}
#' @param num_components numeric for number of components used. Default = 2
#' @param meta_var character for the desired meta variable to be used for 
#'   labelling the scatterplot. Default = "Country".
#' 
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' dat <- dummy_data()
#' dat <- filter_misc(dat = dat)
#' dat <- filter_coverage(dat = dat, min_coverage = 2)
#' dat <- melt_mip_data(dat = dat)
#' dat <- impute_mip_data(dat = dat)
#' pca <- pca_mip_data(dat = dat)
#' plot_pca(pca, num_components = 2, meta_var = "Country")
#' plot_pca(pca, num_components = 3, meta_var = "Study")

plot_pca <- function(pca, num_components = 2, meta_var = "Country") {

  # check meta_var
  if (!is.element(meta_var, colnames(pca$dat))) {
    stop(meta_var, " is not a member of pca$dat")
  }

  # check num_components
  nc <- min(length(colnames(pca$x)), num_components)
  if (nc != num_components) {
    message("Using all the components (", nc, ") that are available.")
  }

  # create color vector
  col_vec <- suppressWarnings(
    brewer.pal(length(unique(pca$dat[, meta_var])), "Set1")
    )

  if (nc == 2) {
  # scatterplot of first 2 principal components
  # 2D scatter
  out <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
          color = pca$dat[, meta_var], type = "scatter", colors = col_vec,
          mode = "markers", marker = list(size = 5))

  }

  if (nc == 3) {
  # scatterplot of first 3 principal components
  # 3D scatter
  out <- plot_ly(as.data.frame(pca$x[, 1:3]), x = ~PC1, y = ~PC2, z = ~PC3,
          color = pca$dat[, meta_var], type = "scatter3d", colors = col_vec,
          mode = "markers", marker = list(size = 3))

  }

  # return invisibly
  print(out)
  invisible(out)
}
