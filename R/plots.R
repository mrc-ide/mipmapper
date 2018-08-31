
#------------------------------------------------
# generate any number of divergent colours
#' @noRd
divergent_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- colorRampPalette(raw_cols)
  
  # simple case if small K
  if (K<=2) {
    return(my_palette(K))
  }
  
  # some logic to choose a palette size and sequence of colours that is
  # consistent across different values of K
  ncol <- 3
  while(ncol<K) {
    ncol <- ncol+(ncol-1)
  }
  dist_mat <- matrix(1:ncol, ncol, ncol)
  dist_mat <- abs(t(dist_mat)-dist_mat)
  x <- rep(FALSE, ncol)
  
  col_index <- 1
  for (i in 2:K) {
    x[col_index] <- TRUE
    s <- apply(dist_mat[which(x),,drop=FALSE], 2, min)
    next_index <- which.max(s)
    col_index <- c(col_index, next_index)
  }
  col_index
  ret <- my_palette(ncol)[col_index]
  
  return(ret)
}

#------------------------------------------------
#' @title Plot raster image of coverage matrix
#'
#' @description Plot raster image of coverage matrix.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#'
#' @export

plot_coverage_matrix <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (is.null(project$data_processed$data_coverage)) {
    stop("no data loaded")
  }
  
  # plot raster image
  mat_raster <- raster(project$data_processed$data_coverage)
  raster::plot(mat_raster, col = plasma(100), axes = FALSE, box = FALSE,
               xlab = "locus", ylab = "sample", legend.args = list(text = "Coverage"))
}

#------------------------------------------------
#' @title Plot raster image of within-sample allele frequencies
#'
#' @description Plot raster image of within-sample allele frequencies.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#'
#' @export

plot_WSAF_matrix <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (is.null(project$data_processed$data_coverage)) {
    stop("no data loaded")
  }
  
  # plot raster image
  mat_raster <- raster(project$data_processed$WSAF)
  raster::plot(mat_raster, col = plasma(100), axes = FALSE, box = FALSE,
               xlab = "locus", ylab = "sample", legend.args = list(text = "WSAF"))
  
}

#------------------------------------------------
#' @title Plot PCA
#'
#' @description Plots either the first 2 or 3 principal components, with the
#'   data points labelled accordng to chosen meta data.
#'
#' @param project a \code{mipmapper_project} containing pca analysis output
#' @param num_components numeric for number of components used
#' @param grouping_var character for the desired meta variable to be used for 
#'   labelling the scatterplot
#'
#' @export

plot_pca <- function(project, num_components = 2, grouping_var = NULL) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$analyses$pca) == 0) {
    stop("project must have a pca analysis already completed")
  }
  assert_in(num_components, c(2,3))
  if (!is.null(grouping_var)) {
    assert_single_string(grouping_var)
    assert_in(grouping_var, names(project$data_processed$data_sample))
    assert_string(project$data_processed$data_sample[,grouping_var])
  }
  
  # check number of components
  pca <- project$analyses$pca
  n <- ncol(pca$x)
  
  # create colour vector
  group <- rep("",n)
  if (!is.null(grouping_var)) {
    group <- project$data_processed$data_sample[,grouping_var]
  }
  col_vec <- divergent_colours(length(unique(group)))
  
  # scatterplot of first 2 principal components
  if (num_components == 2) {
    ret <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
                   color = group, type = "scatter", colors = col_vec,
                   mode = "markers", marker = list(size = 5))
  }
  
  # scatterplot of first 3 principal components
  if (num_components == 3) {
    ret <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3,
                   color = group, type = "scatter3d", colors = col_vec,
                   mode = "markers", marker = list(size = 3))
  }
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Plot variance explained by PCA components
#'
#' @description Plot the variance explained by each component. The number of 
#'   components shown is controlled by \code{num_components}, with up to the
#'   first 10 components shown by default.
#'
#' @param project a \code{mipmapper_project} containing pca analysis output
#' @param num_components number of components to be shown
#'
#' @export

plot_pca_variance <- function(project, num_components = 10) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$analyses$pca) == 0) {
    stop("project must have a pca analysis already completed")
  }
  assert_single_pos_int(num_components, zero_allowed = FALSE)
  
  # check number of components
  pca <- project$analyses$pca
  nc <- min(length(colnames(pca$x)), num_components)
  
  # x has to be factored in order to preserve PC order
  x <- factor(colnames(pca$x)[seq_len(nc)],
              levels = colnames(pca$x)[seq_len(nc)])
  
  ret <- plot_ly(x = x, y = pca$var[seq_len(nc)], type = "bar") %>%
    plotly::layout(title = "Percentage of total variance explaned by PCA",
                   xaxis = list(title = "Principal Component"),
                   yaxis = list(title = "Percentage of variance explained"))
  
  # return
  return(ret)
}
