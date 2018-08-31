
#------------------------------------------------
#' @title Impute missing values
#'
#' @description Missing values in the input matrix are imputed by applying a 
#'   summary function over rows or columns. The summary function is only applied to
#'   non-NA values. The default summary function takes the mean.
#'
#' @param dat a matrix of values
#' @param MARGIN which margin to apply summary function over: 1 = rows, 2 = 
#'   columns
#' @param FUN function to impute missing values. Default = `mean`
#' @param ... other arguments to pass to FUN
#'
#' @export

impute <- function(dat, MARGIN = 2, FUN = mean, ...) {
  
  # check inputs
  assert_matrix(dat)
  assert_in(MARGIN, c(1,2))
  
  # impute values
  dat <- apply(dat, MARGIN, function(x) {
                              x_summary <- FUN(x[!is.na(x)])
                              x[is.na(x)] <- x_summary
                              return(x)
                            })
  
  # rotate dat if needed
  if (MARGIN == 1) {
    dat <- t(dat)
  }
  
  return(dat)
}

#------------------------------------------------
#' @title PCA of MIP Data
#'
#' @description Using the output of \code{impute_mip_data}, principal component
#'   analysis will be condcuted and the resultant components returned, along
#'   with the variance explained by each component and the loadings of each
#'   component.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#' @param plot_on whether to plot PCA results
#'
#' @export

pca_mip_data <- function(project, plot_on = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  assert_single_logical(plot_on)
  
  # impute missing values
  WSAF <- impute(project$data_processed$WSAF, MARGIN = 2, FUN = mean)
  
  # compute PCA
  pca <- prcomp(WSAF)
  
  # compute variance explained
  pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100
  
  # compute loadings from rotations
  pca$loadings <- abs(pca$rotation)
  pca$loadings <- sweep(pca$loadings, 2, colSums(pca$loadings), "/")
  
  # save pca analysis results
  project$analyses$pca <- pca
  
  return(project)
}
