#  pca_mip_data -------------------------------------------------------------
#' @title PCA of MIP Data
#'
#' @description Using the output of \code{impute_mip_data}, principal 
#'   component analysis will be condcuted and the resultant components 
#'   returned, with the variance in the data explained by each component 
#'   and the loadings of each component also returned.
#'  
#' @param dat output of \code{impute_mip_data}
#'
#' @return Invisibly returns a list of class `prcomp` with the following 
#'   components
#'   \itemize{
#'       \item{"sdev"}{ the standard deviations of the principal components
#'       (i.e., the square roots of the eigenvalues of the
#'       covariance/correlation matrix, though the calculation is actually
#'       done with the singular values of the data matrix).}
#'       \item{"rotation"}{ the matrix of variable loadings (i.e., a matrix
#'       whose columns contain the eigenvectors). The function princomp
#'       returns this in the element loadings.}
#'       \item{"center, scale"}{ the centering and scaling used }
#'       \item{"x"}{ the value of the rotated data (the centred data multiplied
#'       by the rotation matrix) is returned. Hence, cov(x) is the diagonal
#'       matrix diag(sdev^2).}
#'       \item{"var"}{ the variance in the data explained by each component }
#'       \item{"loadings"}{ the loadings of each component }
#'       \item{"dat"}{ the raw data used to conduct pca }
#'       }
#' @importFrom stats prcomp
#' @export
#' @examples
#' 
#' dat <- dummy_data()
#' dat <- filter_misc(dat = dat)
#' dat <- filter_coverage(dat = dat, min_coverage = 2)
#' dat <- melt_mip_data(dat = dat)
#' dat <- impute_mip_data(dat = dat)
#' pca <- pca_mip_data(dat = dat)
#' 

pca_mip_data <- function(dat) {

  # get which columns refer to loci and which refer to variable names
  dat_names <- names(dat)[which(!grepl("^chr.*", names(dat)))]
  locus_names <- setdiff(names(dat), dat_names)

  # compute PCA
  pca <- prcomp(dat[, locus_names])

  # compute variance explained
  pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100

  # compute loadings from rotations
  pca$loadings <- abs(pca$rotation)
  pca$loadings <- sweep(pca$loadings, 2, colSums(pca$loadings), "/")

  # and attach the data
  pca$dat <- dat

  # return this invisibly
  invisible(pca)

}
