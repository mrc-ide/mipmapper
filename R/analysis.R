
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
#' @title Perform PCA analysis
#'
#' @description Principal component analysis will be condcuted and the resultant
#'   components returned, along with the variance explained by each component
#'   and the loadings of each component.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#'
#' @export

perform_pca <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  assert_single_logical(plot_on)
  
  # impute missing values
  WSAF <- impute(project$data_processed$data_WSAF, MARGIN = 2, FUN = mean)
  
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

#------------------------------------------------
#' @title Estimate pairwise inbreeding coefficient F by maximum likelihood
#'
#' @description The probability of seeing the same or different alleles at a 
#'   locus can be written in terms of the global allele frequency p and the 
#'   inbreeding coefficient f. For example, the probability of seeing the same 
#'   REF allele is \code{(1-f)*p^2 + f*p}. This formula can be multiplied over 
#'   all loci to arrive at the overall likelihood of each value of f, which can 
#'   then be chosen by maximum likelihood. This function carries out this
#'   comparison between all pairwise samples, passed in as a matrix.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#' @param x a matrix with samples in rows and loci in columns. The value in each
#'   cell should be 1 or 0, or \code{NA} for missing data
#' @param f_breaks the number of values of f that are explored, distributed
#'   evenly in the [0,1] interval
#' @param report_progress this analysis can take a long time for large datasets,
#'   therefore this option prints running progress to the console
#'
#' @export

estimate_f <- function(project, samples = NULL, f_breaks = 11, report_progress = FALSE) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  samples <- define_default(samples, 1:nrow(project$data_processed$data_sample))
  assert_pos_int(samples, zero_allowed = FALSE)
  assert_noduplicates(samples)
  assert_leq(max(samples), nrow(project$data_processed$data_sample))
  assert_single_pos_int(f_breaks, zero_allowed = FALSE)
  assert_single_logical(report_progress)
  
  # call major allele at each locus
  WSAF_major <- round(project$data_processed$data_WSAF[samples,])
  
  # get global allele frequencies
  p <- colMeans(WSAF_major, na.rm = TRUE)
  
  # convert NA to -1 before passing to C++
  WSAF_major[is.na(WSAF_major)] <- -1
  
  # run efficient C++ function
  args <- list(x = mat_to_rcpp(WSAF_major), f_breaks = f_breaks, p = p, report_progress = report_progress)
  output_raw <- estimate_f_cpp(args)
  
  # process output
  ret <- rcpp_to_mat(output_raw$ret)
  diag(ret) <- 1
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate proportion identical by state (IBS)
#'
#' @description Calculate the proportion of identity by state (IBS) between all
#'   pairwise samples. IBS is calculated as the proportion of sites that are
#'   identical.
#'
#' @param x a matrix with samples in rows and loci in columns. The value in each
#'   cell should be 1, 0, NA for missing data
#' @param report_progress this analysis can take a long time for large datasets,
#'   therefore this option prints running progress to the console
#'
#' @export

calculate_IBS <- function(x, report_progress = FALSE) {
  
  # check inputs
  assert_matrix(x)
  assert_in(as.vector(x), c(0,1,NA,NaN))
  assert_single_logical(report_progress)
  
  # convert NA to -1 before passing to C++
  x[is.na(x)] <- -1
  
  # run efficient C++ function
  args <- list(x = mat_to_rcpp(x), report_progress = report_progress)
  output_raw <- calculate_IBS_cpp(args)
  
  # process output
  ret <- rcpp_to_mat(output_raw$ret)
  diag(ret) <- 1
  
  return(ret)
}
