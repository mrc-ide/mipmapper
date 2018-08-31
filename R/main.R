
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @importFrom raster raster plot
#' @importFrom viridis plasma
#' @importFrom dplyr group_by summarize
#' @importFrom data.table rbindlist fread
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot abline hist par
#' @importFrom grDevices colorRampPalette grey
#' @importFrom stats prcomp runif
NULL

#------------------------------------------------
#' @title Bind data to project
#'   
#' @description Load data into a \code{mipmapper_project} prior to analysis. See
#'   details for the required data format.
#'   
#' @details Data must be read in in four parts. Although this may seem 
#'   cumbersome, it provides an efficient way of storing sample-locus data 
#'   without duplication of elements (as tends to be found in "long" format).
#'   This reduces the amount of memory required to store the data.
#' \itemize{
#'       \item{\code{data_coverage}}{ a matrix giving the coverage (i.e. the
#'       total number of reads) for each sample-locus combination. Samples are
#'       in rows and loci in columns }
#'       
#'       \item{\code{data_counts}}{ a matrix giving the barcode counts (i.e. the
#'       number of reads matched to the Ref allele) for each sample-locus
#'       combination }
#'       
#'       \item{\code{data_sample}}{ a dataframe giving details of each sample. 
#'       This must have the same number of rows as \code{data_coverage}, and may
#'       have any number of columns as long as it contains the column
#'       "Sample_ID" with a unique ID for each sample }
#'       
#'       \item{\code{data_loci}}{ a dataframe giving details of each locus. This
#'       must have the same number of rows as \code{data_coverage} has columns, 
#'       and may have any number of columns as long as it contains the columns 
#'       "Chrom" for the chromosome, "Pos" for the genomic position on the 
#'       chromosome, "Ref" for the reference allele, and "Alt" for the
#'       alternative allele }
#'       }
#'   
#' @param project a \code{mipmapper_project}, as produced by the function 
#'   \code{mipmapper_project()}
#' @param data_coverage a matrix giving the number of reads for each
#'   sample-locus combination (see details)
#' @param data_counts a matrix giving the number of reads matched to the Ref
#'   allele for each sample-locus combination (see details)
#' @param data_sample a dataframe giving details of each sample (see details)
#' @param data_locus a dataframe giving details of each locus (see details)
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data/outputs
#'
#' @export

bind_data <- function(project,
                      data_coverage,
                      data_counts,
                      data_sample,
                      data_locus,
                      name = NULL,
                      check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  assert_matrix(data_coverage)
  apply(data_coverage, 2, assert_pos)
  assert_matrix(data_counts)
  apply(data_counts, 2, assert_pos)
  assert_dim(data_counts, dim(data_coverage))
  assert_dataframe(data_sample)
  assert_in("Sample_ID", names(data_sample))
  assert_string(data_sample$Sample_ID)
  assert_noduplicates(data_sample$Sample_ID)
  assert_dataframe(data_locus)
  assert_in(c("Chrom", "Pos", "Ref", "Alt"), names(data_locus))
  assert_pos_int(data_locus$Chrom)
  assert_pos_int(data_locus$Pos)
  assert_string(data_locus$Ref)
  assert_string(data_locus$Alt)
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (check_delete_output && length(project$data_raw) > 0) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("Re-loading data will erase all existing data and analyses. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- mipmapper_project()
  }
  
  # add data to project
  project$data_raw <- list(data_coverage = data_coverage,
                           data_counts = data_counts,
                           data_sample = data_sample,
                           data_locus = data_locus)
  
  project$data_processed <- project$data_raw
  project$data_processed$WSAF <- data_counts / data_coverage
  
  message("data loaded")
  
  invisible(project)
}

#------------------------------------------------
#' @title Filter miscellaneous values
#'
#' @description Apply miscellaneous filters to data.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#' @param filter_SNP filter to single nucleotide polymorphism (SNP) only. Remove
#'   all loci where Ref or Alt allele is longer than one character
#' @param integer_counts re-code as missing all genotype calls where coverage or
#'   barcode counts are non-integer
#' @param coverage_greq_counts re-code as missing all genotype calls where 
#'   coverage is not greater than or equal to barcode counts
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data/outputs
#'
#' @export

filter_misc <- function(project,
                        filter_SNP = TRUE,
                        integer_counts = TRUE,
                        coverage_greq_counts = TRUE,
                        check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$data_raw) == 0) {
    stop("do data loaded")
  }
  assert_single_logical(filter_SNP)
  assert_single_logical(integer_counts)
  assert_single_logical(coverage_greq_counts)
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (check_delete_output & length(project$analyses) > 0) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("Re-filtering data will erase all existing analyses. Continue? (Y/N): ")) {
      return(project)
    }
    
    # erase or overwrite elements
    project$analyses <- list()
  }
  
  # filter to SNP only
  if (filter_SNP) {
    
    # find non-SNPs
    Ref_single_character <- (nchar(project$data_processed$data_locus$Ref) == 1)
    Alt_single_character <- (nchar(project$data_processed$data_locus$Ref) == 1)
    both_single_character <- (Ref_single_character & Alt_single_character)
    
    # report to console
    message(sprintf("%s loci found to be non-SNPs (%s%% of loci)",
                    sum(!both_single_character, na.rm = TRUE),
                    round(mean(!both_single_character, na.rm = TRUE)*100, 2)))
    
    # filter
    w <- both_single_character
    project$data_processed$data_coverage <- project$data_processed$data_coverage[,w]
    project$data_processed$data_counts <- project$data_processed$data_counts[,w]
    project$data_processed$data_locus <- project$data_processed$data_locus[w,]
  }
  
  # filter to integer counts only
  if (integer_counts) {
    
    # find non-integer values
    is_integer_coverage <- apply(project$data_processed$data_coverage, 2, is_integer)
    is_integer_counts <- apply(project$data_processed$data_counts, 2, is_integer)
    is_integer_both <- (is_integer_coverage & is_integer_counts)
    
    # report to console
    message("")
    message(sprintf("%s coverage calls found to be non-integer (%s%% of calls)",
                    sum(!is_integer_coverage, na.rm = TRUE),
                    round(mean(!is_integer_coverage, na.rm = TRUE)*100, 2)))
    message(sprintf("%s barcode count calls found to be non-integer (%s%% of calls)",
                    sum(!is_integer_counts, na.rm = TRUE),
                    round(mean(!is_integer_counts, na.rm = TRUE)*100, 2)))
    message(sprintf("total %s values re-coded as missing (%s%% of calls)",
                    sum(!is_integer_both, na.rm = TRUE),
                    round(mean(!is_integer_both, na.rm = TRUE)*100, 2)))
    
    # re-code as missing
    w <- which(!is_integer_both)
    project$data_processed$data_coverage[w] <- NA
    project$data_processed$data_counts[w] <- NA
  }
  
  # filter to remove calls where barcode counts exceed coverage
  if (coverage_greq_counts) {
    
    # find counts > coverage
    counts_gr_coverage <- (project$data_processed$data_counts > project$data_processed$data_coverage)
    
    # report to console
    message("")
    message(sprintf("%s barcode counts found to exceed coverage (%s%% of calls)",
                    sum(counts_gr_coverage, na.rm = TRUE),
                    round(mean(counts_gr_coverage, na.rm = TRUE)*100, 2)))
    
    # re-code as missing
    w <- which(counts_gr_coverage)
    project$data_processed$data_coverage[w] <- NA
    project$data_processed$data_counts[w] <- NA
  }
  
  # recalculate WSAF
  project$data_processed$data_WSAF <- project$data_processed$data_counts / project$data_processed$data_coverage
  
  invisible(project)
}

#------------------------------------------------
#' @title Filter data based on coverage
#'
#' @description Filter data based on coverage. First, a minimum coverage level 
#'   is chosen and all genotype calls below this level are re-coded as missing 
#'   data. Second, loci are filtered based on the allowed proportion of missing 
#'   data. Third, samples are filtered in an analogous fashion. By default this
#'   is carried out interactively, although filters can be also be specified as
#'   fixed arguments.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#' @param filter_interactive whether to carry out plotting and data filtering
#'   interactively
#' @param min_coverage minimum required coverage of any genotype call. All calls
#'   with coverage less than this level are re-coded as missing data. Only 
#'   applies if \code{filter_interactive} is \code{FALSE}, otherwise chosen 
#'   interactively
#' @param max_proportion_missing_locus maximum proportion of genotype calls over
#'   all samples at a given locus that are allowed to be missing (see 
#'   \code{min_coverage} argument above). Only applies if
#'   \code{filter_interactive} is \code{TRUE}, otherwise chosen interactively
#' @param max_proportion_missing_sample maximum proportion of genotype calls
#'   over all loci for a given sample that are allowed to be missing (see 
#'   \code{min_coverage} argument above). Only applies if 
#'   \code{filter_interactive} is \code{TRUE}, otherwise chosen interactively
#'
#' @export

filter_coverage <- function(project,
                            filter_interactive = TRUE,
                            min_coverage = 2,
                            max_proportion_missing_locus = 0.5,
                            max_proportion_missing_sample = 0.5) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  assert_single_logical(filter_interactive)
  assert_single_pos_int(min_coverage, zero_allowed = TRUE)
  assert_single_pos(max_proportion_missing_locus, zero_allowed = TRUE)
  assert_bounded(max_proportion_missing_locus, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_pos(max_proportion_missing_sample, zero_allowed = TRUE)
  assert_bounded(max_proportion_missing_sample, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  
  # make copy of data
  data_coverage <- project$data_processed$data_coverage
  data_counts <- project$data_processed$data_counts
  data_sample <- project$data_processed$data_sample
  data_locus <- project$data_processed$data_locus
  
  # choose min_coverage interactively
  if (filter_interactive) {
    
    # setup two plotting rows
    par_save <- par(mfrow = c(2,1))
    on.exit(par(par_save))
    
    # initialise artificial user choice
    min_coverage <- min_coverage + 1
    user_choice <- "a"
    
    # wait for required user input
    while (!user_choice %in% c("s", "x")) {
      
      # re-filter and plot if "a" or "d"
      if (user_choice %in% c("a", "d")) {
        
        # change min_coverage value
        min_coverage <- ifelse(user_choice == "a", min_coverage - 1, min_coverage + 1)
        min_coverage <- max(min_coverage, 0)
        message(sprintf("min_coverage = %s", min_coverage))
        
        # re-filter
        coverage_greq_limit <- (data_coverage >= min_coverage)
        
        # get proportion missing loci and samples
        prop_missing_loci <- 1 - colMeans(coverage_greq_limit, na.rm = TRUE)
        prop_missing_samples <- 1 - rowMeans(coverage_greq_limit, na.rm = TRUE)
        
        # produce plot
        hist(prop_missing_loci, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing samples\nper locus")
        hist(prop_missing_samples, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing loci\nper sample")
      }
      
      # prompt user
      user_choice <- readline("\nUse 'a' and 'd' to decrease/increase the min_coverage level, \n's' to lock in your choce, or 'x' to abort: ")
    }
    
    # return without modification on abort
    if (user_choice == "x") {
      message("aborting without saving changes")
      return(project)
    }
    
    message(sprintf("applying cut-off at min_coverage = %s\n", min_coverage))
  }
  
  # drop low coverage calls
  coverage_greq_limit <- (data_coverage >= min_coverage)
  w <-  which(!coverage_greq_limit)
  data_coverage[w] <- NA
  data_counts[w] <- NA
  prop_missing_loci <- colMeans(is.na(data_coverage))
  
  # choose max_proportion_missing_locus interactively
  if (filter_interactive) {
    
    # initialise artificial user choice
    user_choice <- max_proportion_missing_locus
    
    # wait for required user input
    while (!user_choice %in% c("s", "x")) {
      
      # if numeric value
      if (!is.na(suppressWarnings(as.numeric(user_choice)))) {
        
        # change max_proportion_missing_locus value
        max_proportion_missing_locus <- as.numeric(user_choice)
        message(sprintf("max_proportion_missing_locus = %s", max_proportion_missing_locus))
        
        # re-filter
        w <- which(prop_missing_loci < max_proportion_missing_locus)
        prop_missing_samples <- rowMeans(is.na(data_coverage[, w, drop = FALSE]))
        message(sprintf("remaining loci: %s (%s%%)", length(w),
                        round(length(w)/length(prop_missing_loci)*100, 2)))
        
        # produce plot
        hist(prop_missing_loci, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing samples\nper locus")
        abline(v = max_proportion_missing_locus, col = 2)
        hist(prop_missing_samples, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing loci\nper sample")
        
      }
      
      # prompt user
      user_choice <- readline("\nEnter new max_proportion_missing_locus value,\nor type 's' to lock in your choce, or 'x' to abort: ")
    }
    
    # return without modification on abort
    if (user_choice == "x") {
      message("aborting without saving changes")
      return(project)
    }
    
    message(sprintf("applying cut-off at max_proportion_missing_locus = %s\n", max_proportion_missing_locus))
  }
  
  # drop loci with too much missing data
  w <- which(prop_missing_loci < max_proportion_missing_locus)
  data_coverage <- data_coverage[,w]
  data_counts <- data_counts[,w]
  data_locus <- data_locus[w,,drop = FALSE]
  prop_missing_samples <- rowMeans(is.na(data_coverage))
  
  # choose max_proportion_missing_sample interactively
  if (filter_interactive) {
    
    # initialise artificial user choice
    user_choice <- max_proportion_missing_sample
    
    # wait for required user input
    while (!user_choice %in% c("s", "x")) {
      
      # if numeric value
      if (!is.na(suppressWarnings(as.numeric(user_choice)))) {
        
        # change max_proportion_missing_sample value
        max_proportion_missing_sample <- as.numeric(user_choice)
        message(sprintf("max_proportion_missing_sample = %s", max_proportion_missing_sample))
        
        # re-filter
        w <- which(prop_missing_samples < max_proportion_missing_sample)
        message(sprintf("remaining samples: %s (%s%%)", length(w),
                        round(length(w)/length(prop_missing_samples)*100, 2)))
        
        # produce plot
        hist(prop_missing_loci, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing samples\nper locus")
        abline(v = max_proportion_missing_locus, col = 2)
        hist(prop_missing_samples, breaks = seq(0,1,0.01), col = grey(0.6),
             xlab = "proportion missing", main = "proportion missing loci\nper sample")
        abline(v = max_proportion_missing_sample, col = 2)
        
      }
      
      # prompt user
      user_choice <- readline("\nEnter new max_proportion_missing_sample value,\nor type 's' to lock in your choce, or 'x' to abort: ")
    }
    
    # return without modification on abort
    if (user_choice == "x") {
      message("aborting without saving changes")
      return(project)
    }
    
    message(sprintf("applying cut-off at max_proportion_missing_sample = %s\n", max_proportion_missing_sample))
  }
  
  # drop samples with too much missing data
  w <- which(prop_missing_samples < max_proportion_missing_sample)
  data_coverage <- data_coverage[w,]
  data_counts <- data_counts[w,]
  data_sample <- data_sample[w,,drop = FALSE]
  
  # save back into project
  project$data_processed$data_coverage <- data_coverage
  project$data_processed$data_counts <- data_counts
  project$data_processed$data_sample <- data_sample
  project$data_processed$data_locus <- data_locus
  
  # recalculate WSAF
  project$data_processed$WSAF <- data_counts / data_coverage
  
  # return invisibly
  invisible(project)
}
