#  filter_misc -----------------------------------------------------------------
#' @title Filter data
#'
#' @description Filter out miscellaneous unwanted values from data
#'
#' @details `filter_misc` takes a raw dataset and applies a series of filters
#'   that enable the removal of non-standarad ALT allele calls and long ALT
#'   calls. It also enables the filtering of irregular data, such as data
#'   with non integer counts of barcodes or more barcode counts than coverage.
#'
#' @param dat MIP data. The data must have the following variables:
#'   \itemize{
#'       \item Ref : Character for reference allele nucleotide
#'       \item Alt : Character for alternative allele nucelotide
#'       \item Coverage : The total read coverage as numerics
#'       \item Barcode_Count : The total barocdes recovered
#'       }
#' @param SNP_only Boolean detailing if subsetting should remove all cases
#'   where the ALT allele is longer than 1 nucleotide. Default = `TRUE`
#' @param group_Alt Boolean detailing whether all ALT alleles should be 
#'   grouped together. Default = `TRUE`
#' @param drop_irregular Boolean detailing whether irregular data should
#'   be removed. This includes non integer barcode counts and any instances 
#'   where there is less coverage than barcode counts. Default = `TRUE`
#'   
#' @importFrom dplyr group_by summarize
#' @export
#' @examples
#' dat <- dummy_data()
#' dat2 <- filter_misc(dat)
#' dim(dat2)[1] < dim(dat)[1]


filter_misc <- function(dat,
                        SNP_only = TRUE,
                        group_Alt = TRUE,
                        drop_irregular = TRUE) {

  # TODO - check format of dat
  # Worth a discussion of the above - best to be some check_dat() where we
  # find and standardise names

  # subset to single nucleotide polymorphisms
  if (SNP_only) {
    dat <- subset(dat, nchar(dat$Ref) == 1 & nchar(dat$Alt) == 1)
  }

  # group all ALT alleles together
  if (group_Alt) {
    dat2 <- dat %>%
      group_by(.dots = setdiff(names(dat), c("Alt", "Barcode_Count"))) %>%
      summarize(Alt = paste(Alt, collapse = ","),
                Barcode_Count = sum(Barcode_Count))
    dat2 <- as.data.frame(dat2)
    dat <- dat2[, names(dat)]
  }

  # filter barcode counts. Drop non-int values and those where Barcode_Count
  # exceeds Coverage
  if (drop_irregular) {
    dat <- subset(dat, is.int_vector(Barcode_Count) & is.int_vector(Coverage))
    dat <- subset(dat, Barcode_Count <= Coverage)
  }

  # return invisibly
  invisible(dat)
}


#  filter_coverage -------------------------------------------------------------
#' @title Filter data
#'
#' @description Filter data based on coverage
#'
#' @details The MIP dataset is filtered to remove all loci that have less than 
#'   the minimum required coverage. 
#'
#' @param dat MIP data. The data must have the following variables:
#'   \itemize{
#'       \item Coverage : The total read coverage as numerics
#'       }
#' @param min_coverage Numeric for the minimum required coverage. Default = 2
#'
#' @export
#' @examples
#' dat <- dummy_data()
#' dat <- filter_misc(dat)
#' dat2 <- filter_coverage(dat, min_coverage = 5)
#' dim(dat2)[1] < dim(dat)[1]


filter_coverage <- function(dat, min_coverage = 2) {

  # checks on inputs
  assert_pos_int(min_coverage)

  # drop low coverage barcodes
  dat <- subset(dat, Coverage >= min_coverage)

  # return invisibly
  invisible(dat)
}


#  melt_mip_data ---------------------------------------------------------------
#' @title Melt mip data
#'
#' @description Melt the mip data into a data frame with one column per locus. 
#'   The read depth for each locus, denoted by `chr<x>_<pos>`, is given for 
#'   each sample. An NA is given at loci where there is no read depth.
#'
#' @param dat MIP data. The data must have the following variables:
#'   \itemize{
#'       \item Chrom : The chromosome number of the MIP read
#'       \item Pos : The chromosome position of the MIP read
#'       \item Sample_ID : The name of the sample
#'       \item Coverage : The total read coverage as numerics
#'       \item Barcode_Count : The total barocdes recovered
#'       }
#'
#' @return Invisibly returns the melted data frame
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' 
#' dat <- data.frame(
#' "Sample_ID" = c(rep("a", 3), rep("b", 2)),
#' "Chrom" = c(1, 1, 2, 1, 1),
#' "Pos" = c(100, 200, 50, 100, 200),
#' "Coverage" = c(47, 95, 100, 52, 100),
#' "Barcode_Count" = c(47, 0, 40, 52, 70)
#' )
#' 
#' melt_mip_data(dat = dat)
#' 

melt_mip_data <- function(dat) {

  # make unique identifier for the genome
  dat$ID <- paste0(dat$Chrom, "_", dat$Pos)

  # these are all the genome potential variables
  depth_related <- c("Chrom", "Pos", "Ref", "Alt",
                     "Coverage", "Barcode_Count", "ID")
  meta <- which(!names(dat) %in% depth_related)

  # the total width of our final dataset
  cols <- length(unique(dat$ID)) + length(meta)

  df <- data.frame(matrix(ncol = cols, nrow = 0))
  colnames(df) <- c(names(dat)[meta], unique(dat$ID))

  s <- unique(dat$Sample_ID)
  slist <- list()
  length(slist) <- length(s) + 1
  slist[[1]] <- df

  for (i in seq_len(length(s))){
    d <- dat[dat$Sample_ID == s[i], ]
    df <- data.frame(matrix(ncol = length(d$Barcode_Count),
                            nrow = 1,
                            data = d$Barcode_Count / d$Coverage))
    names(df) <- d$ID
    df <- cbind(d[1, meta, drop = FALSE], df)
    slist[[i + 1]] <- df
  }

  res <- rbindlist(slist, fill = TRUE) %>% as.data.frame()
  invisible(res)

}


#  impute_mip_data -------------------------------------------------------------
#' @title Impute missing mip data
#'
#' @description Using the output of \code{melt_mip_data}, missing values are 
#'  imputed by applying a summary function to the non NA values for a locus. 
#'  The default summary function takes the mean of the non NA samples.
#'  
#'
#' @param dat output of \code{melt_mip_data}
#' @param FUN function to impute missing values. Default = `mean`
#' @param ... other arguments to pass to FUN. 
#'
#' @return Invisibly returns the mip data frame with missing values imputed
#' @export
#' @examples
#' 
#' dat <- data.frame(
#' "Sample_ID" = c(rep("a", 3), rep("b", 2)),
#' "Chrom" = c(1, 1, 2, 1, 1),
#' "Pos" = c(100, 200, 50, 100, 200),
#' "Coverage" = c(47, 95, 100, 52, 100),
#' "Barcode_Count" = c(47, 0, 40, 52, 70)
#' )
#' 
#' dat <- melt_mip_data(dat = dat)
#' dat <- impute_mip_data(dat = dat)
#' 

impute_mip_data <- function(dat, FUN = mean, ...) {

  # get which columns refer to loci and which refer to variable names
  dat_names <- names(dat)[which(!grepl("^chr.*", names(dat)))]
  locus_names <- setdiff(names(dat), dat_names)

  # TODO: Discuss at some point standard variable handling, where PCs will be
  # subset out positive controls


  # impute NAs
  dat_mat <- as.matrix(dat[, locus_names])
  locus_impute <- apply(dat_mat, 2, FUN, na.rm = TRUE, ...)
  locus_impute <- outer(rep(1, nrow(dat_mat)), locus_impute)
  dat_mat[is.na(dat_mat)] <- locus_impute[is.na(dat_mat)]
  dat <- cbind(dat[, dat_names], dat_mat)

  # return invisbly
  invisible(dat)

}
