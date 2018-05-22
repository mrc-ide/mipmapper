
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @importFrom data.table fread
#' @import dplyr
#' @import graphics
NULL

#------------------------------------------------
#' @title Filter data
#'
#' @description Filter out miscellaneous unwanted values from data
#'
#' @details TODO
#'
#' @param dat TODO
#' @param SNP_only TODO
#' @param group_Alt TODO
#' @param drop_irregular TODO
#'
#' @export
#' @examples
#' # TODO

filter_misc <- function(dat, SNP_only = TRUE, group_Alt = TRUE, drop_irregular = TRUE) {

  # TODO - check format of dat

  # subset to single nucleotide polymorphisms
  if (SNP_only) {
    dat <- subset(dat, nchar(dat$Ref)==1 & nchar(dat$Alt)==1)
  }

  # group all ALT alleles together
  if (group_Alt) {
    dat2 <- dat %>%
      group_by(.dots = setdiff(names(dat), c("Alt", "Barcode_Count"))) %>%
      summarize(Alt = paste(Alt, collapse = ","), Barcode_Count = sum(Barcode_Count))
    dat2 <- as.data.frame(dat2)
    dat <- dat2[,names(dat)]
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

#------------------------------------------------
#' @title Filter data
#'
#' @description Filter data based on coverage
#'
#' @details TODO
#'
#' @param dat TODO
#' @param min_coverage TODO
#'
#' @export
#' @examples
#' # TODO

filter_coverage <- function(dat, min_coverage = 2) {

  # checks on inputs
  assert_pos_int(min_coverage)

  # drop low coverage barcodes
  dat <- subset(dat, Coverage >= min_coverage)

  # return invisibly
  invisible(dat)
}

