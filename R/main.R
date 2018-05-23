
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

#------------------------------------------------
#' @title Melt mip data
#'
#' @description Melt the mip data into a data frame with one column per locus. 
#'   The read depth for each locus, denoted by `chr<x>_<pos>`, is given for 
#'   each sample. An NA is given at loci where there is no read depth.
#'
#' @param data MIP data. The data must have the following variables:
#'   \itemize{
#'       \item Chrom : The chromosome number of the MIP read
#'       \item Pos : The chromosome position of the MIP read
#'       \item Sample_ID : The name of the sample
#'       \item Coverage : The total read coverage as numerics
#'       \item Barcode_Count : The total barocdes recovered
#'       }
#'
#' @export
#' @examples
#' 
#' data <- data.frame(
#' "Sample_ID" = c(rep("a", 3), rep("b", 2)),
#' "Chrom" = c(1, 1, 2, 1, 1),
#' "Pos" = c(100, 200, 50, 100, 200),
#' "Coverage" = c(47, 95, 100, 52, 100),
#' "Barcode_Count" = c(47, 0, 40, 52, 70)
#' )
#' 
#' melt_mip_data(data)
#' 

melt_mip_data <- function(data) {
  
  # make unique identifier for the genome
  data$ID <- paste0(data$Chrom, "_", data$Pos)
  
  # these are all the genome potential variables
  depth_related <- c("Chrom", "Pos", "Ref", "Alt", 
                     "Coverage", "Barcode_Count", "ID")
  meta <- which(!names(data) %in% depth_related)
  
  
  df <- data.frame(matrix(ncol = length(unique(data$ID)) + length(meta), nrow = 0))
  colnames(df) <- c(names(data)[meta],unique(data$ID))
  
  s <- unique(data$Sample_ID)
  slist <- list(); length(slist) <- length(s)+1
  slist[[1]] <- df
  
  for(i in 1:length(s)){
    d <- data[data$Sample_ID==s[i],]  
    df <- data.frame(matrix(ncol = length(d$Barcode_Count),
                            nrow=1, 
                            data = d$Barcode_Count / d$Coverage))
    names(df) <- d$ID
    df <- cbind(d[1, meta, drop = FALSE], df)
    slist[[i+1]] <- df
  }
  
  res <- data.table::rbindlist(slist, fill=TRUE) %>% as.data.frame()
  return(res)
  
}

