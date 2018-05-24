#' Pipe operator
#'
#' See \code{\link[magrittr:pipe]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom dplyr %>%
#' @export
#' @usage lhs \%>\% rhs
NULL

#' @title Dummy data
#'
#' @description Create a dummy data set
#'
#' @details Creat than a given number of unique reads (coverage). 
#'
#' @return Invisibly returns a MIP dataset with the following variables:
#'   \itemize{
#'       \item{"Year"}{ Numerics of sample year collect }
#'       \item{"Country"}{ String of country of sample collection }
#'       \item{"Study"}{ String of study name }
#'       \item{"DHS_Cluster"}{ Numeric of DHS cluster}
#'       \item{"Lat"}{ Numeric for latitude }
#'       \item{"Long"}{ Numeric for longitude}
#'       \item{"Sample_ID"}{ String for unique sample ID }
#'       \item{"Chrom"}{ Numeric for chromosome number }
#'       \item{"Pos"}{ Numeric for chromosome position}
#'       \item{"Ref"}{ Reference allele nucleotide}
#'       \item{"Alt"}{ Alternative allele nucleotide}
#'       \item{"Coverage"}{ Numeric for read coverage}
#'       \item{"Barcode_Count"}{ Numeric for count of barcodes matched to Ref}
#'       }
#' @export
#' @examples
#' dat <- dummy_data()

dummy_data <- function() {

    df <- data.frame(
    "Year" = rep(2011, 11),
    "Country" = c(rep("mipmappia", 6), rep(NA, 5)),
    "Study" = c(rep("main", 6), rep("PC", 5)),
    "DHS_Cluster" = c(rep(1, 6), rep(NA, 5)),
    "Lat" = c(rep(26.432513, 6), rep(NA, 5)),
    "Long" = c(rep(-46.536719, 6), rep(NA, 5)),
    "Sample_ID" = c(rep(paste0("samp", 1:2), 3), paste0("PC", 1:5)),
    "Chrom" = c(1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1),
    "Pos" = c(100, 200, 100, 100, 100, 300, 100, 200, 100, 100, 200),
    "Ref" = c("A", "C", "TC", "A", "A", "C", "T", "G", "A", "C", "T"),
    "Alt" = c("G", "A", "T", "T", "T", "A", "TCCAGAGACT", "A", "C", "T", "A"),
    "Coverage" = c(4, 80, 600, 400, 4, 200, 0, 457, 341, NA, -4),
    "Barcode_Count" = c(0, 0, 600, 400, 0, 120, 0, NA, 341, 0, 100),
    stringsAsFactors = FALSE
  )

  invisible(df)

}
