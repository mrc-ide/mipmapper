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

#  dummy_data ------------------------------------------------------------------
#' @title Dummy data
#'
#' @description Create a dummy data set
#' 
#' @details Fixed dummy data set. Mainly used for tests. For a more nuanced 
#'   dummy data set creation see \code{generate_dummy_data}.
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

#  generate_dummy_data ---------------------------------------------------------
#' @title Generate Dummy Data
#'
#' @description Create a dummy data set, with control over the creation.
#'
#' @details Creates a dummy data set, with more control over how the data
#'   set variables are created
#' @param n Number of rows. Default = 1000
#' @param countries Number of different countries. Default = 4
#' @param studies Number of different studies. Default = 3
#' @param samples Number of different samples. Default = 50
#' @param pos Number of different chromosome positions. Default = 20
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
#' @importFrom stats runif
#' @examples
#' dat <- rdhs:::generate_dummy_data()

generate_dummy_data <- function(n = 1000, countries = 4,
                                studies = 3, samples = 50,
                                pos = 20) {


  chroms <- rep(rep(1:10, samples), length.out = n)
  positions <- sort(rep(sample(1000, pos, TRUE), length.out = n))
  refs <- sort(rep(sample(c("A", "C", "G", "T"), pos, TRUE), length.out = n))
  alts <- sort(rep(sample(c("A", "C", "G", "T"), pos, TRUE), length.out = n),
               decreasing = TRUE)
  sample_names <- sort(sample(paste0("samp", 1:samples), n, TRUE))
  country <- rep(sample(paste0("c", seq_len(countries)), samples, TRUE),
                 table(sample_names))
  study <- rep(sample(paste0("st", seq_len(studies)), samples, TRUE),
               table(sample_names))
  dhs <- rep(sample(n / 10, samples, TRUE), table(sample_names))
  lat <- rep(runif(samples, - 50, 50), table(sample_names))
  lon <- rep(runif(samples, - 50, 50), table(sample_names))
  year <- rep(sample(2010:2016, samples, TRUE),
              table(sample_names))

  coverage <- sample(1000, 1000, TRUE)
  barcode_counts <- apply(cbind(coverage, country), 1, function(x) {
    low <- as.numeric(gsub("c([0-9]+)$", "\\1", x[2])) / countries
    sample(seq_len(as.numeric(x[1])), 1,
      prob = seq(low - .25 + 0.00001, low, length.out = as.numeric(x[1]))
    )
  }
  )


  df <- data.frame(
    "Year" = year,
    "Country" = country,
    "Study" = study,
    "DHS_Cluster" = dhs,
    "Lat" = lat,
    "Long" = lon,
    "Sample_ID" = sample_names,
    "Chrom" = chroms,
    "Pos" = positions,
    "Ref" = refs,
    "Alt" = alts,
    "Coverage" = coverage,
    "Barcode_Count" = barcode_counts,
    stringsAsFactors = FALSE
  )

  invisible(df)

}
