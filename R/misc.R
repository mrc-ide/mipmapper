
#------------------------------------------------
#' @title Quickly load data as data.frame
#'
#' @description Quickly load data as data.frame
#'
#' @details TODO
#'
#' @param file String to the path where your mip dataset is saved
#'
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' \dontrun{
#' file <- "mip_files/dataset.csv"
#' dat <- fast_read(file)
#' }

fast_read <- function(file) {
  fread(file, data.table = FALSE)
}

#------------------------------------------------
# check if values can be coerced to integer
# (not exported)
#' @noRd
is.int_vector <- function(x) {
  x == as.integer(x)
}

#  load system file ---------------------------------------------------------------
#' @title Load system file
#'
#' @description Load system file
#'
#' @details Load a file from within the mipmapper package.
#'
#' @param name the file name
#'
#' @export
#' @examples
#' # TODO

mipmapper_file <- function(name) {
  name_full <- system.file("extdata/", name, package='mipmapper', mustWork = TRUE)
  ret <- fast_read(name_full)

  return(ret)
}

