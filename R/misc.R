
#------------------------------------------------
#' @title Quickly load data as data.frame
#'
#' @description Quickly load data as data.frame
#'
#' @details TODO
#'
#' @param file TODO
#'
#' @export
#' @examples
#' # TODO

fast_read <- function(file) {
  fread(file, data.table = FALSE)
}

#------------------------------------------------
# check if values can be coerced to integer
# (not exported)

is.int_vector <- function(x) {
  x == as.integer(x)
}
