
#------------------------------------------------
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

#------------------------------------------------
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

#------------------------------------------------
#' @title Quickly load data as dataframe
#'
#' @description Quickly load data as a dataframe.
#'
#' @param file path to where your data is saved
#'
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

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE
#' @noRd
user_yes_no <- function(x="continue? (Y/N): ") {
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

# -----------------------------------
# test if x is integer-valued (as opposed to being of class "integer")
#' @noRd
is_integer <- function(x) {
  as.integer(x) == x
}
