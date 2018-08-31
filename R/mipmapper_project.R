
#------------------------------------------------
#' @title Create a new mipmapper project
#'
#' @description Create a new mipmapper project.
#'
#' @export

mipmapper_project <- function() {
  
  # initialise project with default values
  ret <- list(data_raw = list(),
              data_processed = list(),
              analyses = list()
              )
  
  # create class and return
  class(ret) <- "mipmapper_project"
  return(ret)
}

#------------------------------------------------
#' @title Custom print function for class mipmapper_project
#'   
#' @description Custom print function for class \code{mipmapper_project},
#'   printing a summary of the key elements (also equivalent to
#'   \code{summary(x)}). To do an ordinary \code{print()} of all elements of the
#'   project, use the \code{print_full()} function.
#'   
#' @param x object of class \code{mipmapper_project}
#' @param ... other arguments (ignored)
#'   
#' @export

print.mipmapper_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class mipmapper_project
#'
#' @description Calling \code{print()} on an object of class
#'   \code{mipmapper_project} results in custom output. This function therefore
#'   stands in for the base \code{print()} function, and is equivalent to
#'   running \code{print(unclass(x))}.
#'
#' @param x object of class \code{mipmapper_project}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "mipmapper_project")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class mipmapper_project
#'   
#' @description Custom summary function for class \code{mipmapper_project}.
#'   
#' @param object object of class \code{mipmapper_project}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.mipmapper_project <- function(object, ...) {
  
  # summarise raw data
  if (length(object$data_raw) > 0) {
    n_samples_raw <- nrow(object$data_raw$data_coverage)
    n_loci_raw <- ncol(object$data_raw$data_coverage)
    
    cat("RAW DATA\n")
    cat(sprintf("  samples = %s\n", n_samples_raw))
    cat(sprintf("  loci = %s\n", n_loci_raw))
  }
  cat("\n")
  
  # summarise processed data
  if (length(object$data_processed) > 0) {
    n_samples_processed <- nrow(object$data_processed$data_coverage)
    n_loci_processed <- ncol(object$data_processed$data_coverage)
    
    cat("PROCESSED DATA\n")
    cat(sprintf("  samples = %s (%s%%)\n", n_samples_processed, round(n_samples_processed/n_samples_raw*100)))
    cat(sprintf("  loci = %s (%s%%)\n", n_loci_processed, round(n_loci_processed/n_loci_raw*100)))
  }
  
}

#------------------------------------------------
#' @title Determine if object is of class mipmapper_project
#'
#' @description Determine if object is of class \code{mipmapper_project}.
#'
#' @param x object of class \code{mipmapper_project}
#'
#' @export

is.mipmapper_project <- function(x) {
  inherits(x, "mipmapper_project")
}

