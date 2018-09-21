
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
#' @details Load a file from within the mipmapper package. Must be one of .csv, .txt, or .rds.
#'
#' @param name the file name
#'
#' @export

mipmapper_file <- function(name) {
  
  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  assert_in(ext, c("txt", "csv", "rds"), message = "file extension not valid")
  
  # get full file path
  name_full <- system.file("extdata/", name, package='mipmapper', mustWork = TRUE)
  
  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <- fast_read(name_full)
  }
  
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

#------------------------------------------------
#' @title Pairwise Great-Circle Distance
#'
#' @description Analogue of \code{dist()} function, but returning great-circle
#'   distances.
#'
#' @param points a matrix with two columns specifying latitude and longitude
#'
#' @export
#' @examples
#' some_points <- cbind(1:10, 1:10)
#' dist_gc(some_points)

dist_gc <- function(points) {
  
  # check inputs
  assert_matrix(points)
  assert_eq(ncol(points), 2)
  
  # calculate distance matrix
  d <- apply(points, 1, function(x) {
    latlon_to_bearing(x[1], x[2], points[,1], points[,2])$gc_dist
    })
  diag(d) <- 0
  
  return(d)
}

#------------------------------------------------
#' @title latlon_to_bearing
#'
#' @description Calculate distance between latitude/longitude points
#'
#' @param origin_lat The origin latitude
#' @param origin_lon The origin longitude
#' @param dest_lat The destination latitude
#' @param dest_lon The destination longitude
#'
#' @export

latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
  
  # convert input arguments to radians
  origin_lat <- origin_lat*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  origin_lon <- origin_lon*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  
  delta_lon <- dest_lon-origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  gc_angle[is.nan(gc_angle)] <- 0
  
  # convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earthRad <- 6371
  gc_dist <- earthRad*gc_angle
  
  return(list(bearing = bearing, gc_dist = gc_dist))
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}
