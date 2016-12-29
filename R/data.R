#' United States county centers.
#'
#' A dataset containing the latitudes and longitudes for United
#' States spatial and population-weighted county centers from
#' 2000 and 2010 census.
#'
#' @format A data frame with 3147 rows and 9 variables:
#' \describe{
#'   \item{fips}{5 digit county FIPS codes}
#'   \item{clon00}{Spatial center longitude, 2000}
#'   \item{clat00}{Spatial center latitude, 2000}
#'   \item{clon10}{Spatial center longitude, 2010}
#'   \item{clat10}{Spatial center latitude, 2010}
#'   \item{pclon00}{Population-weighted spatial center longitude, 2000}
#'   \item{pclat00}{Population-weighted Spatial center latitude, 2000}
#'   \item{pclon10}{Population-weighted Spatial center longitude, 2010}
#'   \item{pclat10}{Population-weighted Spatial center latitude, 2010}
#' }
#' @source \url{https://github.com/btskinner/spatial}
"county_centers"
