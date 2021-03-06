# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute distance between each coordinate pair (many to many)
#' and return matrix.
#'
#' @param xlon Vector of longitudes for starting coordinate pairs
#' @param xlat Vector of latitudes for starting coordinate pairs
#' @param ylon Vector of longitudes for ending coordinate pairs
#' @param ylat Vector of latitudes for ending coordinate pairs
#' @param dist_function String name of distance function: Haversine, Vincenty
#' @return Matrix of distances between each coordinate pair in meters
#' @export
dist_mtom <- function(xlon, xlat, ylon, ylat, dist_function = "Haversine") {
    .Call('_distRcpp_dist_mtom', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat, dist_function)
}

#' Compute distance between corresponding coordinate pairs in data frame.
#'
#' Compute distance between corresponding coordinate pairs and return vector.
#' For use when creating a new data frame or tbl_df column.
#'
#' @param xlon Vector of longitudes for starting coordinate pairs
#' @param xlat Vector of latitudes for starting coordinate pairs
#' @param ylon Vector of longitudes for ending coordinate pairs
#' @param ylat Vector of latitudes for ending coordinate pairs
#' @param dist_function String name of distance function: Haversine, Vincenty
#' @return Vector of distances between each coordinate pair in meters
#' @export
dist_df <- function(xlon, xlat, ylon, ylat, dist_function = "Haversine") {
    .Call('_distRcpp_dist_df', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat, dist_function)
}

#' Compute one to many distances.
#'
#' Compute distances between single starting coordinate and vector of
#' ending coordinates (one to many) and return vector.
#'
#' @param xlon Longitude for starting coordinate pair
#' @param xlat Latitude for starting coordinate pair
#' @param ylon Vector of longitudes for ending coordinate pairs
#' @param ylat Vector of latitudes for ending coordinate pairs
#' @param dist_function String name of distance function: Haversine, Vincenty
#' @return Vector of distances in meters
#' @export
dist_1tom <- function(xlon, xlat, ylon, ylat, dist_function = "Haversine") {
    .Call('_distRcpp_dist_1tom', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat, dist_function)
}

#' Compute one to one distance.
#'
#' Compute distance between two points (one to one) and return single value.
#'
#' @param xlon Longitude for starting coordinate pair
#' @param xlat Latitude for starting coordinate pair
#' @param ylon Longitude for ending coordinate pair
#' @param ylat Latitude for ending coordinate pair
#' @param dist_function String name of distance function: Haversine, Vincenty
#' @return Distance in meters
#' @export
dist_1to1 <- function(xlon, xlat, ylon, ylat, dist_function = "Haversine") {
    .Call('_distRcpp_dist_1to1', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat, dist_function)
}

#' Interpolate population/inverse-distance-weighted measures.
#'
#' Interpolate population/inverse-distance-weighted measures for each \strong{x}
#' coordinate using measures taken at surrounding \strong{y} coordinates.
#' Ending measures are double weighted by population and distance so that
#' surrounding measures taken in nearby areas and those with greater
#' populations are given more weight in final average.
#'
#' @param x_df DataFrame with coordinates that need weighted measures
#' @param y_df DataFrame with coordinates at which measures were taken
#' @param measure_col String name of measure column in y_df
#' @param x_id String name of unique identifer column in x_df
#' @param x_lon_col String name of column in x_df with longitude values
#' @param x_lat_col String name of column in x_df with latitude values
#' @param y_lon_col String name of column in y_df with longitude values
#' @param y_lat_col String name of column in y_df with latitude values
#' @param pop_col String name of column in x_df with population values
#' @param dist_function String name of distance function: "Haversine" (default) or
#' "Vincenty"
#' @param dist_transform String value of distance weight transform: "level" (default)
#' or "log"
#' @param decay Numeric value of distance weight decay: 2 (default)
#' @return Dataframe of population/distance-weighted values
#' @export
popdist_weighted_mean <- function(x_df, y_df, measure_col, x_id = "id", x_lon_col = "lon", x_lat_col = "lat", y_lon_col = "lon", y_lat_col = "lat", pop_col = "pop", dist_function = "Haversine", dist_transform = "level", decay = 2) {
    .Call('_distRcpp_popdist_weighted_mean', PACKAGE = 'distRcpp', x_df, y_df, measure_col, x_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, pop_col, dist_function, dist_transform, decay)
}

#' Interpolate inverse-distance-weighted measures.
#'
#' Interpolate inverse-distance-weighted measures for each \strong{x}
#' coordinate using measures taken at surrounding \strong{y} coordinates.
#' Ending measures are weighted by inverse distance so that
#' surrounding measures taken in nearby areas are given more weight in final
#' average.
#'
#' @param x_df DataFrame with coordinates that need weighted measures
#' @param y_df DataFrame with coordinates at which measures were taken
#' @param measure_col String name of measure column in y_df
#' @param x_id String name of unique identifer column in x_df
#' @param x_lon_col String name of column in x_df with longitude values
#' @param x_lat_col String name of column in x_df with latitude values
#' @param y_lon_col String name of column in y_df with longitude values
#' @param y_lat_col String name of column in y_df with latitude values
#' @param dist_function String name of distance function: "Haversine" (default) or
#' "Vincenty"
#' @param dist_transform String value of distance weight transform: "level" (default)
#' or "log"
#' @param decay Numeric value of distance weight decay: 2 (default)
#' @return Dataframe of distance-weighted values
#' @export
dist_weighted_mean <- function(x_df, y_df, measure_col, x_id = "id", x_lon_col = "lon", x_lat_col = "lat", y_lon_col = "lon", y_lat_col = "lat", dist_function = "Haversine", dist_transform = "level", decay = 2) {
    .Call('_distRcpp_dist_weighted_mean', PACKAGE = 'distRcpp', x_df, y_df, measure_col, x_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function, dist_transform, decay)
}

#' Find minimum distance.
#'
#' Find minimum distance between each starting point in \strong{x} and
#' possible end points, \strong{y}.
#'
#' @param x_df DataFrame with starting coordinates
#' @param y_df DataFrame with ending coordinates
#' @param x_id String name of unique identifer column in x_df
#' @param y_id String name of unique identifer column in y_df
#' @param x_lon_col String name of column in x_df with longitude values
#' @param x_lat_col String name of column in x_df with latitude values
#' @param y_lon_col String name of column in y_df with longitude values
#' @param y_lat_col String name of column in y_df with latitude values
#' @param dist_function String name of distance function: "Haversine" (default) or
#' "Vincenty"
#' @return DataFrame with id of closest point and distance in meters
#' @export
dist_min <- function(x_df, y_df, x_id = "id", y_id = "id", x_lon_col = "lon", x_lat_col = "lat", y_lon_col = "lon", y_lat_col = "lat", dist_function = "Haversine") {
    .Call('_distRcpp_dist_min', PACKAGE = 'distRcpp', x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function)
}

#' Find maximum distance.
#'
#' Find maximum distance between each starting point in \strong{x} and
#' possible end points, \strong{y}.
#'
#' @param x_df DataFrame with starting coordinates
#' @param y_df DataFrame with ending coordinates
#' @param x_id String name of unique identifer column in x_df
#' @param y_id String name of unique identifer column in y_df
#' @param x_lon_col String name of column in x_df with longitude values
#' @param x_lat_col String name of column in x_df with latitude values
#' @param y_lon_col String name of column in y_df with longitude values
#' @param y_lat_col String name of column in y_df with latitude values
#' @param dist_function String name of distance function: "Haversine" (default) or
#' "Vincenty"
#' @return DataFrame with id of farthest point and distance in meters
#' @export
dist_max <- function(x_df, y_df, x_id = "id", y_id = "id", x_lon_col = "lon", x_lat_col = "lat", y_lon_col = "lon", y_lat_col = "lat", dist_function = "Haversine") {
    .Call('_distRcpp_dist_max', PACKAGE = 'distRcpp', x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function)
}

#' Sum inverse distances.
#'
#' Find sum of inverse distances between each starting point in \strong{x}
#' and possible end points, \strong{y}.
#'
#' @param x_df DataFrame with starting coordinates
#' @param y_df DataFrame with ending coordinates
#' @param x_id String name of unique identifer column in x_df
#' @param y_id String name of unique identifer column in y_df
#' @param x_lon_col String name of column in x_df with longitude values
#' @param x_lat_col String name of column in x_df with latitude values
#' @param y_lon_col String name of column in y_df with longitude values
#' @param y_lat_col String name of column in y_df with latitude values
#' @param dist_function String name of distance function: "Haversine" (default) or
#' "Vincenty"
#' @param dist_transform String value of distance transform: "level" (default)
#' or "log"
#' @param decay Numeric value of distance weight decay: 2 (default)
#' @param scale_units Double value to divide return value by (e.g., 1000 == km)
#' @return DataFrame with sum of distances
#' @export
dist_sum_inv <- function(x_df, y_df, x_id = "id", y_id = "id", x_lon_col = "lon", x_lat_col = "lat", y_lon_col = "lon", y_lat_col = "lat", dist_function = "Haversine", dist_transform = "level", decay = 2, scale_units = 1) {
    .Call('_distRcpp_dist_sum_inv', PACKAGE = 'distRcpp', x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function, dist_transform, decay, scale_units)
}

#' Convert degrees to radians
#'
#' @param degree Degree value
#' @return Radian value (double)
#' @export
deg_to_rad <- function(degree) {
    .Call('_distRcpp_deg_to_rad', PACKAGE = 'distRcpp', degree)
}

#' Compute Haversine distance between two points
#'
#' @param xlon Longitude for starting coordinate pair
#' @param xlat Latitude for starting coordinate pair
#' @param ylon Longitude for ending coordinate pair
#' @param ylat Latitude for ending coordinate pair
#' @return Double of distance between coordinate pairs in meters
#' @export
dist_haversine <- function(xlon, xlat, ylon, ylat) {
    .Call('_distRcpp_dist_haversine', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat)
}

#' Compute Vincenty distance between two points
#'
#' @param xlon Longitude for starting coordinate pair
#' @param xlat Latitude for starting coordinate pair
#' @param ylon Longitude for ending coordinate pair
#' @param ylat Latitude for ending coordinate pair
#' @return Double of distance between coordinate pairs in meters
#' @export
dist_vincenty <- function(xlon, xlat, ylon, ylat) {
    .Call('_distRcpp_dist_vincenty', PACKAGE = 'distRcpp', xlon, xlat, ylon, ylat)
}

#' Compute inverse values from vector
#'
#' @param d Vector of values (e.g., distances)
#' @param exp Rate of decay
#' @param transform == "log" if natural log transform
#' @return Vector of inverse weights
#' @export
inverse_value <- function(d, exp, transform) {
    .Call('_distRcpp_inverse_value', PACKAGE = 'distRcpp', d, exp, transform)
}

