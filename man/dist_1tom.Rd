% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dist_1tom}
\alias{dist_1tom}
\title{Compute one to many distances.}
\usage{
dist_1tom(xlon, xlat, ylon, ylat, dist_function = "Haversine")
}
\arguments{
\item{xlon}{Longitude for starting coordinate pair}

\item{xlat}{Latitude for starting coordinate pair}

\item{ylon}{Vector of longitudes for ending coordinate pairs}

\item{ylat}{Vector of latitudes for ending coordinate pairs}

\item{dist_function}{String name of distance function: Haversine, Vincenty}
}
\value{
Vector of distances in meters
}
\description{
Compute distances between single starting coordinate and vector of
ending coordinates (one to many) and return vector.
}
