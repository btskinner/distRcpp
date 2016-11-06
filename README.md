# distRcpp

![](https://img.shields.io/badge/dev-beta-red.svg)
[![GitHub release](https://img.shields.io/github/release/btskinner/distRcpp.svg)](https://github.com/btskinner/distRcpp)

This package uses [Rcpp](http://www.rcpp.org) to quickly compute population/distance-weighted measures. Geodesic distances can be computed using either [Haversine](https://en.wikipedia.org/wiki/Haversine_formula) or [Vincenty](https://en.wikipedia.org/wiki/Vincenty%27s_formulae) formulas. The package also has functions to return raw distance measures. If you are able to [install Rcpp on your machine](https://github.com/RcppCore/Rcpp), you should be able to install this package and use these functions.

Install the latest development version from Github with

```r
devtools::install_github('btskinner/distRcpp')
```

## Available functions

### Weighted measures

#### `dist_weighted_mean()`

Interpolate values for a vector of locations (**X**) that are the inverse-distance-weighted average of measures taken at surrounding locations (**Y**). For each point, *x*, nearby values of the measure taken at **Y** are weighted more heavily than those from locations that are farther away.

#### `popdist_weighted_mean()`

Interpolate values for a vector of locations (**X**) that are the population/inverse-distance-weighted average of measures taken at surrounding locations (**Y**). For each point, *x*, nearby values of the measure taken at **Y** are weighted more heavily than those from locations that are farther away. Measures taken in more heavily populated *y* are given more weight than those with lower populations. This weighting scheme is a compromise between distance and population and is useful for interpolating measures that need to take both into account.

### Distances

#### `dist_1to1()`

Compute and return the geodesic distance between two spatial points. Returns distance in meters.

#### `dist_1tom()`

Compute and return the geodesic distance between one location and a vector of other locations. Returns vector of distances in meters.

#### `dist_mtom()`

Compute and return the geodesic distance between each coordinate pair in two vectors. Returns *n x k* matrix of distances in meters, where *n* = # of locations in first vector and *k* = # of locations in second vector.

#### `dist_df()`

Compute distance between corresponding coordinate pairs and return vector of distances in meters. For use when creating a new `data.frame` or [dplyr](https://CRAN.R-project.org/package=dplyr) `tbl_df()` column.

#### `dist_min()`

Compute minimum distance between each starting point, *x*, and possible end points, **Y**. Returns vector of minimum distances in meters that equals # of starting points (size of **X**).

