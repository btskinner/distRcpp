![](https://img.shields.io/badge/dev-beta-red.svg) [![GitHub
release](https://img.shields.io/github/release/btskinner/distRcpp.svg)](https://github.com/btskinner/distRcpp)
[![Travis-CI Build
Status](https://travis-ci.org/btskinner/distRcpp.svg?branch=master)](https://travis-ci.org/btskinner/distRcpp)

This package uses [Rcpp](http://www.rcpp.org) to quickly compute
population/distance-weighted measures. Geodesic distances can be
computed using either
[Haversine](https://en.wikipedia.org/wiki/Haversine_formula) or
[Vincenty](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)
formulas. The package also has functions to return raw distance
measures. If you are able to [install Rcpp on your
machine](https://github.com/RcppCore/Rcpp), you should be able to
install this package and use these functions.

Install the latest development version from Github with

    devtools::install_github('btskinner/distRcpp')

**NB** This package is still in early beta stages. It does not have much
in the way of error handling. Data must be pre-processed so that no
missing (`NA`) values are given to the functions.

Available functions
-------------------

### Weighted measures

#### `dist_weighted_mean()`

Interpolate values for a vector of locations (**X**) that are the
inverse-distance-weighted average of measures taken at surrounding
locations (**Y**). For each point, *x*, nearby values of the measure
taken at **Y** are weighted more heavily than those from locations that
are farther away.

#### `popdist_weighted_mean()`

Interpolate values for a vector of locations (**X**) that are the
population/inverse-distance-weighted average of measures taken at
surrounding locations (**Y**). For each point, *x*, nearby values of the
measure taken at **Y** are weighted more heavily than those from
locations that are farther away. Measures taken in more heavily
populated *y* are given more weight than those with lower populations.
This weighting scheme is a compromise between distance and population
and is useful for interpolating measures that need to take both into
account.

### Distances

#### `dist_1to1()`

Compute and return the geodesic distance between two spatial points.
Returns distance in meters.

#### `dist_1tom()`

Compute and return the geodesic distance between one location and a
vector of other locations. Returns vector of distances in meters.

#### `dist_mtom()`

Compute and return the geodesic distance between each coordinate pair in
two vectors. Returns *n x k* matrix of distances in meters, where *n* =
\# of locations in first vector and *k* = \# of locations in second
vector.

#### `dist_df()`

Compute distance between corresponding coordinate pairs and return
vector of distances in meters. For use when creating a new `data.frame`
or [dplyr](https://CRAN.R-project.org/package=dplyr) `tbl_df()` column.

#### `dist_min()`

Compute minimum distance between each starting point, *x*, and possible
end points, **Y**. Returns vector of minimum distances in meters that
equals \# of starting points (size of **X**).

Benchmark
---------

Compare speed with base R function when measuring the distance between
every United States population-weighted county centroid as measured in
2010 (N = 3,143 with complete measurements).

### Load data

    ## libraries
    libs <- c('dplyr','microbenchmark','geosphere','distRcpp')
    lapply(libs, require, character.only = TRUE)

    ## read data
    df <- get(data(county_centers))
    df

    ## # A tibble: 3,147 × 9
    ##     fips    clon00   clat00    clon10   clat10   pclon00  pclat00   pclon10  pclat10
    ##    <chr>     <dbl>    <dbl>     <dbl>    <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
    ## 1  01001 -86.57718 32.52328 -86.64449 32.53638 -86.50183 32.50032 -86.49416 32.50039
    ## 2  01003 -87.74826 30.59278 -87.74607 30.65922 -87.76054 30.56538 -87.76238 30.54892
    ## 3  01005 -85.33131 31.85652 -85.40546 31.87067 -85.30675 31.84787 -85.31004 31.84404
    ## 4  01007 -87.12324 33.04005 -87.12715 33.01589 -87.12702 33.02595 -87.12766 33.03092
    ## 5  01009 -86.55477 33.97846 -86.56725 33.97745 -86.58262 33.96260 -86.59149 33.95524
    ## 6  01011 -85.70491 32.09828 -85.71726 32.10176 -85.70278 32.11414 -85.70119 32.11633
    ## 7  01013 -86.66223 31.73588 -86.68197 31.75167 -86.65606 31.77192 -86.65355 31.77354
    ## 8  01015 -85.81754 33.74199 -85.82251 33.77171 -85.82205 33.72213 -85.81944 33.72546
    ## 9  01017 -85.28875 32.89123 -85.39181 32.91794 -85.26586 32.86135 -85.26647 32.86044
    ## 10 01019 -85.62193 34.18416 -85.65424 34.06952 -85.62710 34.17993 -85.62919 34.17933
    ## # ... with 3,137 more rows

    ## subset to 2010 population-weighted centroids (pclon10, pclat10)
    p <- df %>% select(pclon10, pclat10) %>% na.omit %>% data.frame()

### Check for equality

    dist_R <- distm(p)
    dist_Rcpp <- dist_mtom(p[,1],p[,2],p[,1],p[,2],funname="Haversine")

    dist_R[1:5,1:5]

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]      0.0 248335.5 133369.0  83691.8 162207.0
    ## [2,] 248335.5      0.0 274424.4 282744.5 394877.3
    ## [3,] 133369.0 274424.4      0.0 215905.4 263771.5
    ## [4,]  83691.8 282744.5 215905.4      0.0 114301.5
    ## [5,] 162207.0 394877.3 263771.5 114301.5      0.0

    dist_Rcpp[1:5,1:5]

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]      0.0 248335.5 133369.0  83691.8 162207.0
    ## [2,] 248335.5      0.0 274424.4 282744.5 394877.3
    ## [3,] 133369.0 274424.4      0.0 215905.4 263771.5
    ## [4,]  83691.8 282744.5 215905.4      0.0 114301.5
    ## [5,] 162207.0 394877.3 263771.5 114301.5      0.0

    all.equal(dist_R, dist_Rcpp)

    ## [1] TRUE

### Benchmark

Mid-2012 MacBook Air, 2 GHz Intel Core i7, 8 GB 1600 MHz DDR3 SDRAM

    microbenchmark(
        dist_R = distm(p),
        dist_Rcpp = dist_mtom(p[,1],p[,2],p[,1],p[,2],funname="Haversine"),
        times = 100
    )

    ## Unit: milliseconds
    ##       expr       min        lq      mean    median        uq      max neval cld
    ##     dist_R 2549.4251 2727.2370 2848.2109 2844.3020 2924.5414 3653.536   100   b
    ##  dist_Rcpp  798.7527  819.3178  858.9428  844.8968  880.1373 1009.457   100  a
