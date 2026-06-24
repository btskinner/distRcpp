![](https://img.shields.io/badge/dev-beta-red.svg) [![GitHub
release](https://img.shields.io/github/release/btskinner/distRcpp.svg)](https://github.com/btskinner/distRcpp)
[![R build
status](https://github.com/btskinner/distRcpp/actions/workflows/check-standard.yml/badge.svg)](https://github.com/btskinner/distRcpp/actions)

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

    devtools::install_github("btskinner/distRcpp")

**NB** This package is still in beta stages. It does not have much in
the way of error handling. Data must be pre-processed so that no missing
(`NA`) values are given to the functions.

## Available functions

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

#### `dist_max()`

Compute maximum distance between each starting point, *x*, and possible
end points, **Y**. Returns vector of maximum distances in meters that
equals \# of starting points (size of **X**).

#### `dist_nearest_n()`

Compute the distance between each starting point, *x*, and possible end
points, **Y** and return the ids and distances to the *n* nearest points
(default = 10).

## Benchmark

Compare speed with base R function when measuring the distance between
every United States population-weighted county centroid as measured in
2010 (N = 3,143 with complete measurements).

### Load data

    ## libraries
    libs <- c("tidyverse","microbenchmark","geosphere","distRcpp")
    sapply(libs, require, character.only = TRUE)

    ## read data
    df <- get(data(countycentroids))
    df

    ## # A tibble: 3,147 × 9
    ##    fips  clon00 clat00 clon10 clat10 pclon00 pclat00 pclon10 pclat10
    ##    <chr>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 01001  -86.6   32.5  -86.6   32.5   -86.5    32.5   -86.5    32.5
    ##  2 01003  -87.7   30.6  -87.7   30.7   -87.8    30.6   -87.8    30.5
    ##  3 01005  -85.3   31.9  -85.4   31.9   -85.3    31.8   -85.3    31.8
    ##  4 01007  -87.1   33.0  -87.1   33.0   -87.1    33.0   -87.1    33.0
    ##  5 01009  -86.6   34.0  -86.6   34.0   -86.6    34.0   -86.6    34.0
    ##  6 01011  -85.7   32.1  -85.7   32.1   -85.7    32.1   -85.7    32.1
    ##  7 01013  -86.7   31.7  -86.7   31.8   -86.7    31.8   -86.7    31.8
    ##  8 01015  -85.8   33.7  -85.8   33.8   -85.8    33.7   -85.8    33.7
    ##  9 01017  -85.3   32.9  -85.4   32.9   -85.3    32.9   -85.3    32.9
    ## 10 01019  -85.6   34.2  -85.7   34.1   -85.6    34.2   -85.6    34.2
    ## # ℹ 3,137 more rows

    ## subset to 2010 population-weighted centroids (pclon10, pclat10)
    p <- df |> select(pclon10, pclat10) |> drop_na() |> data.frame()

### Check for equality

    dist_R <- matrix(NA_real_, 5, 5)
    px <- p[1:5, ]
    py <- p[6:10,]
    for (i in 1:nrow(px)) {
        for (j in 1:nrow(py)) {
            dist_R[i,j] <- geosphere::distHaversine(px[i,], py[j,],
                r = 6371008.7714150594547
            )
        }
    }
    dist_Rcpp <- distRcpp::dist_mtom(px[,1],px[,2],py[,1],py[,2])

    dist_R

    ##           [,1]     [,2]      [,3]     [,4]      [,5]
    ## [1,]  85892.47  82203.6 150016.92 121676.7 203245.07
    ## [2,] 262098.33 172259.4 397755.23 348992.9 450634.38
    ## [3,]  47726.63 127196.0 214555.42 113093.1 261372.69
    ## [4,] 167952.38 146728.2 143942.96 174702.7 188582.34
    ## [5,] 220675.60 242664.0  75744.09 173047.7  92074.26

    dist_Rcpp

    ##           [,1]     [,2]      [,3]     [,4]      [,5]
    ## [1,]  85892.47  82203.6 150016.92 121676.7 203245.07
    ## [2,] 262098.33 172259.4 397755.23 348992.9 450634.38
    ## [3,]  47726.63 127196.0 214555.42 113093.1 261372.69
    ## [4,] 167952.38 146728.2 143942.96 174702.7 188582.34
    ## [5,] 220675.60 242664.0  75744.09 173047.7  92074.26

    all.equal(dist_R, dist_Rcpp)

    ## [1] TRUE

### Benchmark

2024 MacBookPro, Apple M4 Pro, 48 GB Memory

    x <- 1:1000
    y <- 1001:2000
    microbenchmark(
        dist_R = geosphere::distm(p[x,], p[y,], fun = distHaversine),
        dist_Rcpp = distRcpp::dist_mtom(p[x,1],p[x,2],p[y,1],p[y,2]),
        times = 100
    )

    ## Unit: milliseconds
    ##       expr       min        lq      mean    median        uq      max neval
    ##     dist_R 172.34071 180.40715 194.31214 187.11603 201.72851 274.3383   100
    ##  dist_Rcpp  42.59912  42.90125  44.57169  43.26931  44.65128  93.9795   100

### Big file

    ## get census block group centers of population
    file_url <- file.path(
        "https://www2.census.gov/geo/docs/reference",
        "cenpop2010/blkgrp/CenPop2010_Mean_BG.txt"
    )
    bg <- readr::read_csv(file_url, show_col_types = FALSE) |> 
        rename_all(tolower) |> 
        filter(statefp < 56) |> 
        mutate(id = paste0(statefp, countyfp, tractce, blkgrpce),
               lon = longitude,
               lat = latitude) |> 
        select(id, lon, lat) |> 
        drop_na()

    ct <- get(data(countycentroids)) |> 
        rename(id = fips,
               lon = pclon10,
               lat = pclat10) |> 
        drop_na()
    bg

    ## # A tibble: 217,330 × 3
    ##    id             lon   lat
    ##    <chr>        <dbl> <dbl>
    ##  1 010010201001 -86.5  32.5
    ##  2 010010201002 -86.5  32.5
    ##  3 010010202001 -86.5  32.5
    ##  4 010010202002 -86.5  32.5
    ##  5 010010203001 -86.5  32.5
    ##  6 010010203002 -86.5  32.5
    ##  7 010010204001 -86.4  32.5
    ##  8 010010204002 -86.4  32.5
    ##  9 010010204003 -86.4  32.5
    ## 10 010010204004 -86.4  32.5
    ## # ℹ 217,320 more rows

    ct

    ## # A tibble: 3,137 × 9
    ##    id    clon00 clat00 clon10 clat10 pclon00 pclat00   lon   lat
    ##    <chr>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl> <dbl>
    ##  1 01001  -86.6   32.5  -86.6   32.5   -86.5    32.5 -86.5  32.5
    ##  2 01003  -87.7   30.6  -87.7   30.7   -87.8    30.6 -87.8  30.5
    ##  3 01005  -85.3   31.9  -85.4   31.9   -85.3    31.8 -85.3  31.8
    ##  4 01007  -87.1   33.0  -87.1   33.0   -87.1    33.0 -87.1  33.0
    ##  5 01009  -86.6   34.0  -86.6   34.0   -86.6    34.0 -86.6  34.0
    ##  6 01011  -85.7   32.1  -85.7   32.1   -85.7    32.1 -85.7  32.1
    ##  7 01013  -86.7   31.7  -86.7   31.8   -86.7    31.8 -86.7  31.8
    ##  8 01015  -85.8   33.7  -85.8   33.8   -85.8    33.7 -85.8  33.7
    ##  9 01017  -85.3   32.9  -85.4   32.9   -85.3    32.9 -85.3  32.9
    ## 10 01019  -85.6   34.2  -85.7   34.1   -85.6    34.2 -85.6  34.2
    ## # ℹ 3,127 more rows

    system.time(dist_Rcpp <- distRcpp::dist_min(x_df = ct, y_df = bg))

    ##    user  system elapsed 
    ##  27.978   0.046  28.027

    dist_Rcpp |> tibble()

    ## # A tibble: 3,137 × 3
    ##    id_start id_end       meters
    ##    <chr>    <chr>         <dbl>
    ##  1 01001    010010201002  2114.
    ##  2 01003    010030109051  1872.
    ##  3 01005    010059505002  6413.
    ##  4 01007    010070100012  4920.
    ##  5 01009    010090502002  1965.
    ##  6 01011    010119522002  1662.
    ##  7 01013    010139531003  2417.
    ##  8 01015    010150007003  1304.
    ##  9 01017    010179542001  1397.
    ## 10 01019    010199560001  2065.
    ## # ℹ 3,127 more rows
