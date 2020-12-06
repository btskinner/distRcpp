![](https://img.shields.io/badge/dev-beta-red.svg) [![GitHub
release](https://img.shields.io/github/release/btskinner/distRcpp.svg)](https://github.com/btskinner/distRcpp)
[![R build
status](https://github.com/btskinner/distRcpp/workflows/R-CMD-check/badge.svg)](https://github.com/btskinner/distRcpp/actions)

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

**NB** This package is still in early beta stages. It does not have much
in the way of error handling. Data must be pre-processed so that no
missing (`NA`) values are given to the functions.

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

## Benchmark

Compare speed with base R function when measuring the distance between
every United States population-weighted county centroid as measured in
2010 (N = 3,143 with complete measurements).

### Load data

    ## libraries
    libs <- c("tidyverse","microbenchmark","geosphere","distRcpp")
    sapply(libs, require, character.only = TRUE)

    ## read data
    df <- get(data(county_centers))
    df

    ## # A tibble: 3,147 x 9
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
    ## # … with 3,137 more rows

    ## subset to 2010 population-weighted centroids (pclon10, pclat10)
    p <- df %>% select(pclon10, pclat10) %>% drop_na() %>% data.frame()

### Check for equality

    dist_R <- geosphere::distm(p, p, fun = distHaversine)
    dist_Rcpp <- distRcpp::dist_mtom(p[,1],p[,2],p[,1],p[,2])

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

2018 MacBookPro, 2.9 GHz Intel Core i9, 32 GB 2400 MHz DDR4 SDRAM

    microbenchmark(
        dist_R = geosphere::distm(p, p, fun = distHaversine),
        dist_Rcpp = distRcpp::dist_mtom(p[,1],p[,2],p[,1],p[,2]),
        times = 100
    )

    ## Unit: milliseconds
    ##       expr       min        lq      mean    median        uq       max neval
    ##     dist_R 1930.1767 2047.0935 2228.0031 2234.7130 2361.3778 2673.8354   100
    ##  dist_Rcpp  466.6231  480.8235  511.5361  501.9158  526.1218  684.4163   100

### Big file

    ## get census block group centers of population
    bg <- readr::read_csv("https://www2.census.gov/geo/docs/reference/cenpop2010/blkgrp/CenPop2010_Mean_BG.txt") %>%
        setNames(tolower(names(.))) %>%
        filter(statefp < 56) %>%
        mutate(id = paste0(statefp, countyfp, tractce, blkgrpce),
               lon = longitude,
               lat = latitude) %>%
        select(id, lon, lat) %>%
        drop_na()

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────
    ## cols(
    ##   STATEFP = col_character(),
    ##   COUNTYFP = col_character(),
    ##   TRACTCE = col_character(),
    ##   BLKGRPCE = col_double(),
    ##   POPULATION = col_double(),
    ##   LATITUDE = col_double(),
    ##   LONGITUDE = col_double()
    ## )

    ct <- get(data(county_centers)) %>%
        rename(id = fips,
               lon = pclon10,
               lat = pclat10) %>%
        drop_na()
    bg

    ## # A tibble: 217,330 x 3
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
    ## # … with 217,320 more rows

    ct

    ## # A tibble: 3,137 x 9
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
    ## # … with 3,127 more rows

    system.time(dist_Rcpp <- distRcpp::dist_min(x_df = ct, y_df = bg))

    ##    user  system elapsed 
    ##  34.304   1.526  35.888
