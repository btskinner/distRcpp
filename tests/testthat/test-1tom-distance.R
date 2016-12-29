context("Check one to many distance functions")

df = data.frame(
    id = c(100654, 100663, 100690, 100706, 100724,
           100733, 100751, 100760, 100812, 100830),
    lon = c(86.56850, 86.80917, 86.17401, 86.63842, 86.29568,
            87.52943, 87.54577, 85.94653, 86.96514, 86.17735),
    lat = c(34.78337, 33.50223, 32.36261, 34.72282, 32.36432,
            33.20663, 33.21440, 32.92443, 34.80562, 32.36994)
)

hav_vec = dist_1tom(df$lon[1], df$lat[1], df$lon, df$lat, 'Haversine')
vin_vec = dist_1tom(df$lon[1], df$lat[1], df$lon, df$lat, 'Vincenty')

test_that("One to many distance function is vector of length 10", {
    expect_equal(length(hav_vec), 10)
    expect_equal(length(vin_vec), 10)
})

test_that("One to many distance function works (Haversine)", {
    expect_equal(hav_vec[1], 0)
    expect_equal(hav_vec[2], 144329.1340255753893871)
    expect_equal(hav_vec[4], 9291.3474984819404199)
    expect_equal(sum(hav_vec), 1611479.2332509944681)
})

test_that("One to many distance function works (Vincenty)", {
    expect_equal(vin_vec[1], 0)
    expect_equal(vin_vec[2], 143833.3494979245297145)
    expect_equal(vin_vec[4], 9279.3240185879749333)
    expect_equal(sum(vin_vec), 1606383.1903742046561)
})
