context("Check many to many distance functions")

df = data.frame(
    id = c(100654, 100663, 100690, 100706, 100724,
           100733, 100751, 100760, 100812, 100830),
    lon = c(86.56850, 86.80917, 86.17401, 86.63842, 86.29568,
            87.52943, 87.54577, 85.94653, 86.96514, 86.17735),
    lat = c(34.78337, 33.50223, 32.36261, 34.72282, 32.36432,
            33.20663, 33.21440, 32.92443, 34.80562, 32.36994)
)

hav_mat = dist_mtom(df$lon, df$lat, df$lon, df$lat, 'Haversine')
vin_mat = dist_mtom(df$lon, df$lat, df$lon, df$lat, 'Vincenty')

test_that("Many to many distance function is matrix", {
    expect_is(hav_mat, 'matrix')
    expect_is(vin_mat, 'matrix')
})

test_that("Many to many distance matrix diagonal is 0", {
    expect_equal(diag(hav_mat), rep(0,10))
    expect_equal(diag(vin_mat), rep(0,10))
})

test_that("Many to many distance function works (Haversine)", {
    expect_equal(hav_mat[1,2], 144329.13402557538939)
    expect_equal(hav_mat[4,8], 210170.68224844287033)
    expect_equal(sum(hav_mat), 13759106.810155073181)
})

test_that("Many to many distance function works (Vincenty)", {
    expect_equal(vin_mat[1,2], 143833.34949792452971)
    expect_equal(vin_mat[4,8], 209505.45701617305167)
    expect_equal(sum(vin_mat), 13724008.648433696479)
})
