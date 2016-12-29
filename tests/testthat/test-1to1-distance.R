context("Check one to one distance function")

xlon = 86.56850
ylon = 86.80917
xlat = 34.78337
ylat = 33.50223

hav = dist_1to1(xlon, xlat, ylon, ylat, 'Haversine')
vin = dist_1to1(xlon, xlat, ylon, ylat, 'Vincenty')

test_that("One to one distance function is numeric", {
    expect_is(hav, 'numeric')
    expect_is(vin, 'numeric')
})

test_that("One to one distance function works (Haversine)", {
    expect_equal(hav, 144329.13402557538939)
})

test_that("One to one distance function works (Vincenty)", {
    expect_equal(vin, 143833.34949792452971)
})
