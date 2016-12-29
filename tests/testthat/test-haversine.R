context("Check Haversine Distances")

xlon = 86.56850
xlat = 34.78337
ylon = 86.80917
ylat = 33.50223

test_that("Haversine distance is correct", {
    expect_equal(dist_haversine(xlon, xlat, ylon, ylat), 144329.13402557538939)
})
