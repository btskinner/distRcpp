context("Check Vincenty Distances")

xlon = 86.56850
xlat = 34.78337
ylon = 86.80917
ylat = 33.50223

test_that("Vincenty distance is correct", {
    expect_equal(dist_vincenty(xlon, xlat, ylon, ylat), 143833.34949792452971)
})
