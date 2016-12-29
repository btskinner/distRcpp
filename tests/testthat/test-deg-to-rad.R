context("Conversion of degrees to radians")

test_that("Degrees to radians", {
    expect_equal(deg_to_rad(0), 0)
    expect_equal(deg_to_rad(90), 1.5707963268)
    expect_equal(deg_to_rad(-90), -1.5707963268)
})
