context("Check inverse value function")

x = c(100, 200, 50, 150, 200, 10, 350, 400)

test_that("Inverse value (level): exp == 1", {
    expect_equal(inverse_value(x, 1, ''), 1 / x)
    expect_equal(sum(inverse_value(x, 1, '')), 0.15202380952380953327)
})

test_that("Inverse value (level): exp == 2", {
    expect_equal(inverse_value(x, 2, ''), 1 / x^2)
    expect_equal(sum(inverse_value(x, 2, '')), 0.010608857709750566661)
})

test_that("Inverse value (log): exp == 1", {
    expect_equal(inverse_value(x, 1, 'log'), 1 / log(x))
    expect_equal(sum(inverse_value(x, 1, 'log')), 1.8217305384628765808)
})

test_that("Inverse value (log): exp == 2", {
    expect_equal(inverse_value(x, 2, 'log'), 1 / log(x)^2)
    expect_equal(sum(inverse_value(x, 2, 'log')), 0.46918109205878344437)
})
