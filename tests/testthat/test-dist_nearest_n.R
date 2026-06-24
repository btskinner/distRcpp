context("Check nearest n distance function")

df_x = data.frame(
    id = c(100654, 100663, 100690, 100706, 100724),
    lon = c(86.56850, 86.80917, 86.17401, 86.63842, 86.29568),
    lat = c(34.78337, 33.50223, 32.36261, 34.72282, 32.36432)
)

df_y = data.frame(
    id = c(100733, 100751, 100760, 100812, 100830),
    lon = c(87.52943, 87.54577, 85.94653, 86.96514, 86.17735),
    lat = c(33.20663, 33.21440, 32.92443, 34.80562, 32.36994)
)

dff = dist_nearest_n(df_x, df_y, num_nearest = 3)

test_that("Nearest N distance function output returns data frame", {
    expect_is(dff, 'data.frame')
})

test_that("Data frame has correct number of rows", {
    expect_equal(nrow(dff), 15)
})

test_that("Nearest N distance function is correct", {
    expect_equal(dff[1,"id_x"], "100654")
    expect_equal(dff[2,"id_x"], "100654")
    expect_equal(dff[3,"id_x"], "100654")
    expect_equal(dff[1,"id_y"], "100812")
    expect_equal(dff[2,"id_y"], "100751")
    expect_equal(dff[3,"id_y"], "100733")
    expect_equal(dff[1,"meters"], 36303.105)
    expect_equal(dff[2,"meters"], 196346.573)
    expect_equal(dff[3,"meters"], 196432.676)
    expect_equal(sum(dff["meters"]), 1542764.486778593855)
})
