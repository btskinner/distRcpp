context("Check data frame distance functions")

df = data.frame(
    xid = c(100654, 100663, 100690, 100706, 100724),
    yid = c(100733, 100751, 100760, 100812, 100830),
    xlon = c(86.56850, 86.80917, 86.17401, 86.63842, 86.29568),
    xlat = c(34.78337, 33.50223, 32.36261, 34.72282, 32.36432),
    ylon = c(87.52943, 87.54577, 85.94653, 86.96514, 86.17735),
    ylat = c(33.20663, 33.21440, 32.92443, 34.80562, 32.36994)
)

df['dist_hav'] = dist_df(df$xlon, df$xlat, df$ylon, df$ylat, 'Haversine')
df['dist_vin'] = dist_df(df$xlon, df$xlat, df$ylon, df$ylat, 'Vincenty')

test_that("Data frame distance function is data frame", {
    expect_is(df['dist_hav'], 'data.frame')
    expect_is(df['dist_vin'], 'data.frame')
})

test_that("Data frame distance function works (Haversine)", {
    expect_equal(df[1,'dist_hav'], 196652.45572134392569)
    expect_equal(df[2,'dist_hav'], 75612.702710177545669)
    expect_equal(sum(df['dist_hav']), 380752.94093842280563)
})

test_that("Data frame distance function works (Vincenty)", {
    expect_equal(df[1,'dist_vin'], 196135.76363466429757)
    expect_equal(df[2,'dist_vin'], 75625.862247905286495)
    expect_equal(sum(df['dist_vin']), 380064.45919823565055)
})
