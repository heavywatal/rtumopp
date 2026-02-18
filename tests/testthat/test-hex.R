test_that("hex transformation works", {
  df = data.frame
  df0 = df(x = numeric(0), y = numeric(0))
  expect_identical(trans_coord_hex(df0), df0)
  v = seq(-2, 2)
  df2 = tidyr::expand_grid(x = v, y = v, z = v)
  expect_identical(trans_coord_hex(df2), trans_coord_fcc(df2))

  expect_identical(
    trans_coord_fcc(df(x = 1, y = 0, z = 0)),
    df(x = sqrt(3) / 2, y = 0.5, z = 0)
  )
  expect_identical(
    trans_coord_fcc(df(x = 0, y = 1, z = 0)),
    df(x = 0, y = 1, z = 0)
  )
  expect_identical(
    trans_coord_fcc(df(x = 0, y = 0, z = 1)),
    df(x = 1 / sqrt(3), y = 0, z = sqrt(2 / 3))
  )
  expect_identical(
    trans_coord_fcc(df(x = 0, y = 0, z = -1)),
    df(x = -1 / sqrt(3), y = 0, z = -sqrt(2 / 3))
  )
  expect_equal(
    trans_coord_hcc(df(x = 0, y = 0, z = -1)),
    df(x = 1 / sqrt(3), y = 0, z = -sqrt(2 / 3)),
    tolerance = 1e-9
  )

  expect_identical(
    trans_coord_hex_xy(df(x = 1, y = 0)),
    df(x = sqrt(3) / 2, y = 0.5)
  )
  expect_identical(trans_coord_hex_xy(df(x = 0, y = 1)), df(x = 0, y = 1))
  expect_identical(
    trans_coord_hex_xy(df(x = 1, y = 1)),
    df(x = sqrt(3) / 2, y = 1.5)
  )

  expect_identical(revert_coord_hex(trans_coord_hex(df0)), df0)
  expect_identical(revert_coord_hex(trans_coord_hex(df2)), df2)
})
