test_that("multiplication works", {
  expect_identical(max_abs_xyz(data.frame(x = 2L, y = -3L, z = 1L)), 3L)
  expect_identical(dist_euclidean(data.frame(x = 1L, y = -1L, z = 1L)), sqrt(3))
  expect_equal(
    transform_rotate(data.frame(x = 1L, y = 1L, z = 1L), pi / 2, "z"),
    data.frame(x = -1, y = 1, z = 1),
    tolerance = 1e-9
  )
  expect_equal(
    transform_rotate(data.frame(x = 1L, y = 1L, z = 1L), pi / 2, "x"),
    data.frame(x = 1, y = -1, z = 1),
    tolerance = 1e-9
  )
  expect_equal(
    transform_rotate(data.frame(x = 1L, y = 1L, z = 1L), pi / 2, "y"),
    data.frame(x = 1, y = 1, z = -1),
    tolerance = 1e-9
  )
})
