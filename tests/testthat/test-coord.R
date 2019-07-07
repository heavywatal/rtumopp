test_that("multiplication works", {
  expect_equal(max_abs_xyz(data.frame(x = 2, y = -3, z = 1)), 3)
  expect_equal(dist_euclidean(data.frame(x = 1, y = -1, z = 1)), sqrt(3))
  expect_equal(
    transform_rotate(data.frame(x = 1, y = 1, z = 1), pi / 2, "z"),
    data.frame(x = -1, y = 1, z = 1)
  )
  expect_equal(
    transform_rotate(data.frame(x = 1, y = 1, z = 1), pi / 2, "x"),
    data.frame(x = 1, y = -1, z = 1)
  )
  expect_equal(
    transform_rotate(data.frame(x = 1, y = 1, z = 1), pi / 2, "y"),
    data.frame(x = 1, y = 1, z = -1)
  )
})
