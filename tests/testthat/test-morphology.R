test_that("add_surface works", {
  result = tumopp("-D3 -Chex -N256 -Llinear")
  surface_df = result$population[[1L]] |>
    filter_extant() |>
    add_surface(result$coord, result$dimensions)
  surface_df |> dplyr::count(surface) |> dplyr::pull(n) |> expect_length(2L)
})

test_that("add_phi works", {
  result = tumopp("-D3 -Cmoore -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))

  result = tumopp("-D3 -Cneumann -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))

  result = tumopp("-D3 -Chex -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))

  result = tumopp("-D2 -Cmoore -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))

  result = tumopp("-D2 -Cneumann -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))

  result = tumopp("-D2 -Chex -N32")
  pop = result$population[[1L]] |> filter_extant() |> add_phi(result$coord, result$dimensions)
  max_phi = sum(structuring_element(result$coord, result$dimensions)) - 1L
  expect_true(all(pop[["phi"]] %in% seq.int(0L, max_phi)))
})

test_that("detect_surface works", {
  cells = expand_xyz(seq.int(-1L, 1L)) |> tibble::new_tibble()
  cuboid = as_cuboid(cells, expand = 1L)
  se = structuring_element("moore", 3L)
  surface = filter_surface(cuboid, se) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(5L ** 3L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(sum(surface), sum(cuboid) - 1L)

  # mmand::erode() does not erodes from edges
  as_cuboid(cells, expand = 0L) |>
    filter_surface(se) |>
    expect_setequal(0L)

  surface_df = cells |> detect_surface(se) |>
    expect_s3_class("data.frame")
  surface_df |> dplyr::count(surface)
  cells |>
    dplyr::slice(0L) |>
    detect_surface(se) |>
    expect_s3_class("data.frame")
})

test_that("structuring_element works", {
  se = structuring_element("moore", 3L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(27L) |>
    expect_setequal(1L)
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 3L))
  expect_identical(sum(se), 27L)

  se = structuring_element("neumann", 3L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(27L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 3L))
  expect_identical(sum(se), 7L)

  se = structuring_element("hex", 3L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(27L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 3L))
  expect_identical(sum(se), 13L)

  se = structuring_element("moore", 2L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(9L) |>
    expect_setequal(1L)
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 1L))
  expect_identical(sum(se), 9L)

  se = structuring_element("neumann", 2L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(9L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 1L))
  expect_identical(sum(se), 5L)

  se = structuring_element("hex", 2L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_type("integer") |>
    expect_length(9L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(dim(se), c(x = 3L, y = 3L, z = 1L))
  expect_identical(sum(se), 7L)
})

test_that("cuboid class works", {
  grid = expand_xyz(seq_len(4L), seq_len(3L), seq_len(2L)) |>
    expect_s3_class("data.frame")
  expect_identical(dim(grid), c(24L, 3L))

  cuboid = as_cuboid(grid) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_length(24L) |>
    expect_setequal(1L)
  expect_identical(sum(cuboid), 24L)

  expanded = as_cuboid(grid, expand = 1L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_length(120L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(sum(expanded), sum(cuboid))
  center = expanded[seq_len(4L) + 1L, seq_len(3L) + 1L, seq_len(2L) + 1L]
  expect_identical(center, cuboid, ignore_attr = TRUE)

  cuboid_df = expanded |> as.data.frame() |>
    expect_s3_class("data.frame")
  expect_identical(dim(cuboid_df), c(120L, 4L))
  expect_identical(cuboid_df |> as.array(), expanded)
  expect_identical(cuboid_df |> as.data.frame(), cuboid_df)
  cuboid_df |> dplyr::filter(state > 0L) |> dplyr::select(x:z) |>
    expect_identical(grid, ignore_attr = TRUE)
})

test_that("cuboid class 2D works", {
  grid = expand_xyz(seq_len(4L), seq_len(3L), 0L) |>
    expect_s3_class("data.frame")
  expect_identical(dim(grid), c(12L, 3L))
  cuboid = as_cuboid(grid) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_length(12L) |>
    expect_setequal(1L)
  expect_identical(sum(cuboid), 12L)

  expanded = as_cuboid(grid, expand = 1L) |>
    expect_s3_class(c("cuboid", "array")) |>
    expect_length(30L) |>
    expect_setequal(c(0L, 1L))
  expect_identical(sum(expanded), sum(cuboid))
  center = expanded[seq_len(4L) + 1L, seq_len(3L) + 1L, 1L, drop = FALSE]
  expect_identical(center, cuboid, ignore_attr = TRUE)

  cuboid_df = expanded |> as.data.frame() |>
    expect_s3_class("data.frame")
  expect_identical(dim(cuboid_df), c(30L, 4L))
  expect_identical(cuboid_df |> as.array(), expanded)
  expect_identical(cuboid_df |> as.data.frame(), cuboid_df)
  cuboid_df |> dplyr::filter(state > 0L) |> dplyr::select(x:z) |>
    expect_identical(grid, ignore_attr = TRUE)
})
