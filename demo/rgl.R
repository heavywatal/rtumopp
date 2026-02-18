# https://github.com/dmurdoch/rgl/issues/488
options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)

## transformation

rgl::close3d()
rgl::open3d(windowRect = c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
rgl::axes3d()
.r = 2
tibble::tibble(x = seq(-.r, .r), y = x, z = x) |>
  tidyr::expand(x, y, z) |>
  tumopp:::trans_coord_hcc() |>
  # tumopp:::trans_coord_fcc() |>
  dplyr::mutate(r = sqrt(x * x + y * y + z * z)) |>
  with(rgl::spheres3d(x, y, z, color = "#009999", radius = 0.51, alpha = 0.6))
rgl::title3d("", "", "x", "y", "z")

######## 1#########2#########3#########4#########5#########6#########7#########
## minimum

.hex_xy = readr::read_csv(I("x,y,z\n0,0,0\n1,0,0\n0,1,0\n1,0,-1"))

rgl::close3d()
rgl::open3d(windowRect = c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
rgl::axes3d()
.hex_xy |>
  tumopp::trans_coord_hex() |>
  with(rgl::spheres3d(x, y, z, color = "#009999", radius = 0.51, alpha = 0.6))
rgl::title3d("", "", "x", "y", "z")

######## 1#########2#########3#########4#########5#########6#########7#########
## neighbors

.hex_xy = readr::read_csv(I(
  "x,y,z\n0,0,0\n0,1,0\n0,-1,0\n-1,0,0\n-1,1,0\n1,0,0\n1,-1,0"
))

rgl::close3d()
rgl::open3d(windowRect = c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
rgl::axes3d()
.hex_xy |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::filter(x > 0 | (x == 0 & y == 0)) |> dplyr::mutate(z = -1)
  ) |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::filter(x < 0 | (x == 0 & y == 0)) |> dplyr::mutate(z = 1)
  ) |>
  print() |>
  tumopp:::trans_coord_fcc() |>
  with(rgl::spheres3d(x, y, z, color = "#009999", radius = 0.51, alpha = 0.6))
rgl::title3d("", "", "x", "y", "z")


rgl::close3d()
rgl::open3d(windowRect = c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(40, 20, 60)
rgl::axes3d()
.hex_xy |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::mutate(z = -1) |> dplyr::filter(x < 0 | (x == 0 & y == 0))
  ) |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::mutate(z = 1) |> dplyr::filter(x < 0 | (x == 0 & y == 0))
  ) |>
  print() |>
  dplyr::bind_rows(.hex_xy |> dplyr::mutate(z = 5)) |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::mutate(z = 4) |> dplyr::filter(x > 0 | (x == 0 & y == 0))
  ) |>
  dplyr::bind_rows(
    .hex_xy |> dplyr::mutate(z = 6) |> dplyr::filter(x > 0 | (x == 0 & y == 0))
  ) |>
  tumopp:::trans_coord_hcc() |>
  with(rgl::spheres3d(x, y, z, color = "#009999", radius = 0.51, alpha = 0.6))
rgl::title3d("", "", "x", "y", "z")
