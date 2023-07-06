tumopp::tumopp("-h")
result = tumopp::tumopp("-D3 -Chex -N10000 -k1e9 -Lstep -Pmindrag")

result$population[[1]] |>
  tumopp::filter_extant() |>
  dplyr::filter(z == 0) |>
  tumopp::plot_lattice2d()

result$population[[1]] |>
  tumopp::filter_extant() |>
  tumopp::add_surface(result$coord, result$dimensions) |>
  dplyr::filter(z == 0) |>
  dplyr::filter(surface) |>
  tumopp::plot_lattice2d()

if (result[["dimensions"]] > 2L) {
  rgl::close3d()
  rgl::open3d(windowRect = c(0, 0, 600, 600))
  result$population[[1]] |>
    tumopp::filter_extant() |>
    tumopp::add_surface(result$coord, result$dimensions) |>
    dplyr::filter(surface) |>
    tumopp::plot_tumor3d()
  rgl::title3d("", "", "x", "y", "z")

  .outfile = tumopp::snapshot_surface(result$population[[1]])
  system(sprintf("open %s", .outfile))

  rgl::writeWebGL(".", "rgl.html", snapshot = FALSE, width = 600, height = 600)
} # fi 3D
