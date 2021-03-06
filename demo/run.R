library(tidyverse)
library(tumopp)

tumopp("-h")
result = tumopp("-D3 -Chex -N10000 -k1e9 -Lstep -Pmindrag")

result$population[[1]] %>%
  filter_extant() %>%
  dplyr::filter(z == 0) %>%
  plot_lattice2d()

result$population[[1]] %>%
  filter_extant() %>%
  add_surface(result$coord, result$dimensions) %>%
  dplyr::filter(z == 0) %>%
  dplyr::filter(surface) %>%
  plot_lattice2d()

if (result[["dimensions"]] > 2L) {
  library(rgl)

  if (rgl.cur()) {
    rgl.close()
  }
  rgl::open3d(windowRect = c(0, 0, 600, 600))
  result$population[[1]] %>%
    add_surface(result$coord, result$dimensions) %>%
    dplyr::filter(surface) %>%
    plot_tumor3d()
  title3d("", "", "x", "y", "z")

  .outfile = snapshot_surface(result$population[[1]])
  system(sprintf("open %s", .outfile))

  writeWebGL(".", "rgl.html", snapshot = FALSE, width = 600, height = 600)
} # fi 3D
