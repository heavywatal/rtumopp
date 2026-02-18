#' 3D plotting with rgl
#'
#' @details
#' `plot_tumor3d` plots tumor in 3D with rgl.
#' @param .tbl data.frame with (x, y, z, col)
#' @param limits passed to xlim, ylim, zlim of [rgl::plot3d()]
#' @rdname plot-rgl
#' @export
plot_tumor3d = function(.tbl = NULL, limits = NULL) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("ERROR: rgl is not installed", call. = FALSE)
  }
  on.exit({
    rgl::box3d()
    rgl::view3d(15, 15, 15, 0.9)
  })
  if (nrow(.tbl) == 0L) {
    return(invisible(.tbl))
  }
  ..col = purrr::pluck(
    .tbl,
    "col",
    .default = getOption("tumopp.default_color", "#666666")
  )
  rgl::plot3d(
    .tbl$x,
    .tbl$y,
    .tbl$z,
    xlab = "",
    ylab = "",
    zlab = "",
    axes = FALSE,
    type = "s",
    col = ..col,
    alpha = 1,
    radius = 1,
    aspect = TRUE,
    xlim = limits,
    ylim = limits,
    zlim = limits
  )
}

#' @details
#' `snapshot_surface` is a shortcut of [plot_tumor3d()] and [rgl::snapshot3d()].
#' @param filename string
#' @param ... passed to [plot_tumor3d()].
#' @rdname plot-rgl
#' @export
snapshot_surface = function(.tbl, filename = tempfile("rgl_", fileext = ".png"), ...) {
  on.exit(rgl::close3d())
  rgl::open3d()
  dplyr::filter(.tbl, .data$surface) |> plot_tumor3d(...)
  rgl::snapshot3d(filename)
  filename
}

#' @details
#' `add_col` adds a column for color in [plot_tumor3d()].
#' @param column column name to colorcode
#' @param palette name of palette in grDevices or viridis.
#' @rdname plot-rgl
#' @export
add_col = function(.tbl, column = "clade", palette = "turbo", ...) {
  .column = .tbl[[column]]
  if (!is.factor(.column)) {
    .column = as.factor(.column)
  }
  .levels = levels(.column)
  n = length(.levels)
  .map = if (n > 1L) {
    palette_colors(n, palette = palette, ...)
  } else {
    getOption("tumopp.default_color", "#666666")
  }
  dplyr::mutate(.tbl, col = .map[.column])
}

palette_colors = function(n = NULL, palette = "Okabe-Ito", ...) {
  tryCatch(
    grDevices::palette.colors(n + 1, palette, ...)[-1],
    error = \(e) scales::pal_viridis(option = palette, ...)(n)
  )
}
