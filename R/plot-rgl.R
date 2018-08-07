#' 3D plotting with rgl
#'
#' @description
#' `plot_tumor3d` plots tumor in 3D with rgl.
#' @param .tbl data.frame with (x, y, z, col)
#' @param limits passed to xlim, ylim, zlim of `rgl::plot3d()`
#' @rdname plot-rgl
#' @export
plot_tumor3d = function(.tbl = NULL, limits = NULL) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("ERROR: rgl is not installed")
  }
  on.exit({
    rgl::box3d()
    rgl::view3d(15, 15, 15, 0.9)
  })
  if (NROW(.tbl) == 0L) return(invisible(.tbl))
  col = purrr::pluck(.tbl, "col", .default = "#333333")
  rgl::plot3d(
    .tbl$x, .tbl$y, .tbl$z,
    xlab = "", ylab = "", zlab = "", axes = FALSE,
    type = "s", col = col, alpha = 1, radius = 1, aspect = TRUE,
    xlim = limits, ylim = limits, zlim = limits
  )
}

#' @description
#' `snapshot_surface` is a shortcut of `plot_tumor3d()` and `rgl::snapshot3d()`.
#' @param filename string
#' @param ... passed to `plot_tumor3d()`
#' @rdname plot-rgl
#' @export
snapshot_surface = function(.tbl, filename = tempfile("rgl_", fileext = ".png"), ...) {
  on.exit(rgl::rgl.close())
  rgl::open3d(useNULL = FALSE)
  dplyr::filter(.tbl, .data$surface) %>% plot_tumor3d(...)
  rgl::snapshot3d(filename)
  filename
}

#' @description
#' `add_col` adds a column for color in `plot_tumor3d()`.
#' @param column column name to colorcode
#' @param palette name for RColorBrewer::brewer.pal()
#' @param reverse logical for order of color vector
#' @rdname plot-rgl
#' @export
add_col = function(.tbl, column="clade", palette="Spectral", reverse=FALSE) {
  .column = .tbl[[column]]
  if (!is.factor(.column)) .column = as.factor(.column)
  .levels = levels(.column)
  .map = wtl::brewer_palette(palette, length(.levels))
  if (reverse) .map = reverse(.map)
  dplyr::mutate(.tbl, col = .map[.column])
}
