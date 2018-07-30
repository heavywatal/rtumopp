#' Utilities for coordinate values
#'
#' `max_abs_xyz` extracts max values for plot limits.
#' @param .tbl data.frame with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @rdname coord
#' @export
#' @examples
#' max_abs_xyz(data.frame(x=2, y=-3, z=4))
max_abs_xyz = function(.tbl) {
  max(abs(.tbl[c("x", "y", "z")]))
}

#' @description
#' `dist_euclidean` calculates distance from a specified cell.
#' @param point named vector or tibble
#' @return numeric vector
#' @rdname coord
#' @export
dist_euclidean = function(.tbl, point = c(x = 0, y = 0, z = 0)) {
  sqrt((.tbl[["x"]] - point[["x"]])**2 + (.tbl[["y"]] - point[["y"]])**2 + (.tbl[["z"]] - point[["z"]])**2)
}

#' @description
#' `rotate` modifies coordinates centering on a specified axis.
#' @param theta radian angle
#' @param axis a string
#' @return modified data.frame
#' @rdname coord
#' @export
rotate = function(.tbl, theta, axis = c("z", "x", "y")) {
  axis = match.arg(axis)
  .x = .tbl[["x"]]
  .y = .tbl[["y"]]
  .z = .tbl[["z"]]
  .sin = sin(theta)
  .cos = cos(theta)
  if (axis == "z") {
    dplyr::mutate(
      .tbl,
      x = .x * .cos - .y * .sin,
      y = .x * .sin + .y * .cos
    )
  } else if (axis == "x") {
    dplyr::mutate(
      .tbl,
      y = .y * .cos - .z * .sin,
      z = .y * .sin + .z * .cos
    )
  } else { # y
    dplyr::mutate(
      .tbl,
      x = .x * .cos + .z * .sin,
      z = -.x * .sin + .z * .cos
    )
  }
}
