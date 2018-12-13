#' Utilities for coordinate values
#'
#' `max_abs_xyz` extracts max values for plot limits.
#' @param .tbl data.frame with (x, y, z) columns
#' @rdname coord
#' @export
#' @examples
#' max_abs_xyz(data.frame(x = 2, y = -3, z = 4))
max_abs_xyz = function(.tbl) {
  max(abs(.tbl[c("x", "y", "z")]))
}

#' @details
#' `dist_euclidean` calculates distance from a specified cell.
#' @param point named vector or tibble
#' @rdname coord
#' @export
dist_euclidean = function(.tbl, point = c(x = 0, y = 0, z = 0)) {
  sqrt((.tbl[["x"]] - point[["x"]])**2 + (.tbl[["y"]] - point[["y"]])**2 + (.tbl[["z"]] - point[["z"]])**2)
}

#' @details
#' `rotate` modifies coordinates centering on a specified axis.
#' @param theta radian angle
#' @param axis a string
#' @rdname coord
#' @export
transform_rotate = function(.tbl, theta, axis = c("z", "x", "y")) {
  axis = match.arg(axis)
  name1 = switch(axis, x = "y", y = "z", z = "x")
  name2 = switch(axis, x = "z", y = "x", z = "y")
  v1 = .tbl[[name1]]
  v2 = .tbl[[name2]]
  sin_theta = sin(theta)
  cos_theta = cos(theta)
  dplyr::mutate(
    .tbl,
    !!name1 := v1 * cos_theta - v2 * sin_theta,
    !!name2 := v1 * sin_theta + v2 * cos_theta
  )
}
