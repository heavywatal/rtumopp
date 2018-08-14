#' Transform coordinates into hexagonal grid
#'
#' @param .tbl data.frame with (x, y, z) columns
#' @rdname hex
#' @export
trans_coord_hex = function(.tbl) {
  if (nrow(.tbl) == 0L) return(.tbl)
  trans_coord_fcc(.tbl)
}

# face centered cubic, cubic close packed
trans_coord_fcc = function(.tbl) {
  dplyr::mutate(
    trans_coord_hex_xy(.tbl),
    x = .data$x + .data$z / sqrt(3),
    z = .data$z * sqrt(2 / 3)
  )
}

# hexagonal close packed
trans_coord_hcc = function(.tbl) {
  dplyr::mutate(
    trans_coord_hex_xy(.tbl),
    x = .data$x + ifelse(.data$z %% 2L == 1L, sqrt(3) / 3, 0),
    z = .data$z * sqrt(2 / 3)
  )
}

# 2D transformation into hexagonal grid
trans_coord_hex_xy = function(.tbl) {
  dplyr::mutate(
    .tbl,
    y = .data$y + .data$x * 0.5,
    x = .data$x * sqrt(3 / 4)
  )
}

revert_coord_hex = function(.tbl) {
  if (nrow(.tbl) == 0L) return(.tbl)
  revert_coord_fcc(.tbl)
}

revert_coord_fcc = function(.tbl) {
  dplyr::mutate(
    .tbl,
    z = round(.data$z * sqrt(3 / 2)),
    x = round(.data$x - .data$z / sqrt(3))
  ) %>% revert_coord_hex_xy()
}

revert_coord_hex_xy = function(.tbl) {
  dplyr::mutate(
    .tbl,
    x = round(.data$x * sqrt(4 / 3)),
    y = round(.data$y - .data$x * 0.5)
  )
}
