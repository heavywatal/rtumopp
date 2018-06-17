#' Transform coordinates into hexagonal grid
#'
#' @param .tbl data.frame with (x, y, z) columns
#' @rdname hex
#' @export
trans_coord_hex = function(.tbl) {
  trans_coord_fcc(.tbl)
}

# face centered cubic, cubic close packed
trans_coord_fcc = function(.tbl) {
  trans_coord_hex_xy(.tbl) %>%
    dplyr::mutate(x = .data$x + .data$z / sqrt(3)) %>%
    dplyr::mutate(z = .data$z * sqrt(2 / 3))
}

# hexagonal close packed
trans_coord_hcc = function(.tbl) {
  trans_coord_hex_xy(.tbl) %>%
    dplyr::mutate(x = .data$x + ifelse(.data$z %% 2L == 1L, sqrt(3) / 3, 0)) %>%
    dplyr::mutate(z = .data$z * sqrt(2 / 3))
}

# 2D transformation into hexagonal grid
trans_coord_hex_xy = function(.tbl) {
  dplyr::mutate(.tbl, y = .data$y + .data$x * 0.5) %>%
    dplyr::mutate(x = .data$x * sqrt(3 / 4))
}
