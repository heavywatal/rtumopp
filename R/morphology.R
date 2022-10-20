#' Extract surface cells with mathematical morphology
#'
#' `add_surface()` adds a binary column to the data.frame.
#' @param population a data.frame with (x, y, z)
#' @param coord a string
#' @param dimensions an integer
#' @rdname morphology
#' @export
add_surface = function(population, coord, dimensions) {
  extant = filter_extant(population)
  kernel = structuring_element(coord, dimensions)
  col_surface = detect_surface(extant, kernel) |>
    dplyr::select("id", "surface")
  population |>
    dplyr::left_join(col_surface, by = "id")
}

#' @description
#' `add_phi()` counts empty neighbors of each cell and adds an integer column.
#' @rdname morphology
#' @export
add_phi = function(population, coord, dimensions) {
  extant = filter_extant(population)
  cells = extant[c("x", "y", "z", "id")]
  if (!is.integer(cells$x)) cells = revert_coord_hex(cells)
  kernel_df = structuring_element(coord, dimensions) |>
    as.data.frame() |>
    dplyr::filter(.data$state > 0L)
  counts = purrr::pmap_dfr(kernel_df, ~ {
    dplyr::mutate(cells, x = .data$x + ..1, y = .data$y + ..2, z = .data$z + ..3)
  }) |>
    dplyr::count(.data$x, .data$y, .data$z)
  info = cells |>
    dplyr::left_join(counts, by = c("x", "y", "z")) |>
    dplyr::transmute(.data$id, phi = nrow(kernel_df) - .data$n)
  population |>
    dplyr::left_join(info, by = "id")
}

detect_surface = function(.tbl, se) {
  if (!is.integer(.tbl$x)) .tbl = revert_coord_hex(.tbl)
  product = as_cuboid(.tbl, expand = 1L) |>
    filter_surface(se) |>
    as.data.frame() |>
    dplyr::mutate(surface = .data$state > 0L, state = NULL)
  dplyr::left_join(.tbl, product, by = c("x", "y", "z"))
}

# Filter img by surface
filter_surface = function(x, kernel) {
  y = x - as.integer(mmand::erode(x, kernel))
  attributes(y) = attributes(x)
  y
}

# Construct a structuring element (kernel)
structuring_element = function(coord = c("moore", "neumann", "hex"), dimensions = 3L) {
  coord = match.arg(coord)
  df = expand_xyz(c(-1L, 0L, 1L))
  if (coord == "neumann") {
    df = dplyr::filter(df, abs(.data$x) + abs(.data$y) + abs(.data$z) < 2)
  } else if (coord == "hex") {
    df = dplyr::filter(df, (trans_coord_hex(df) |> dist_euclidean()) < 1.1)
  }
  if (dimensions < 3L) {
    df = dplyr::filter(df, .data$z == 0L)
  }
  as_cuboid(df)
}

#' Binary expression of cell existence
#'
#' `as_cuboid()` expand a coordinate set to a binary cuboid with zero-padding.
#' The data.frame representation is ordered in the same way as the array representation.
#' They have some extra attributes to ensure round-trip conversion.
#' @param x a data.frame with (x, y, z) or 3D array.
#' @param expand an integer.
#' @rdname cuboid
#' @export
as_cuboid = function(x, expand = 0L) {
  if (nrow(x) == 0L) {
    res = array(integer(0L), c(0L, 0L, 0L))
    attr(res, "start") = c(c(0L, 0L, 0L))
    class(res) = c("cuboid", "array")
    return(res)
  }
  .tbl = x[c("x", "y", "z")]
  if (all(.tbl[["z"]] == 0L)) {
    expand = c(expand, expand, 0L)
  }
  start = purrr::map_int(.tbl, min) - expand
  .dim = purrr::map_int(.tbl, max) - start + 1L + expand
  idx = purrr::map2(.tbl, start, ~.x - .y)
  serial = 1L + idx[["x"]] + .dim[1L] * idx[["y"]] + .dim[1L] * .dim[2L] * idx[["z"]]
  res = array(integer(prod(.dim)), .dim)
  res[serial] = 1L
  attr(res, "start") = start
  class(res) = c("cuboid", "array")
  res
}

# Convert binary array to data.frame with (x, y, z) columns
#' @param row.names,optional,... ignored; to suppress warnings for S3 methods.
#' @rdname cuboid
#' @export
as.data.frame.cuboid = function(x, row.names = NULL, optional = FALSE, ...) {
  if (is.data.frame(x)) return(x)
  .dim = dim(x)
  start = attr(x, "start")
  vx = seq.int(start[[1L]], length.out = .dim[[1L]])
  vy = seq.int(start[[2L]], length.out = .dim[[2L]])
  vz = seq.int(start[[3L]], length.out = .dim[[3L]])
  res = expand_xyz(vx, vy, vz)
  res[["state"]] = c(x)
  attr(res, "dimarray") = .dim
  tibble::new_tibble(res, class = "cuboid")
}

# Convert data.frame to binary array
#' @rdname cuboid
#' @export
as.array.cuboid = function(x, ...) {
  if (is.array(x)) return(x)
  res = array(x[["state"]], dim = attr(x, "dimarray"))
  attr(res, "start") = x[c("x", "y", "z")] |> purrr::map_int(min)
  class(res) = c("cuboid", "array")
  res
}

expand_xyz = function(x, y = x, z = x) {
  # The order is consistent with array.
  # tidyr::expand_grid() sorts differently.
  expand.grid(x = x, y = y, z = z, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}
