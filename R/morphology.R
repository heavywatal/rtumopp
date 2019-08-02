#' Extract surface cells with mathematical morphology
#'
#' @details
#' `add_surface` adds a binary column to the data.frame
#' @param population a data.frame with (x, y, z)
#' @param coord a string
#' @param dimensions an integer
#' @rdname morphology
#' @export
add_surface = function(population, coord, dimensions) {
  extant = filter_extant(population)
  strelem = get_se(coord, dimensions) %>% df2img()
  col_surface = detect_surface(extant, strelem) %>%
    dplyr::select(.data$id, .data$surface)
  population %>%
    dplyr::left_join(col_surface, by = "id")
}

#' @details
#' `add_phi` counts empty neighbors of each cell and adds an integer column
#' @rdname morphology
#' @export
add_phi = function(population, coord, dimensions) {
  extant = filter_extant(population)
  coords = extant[c("x", "y", "z", "id")]
  kernel = get_se(coord, dimensions)
  counts = purrr::pmap_dfr(kernel, ~ {
    dplyr::mutate(coords, x = .data$x + ..1, y = .data$y + ..2, z = .data$z + ..3)
  }) %>%
    dplyr::count(.data$x, .data$y, .data$z)
  info = coords %>%
    dplyr::left_join(counts, by = c("x", "y", "z")) %>%
    dplyr::transmute(.data$id, phi = nrow(kernel) - .data$n)
  population %>%
    dplyr::left_join(info, by = c("id"))
}

detect_surface = function(.tbl, se) {
  if (nrow(.tbl) == 0L) {
    return(dplyr::mutate(.tbl, surface = logical(0L)))
  }
  if (!is.integer(.tbl$x)) .tbl = revert_coord_hex(.tbl)
  axes = c("x", "y", "z")
  mins = dplyr::summarise_at(.tbl, axes, dplyr::funs(min))
  img = df2img(.tbl) %>% filter_surface(se)
  product = img2df(img) %>% dplyr::transmute(
    x = .data$x + mins$x - 1L,
    y = .data$y + mins$y - 1L,
    z = .data$z + mins$z - 1L,
    surface = .data$value > 0L
  )
  dplyr::left_join(.tbl, product, by = axes)
}

# Construct a structuring element (kernel)
get_se = function(coord = c("moore", "neumann", "hex"), dimensions = 3L) {
  coord = match.arg(coord)
  v = c(-1L, 0L, 1L)
  df = tidyr::crossing(x = v, y = v, z = v)
  if (coord == "neumann") {
    df = dplyr::filter(df, abs(.data$x) + abs(.data$y) + abs(.data$z) < 2)
  } else if (coord == "hex") {
    df = dplyr::filter(df, (trans_coord_hex(df) %>% dist_euclidean()) < 1.1)
  }
  if (dimensions < 3L) {
    df = dplyr::filter(df, .data$z == 0L)
  }
  df
}

# Filter img by surface
filter_surface = function(img, se) {
  img - mmand::erode(img, se)
}

# Convert data.frame to binary array
df2img = function(.tbl) {
  vars = c("x", "y", "z")
  .tbl = .tbl[vars]
  .summary = dplyr::summarise_at(.tbl, vars, dplyr::funs(min, max))
  .grid = tidyr::crossing(
    x = seq(.summary$x_min, .summary$x_max),
    y = seq(.summary$y_min, .summary$y_max),
    z = seq(.summary$z_min, .summary$z_max)
  )
  joined = dplyr::left_join(.grid, dplyr::mutate(.tbl, v = 1L), by = vars)
  joined = tidyr::replace_na(joined, list(v = 0L))
  arr = reshape2::acast(joined, x ~ y ~ z, `[`, 1L, value.var = "v", fill = 0L)
  dim(arr) = c(dim(arr), 1L)
  arr
}

# Convert binary array to data.frame with (x, y, z) columns
img2df = function(img) {
  dim(img) = utils::head(dim(img), 3L)
  img %>%
    reshape2::melt(c("x", "y", "z")) %>%
    tibble::as_tibble()
}
