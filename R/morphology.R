#' Extract surface cells with mathematical morphology
#'
#' @description
#' `detect_surface` adds a binary column.
#' @param .tbl a data.frame with (x, y, z) columns
#' @param se structuring element as a binary array
#' @return tibble with surface column
#' @rdname morphology
#' @export
detect_surface = function(.tbl, se) {
  if (nrow(.tbl) == 0L) {
    return(dplyr::mutate(.tbl, surface = logical(0L)))
  }
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

#' @description
#' `get_se` construct a structuring element (kernel) as a binary array.
#' @param coord a string
#' @param dimensions an integer
#' @rdname morphology
#' @export
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
  df2img(df)
}

# Filter img by surface
filter_surface = function(img, se) {
  img - mmand::erode(img, se)
}

# Convert data.frame to binary array
df2img = function(.tbl) {
  vars = c("x", "y", "z")
  .tbl = .tbl[vars]
  .grid = dplyr::summarise_at(.tbl, vars, dplyr::funs(min, max)) %>% {
    expand.grid(
      x = seq(.$x_min, .$x_max),
      y = seq(.$y_min, .$y_max),
      z = seq(.$z_min, .$z_max)
    )
  } %>% tibble::as_tibble()
  joined = dplyr::left_join(.grid, dplyr::mutate(.tbl, v = 1L), by = vars)
  joined = tidyr::replace_na(joined, list(v = 0L))
  arr = reshape2::acast(joined, x ~ y ~ z, `[`, 1L, value.var = "v", fill = 0L)
  dim(arr) = c(dim(arr), 1L)
  arr
}

# Convert binary array to data.frame with (x, y, z) columns
img2df = function(img) {
  img %>%
    {
      dim(.) = utils::head(dim(.), 3L)
      .
    } %>%
    reshape2::melt(c("x", "y", "z")) %>%
    tibble::as_tibble()
}
