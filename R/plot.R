#' ggplot for frequency spectrum
#' @param freqs numeric vector
#' @rdname plot
#' @export
histogram_freqspec = function(freqs) {
  tibble::tibble(x = freqs) |>
    ggplot2::ggplot(ggplot2::aes_(~x, ~..density..)) +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "frequency of alleles (or living descendants)")
}

#' ggplot for 2D lattice
#' @param .tbl tbl with extant cells
#' @param color column name to colorcode
#' @param size relative size of points
#' @param alpha opacity `[0, 1]`
#' @param limit for value range
#' @rdname plot
#' @export
plot_lattice2d = function(.tbl, color = "z", alpha = 1, size = 1, limit = max_abs_xyz(.tbl)) {
  size = size * 96 / (limit - 0.5)
  ggplot2::ggplot(.tbl, ggplot2::aes_(~x, ~y)) +
    ggplot2::geom_point(ggplot2::aes_string(color = color), alpha = alpha, size = size) +
    ggplot2::coord_equal(xlim = limit * c(-1, 1), ylim = limit * c(-1, 1))
}

# #######1#########2#########3#########4#########5#########6#########7#########

#' Plot capture_rate ~ nsam of biopsy
#' @param data tbl from [summarize_capture_rate()]` or [evaluate_mrs()]
#' @param point size of a data point
#' @param alpha opacity of a data point
#' @param errorbar logical
#' @rdname plot-biopsy
#' @export
plot_capture_rate = function(data, point = 1, alpha = 0.3, errorbar = TRUE) {
  # scales::viridis_pal()(5L) |> str_replace("FF$", "") |> tail(2L)
  p = ggplot2::ggplot(data, ggplot2::aes_(~nsam, ~capture_rate)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.0, fill = "#5DC863", alpha = 0.8) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 0.9, fill = "#FDE725", alpha = 0.8) +
    ggplot2::stat_summary(fun.y = mean, geom = "bar", fill = "#777777") +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 0.8, 0.9, 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1))
  if (point > 0 && alpha > 0) {
    p = p + ggplot2::geom_jitter(size = point, alpha = alpha, width = 0.23, height = 0)
  }
  if (errorbar) {
    p = p + ggplot2::stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2)
  }
  p + ggplot2::labs(x = "# of samples", y = "Capture rate")
}

mean_sd = function(x, mult = 1.96) {
  x = stats::na.omit(x)
  div = mult * stats::sd(x)
  mu = mean(x)
  data.frame(y = mu, ymin = mu - div, ymax = mu + div)
}

# #######1#########2#########3#########4#########5#########6#########7#########

#' Plot serial sections of 3D tumor
#' @param .tbl tbl with extant cells
#' @inheritParams ggplot2::ggsave
#' @param ... passed to [plot_lattice2d()]
#' @rdname plot-section
#' @export
save_serial_section = function(.tbl, filename = "png/section_%03d.png", scale = 6, dpi = 72, ...) {
  .lim = max_abs_xyz(.tbl)
  tidyr::nest(.tbl, -.data$z) |>
    dplyr::arrange(.data$z) |>
    dplyr::mutate(i = dplyr::row_number()) |>
    purrr::pwalk(function(z, data, i) {
      .outfile = sprintf(filename, i)
      # TODO: fix color for each lineage
      .p = plot_lattice2d(data, ..., limit = .lim) +
        ggplot2::geom_hline(yintercept = z[1L], color = "#999999", size = 1.5) +
        ggplot2::labs(title = sprintf("z =%4.1f", z)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        )
      ggplot2::ggsave(.outfile, .p, width = 1, height = 1, scale = scale, dpi = dpi)
    })
}
