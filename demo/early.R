tumopp::tumopp("-h")
.const = list(D = 2L, k = 100, N = 256L, R = 256L)
.alt = list(
  C = c("moore", "hex"),
  L = c("const", "linear", "step")
)
args_table = tumopp::make_args(alt = .alt, const = .const, each = 1L) |> print()
results = tumopp::tumopp(args_table)
# tumopp::write_results(results)

.add_clade_to_snapshots = function(population, graph, snapshots, ...) {
  clade_info = tumopp:::sort_clades(population, graph, 4L)
  dplyr::left_join(snapshots, clade_info, by = "id")
}
# .tbl = (dplyr::slice(results, 1L) |> purrr::pmap(.add_clade_to_snapshots))[[1L]]

.plot_snapshot = function(data, limit) {
  okabe_ito = grDevices::palette.colors(palette = "Okabe-Ito")[-1]
  tumopp::plot_lattice2d(data, color = .data$clade, alpha = 1.0, limit = limit) +
    ggplot2::theme_void() +
    ggplot2::theme(palette.colour.discrete = okabe_ito, legend.position = "none")
}

.plot_snapshots = function(.tbl) {
  .lim = tumopp::max_abs_xyz(.tbl)
  tidyr::nest(.tbl, data = !"time")$data |>
    parallel::mclapply(.plot_snapshot, limit = .lim)
}

magick_gif_animation = function(
  infiles,
  outfile = "animation.gif",
  delay = 15,
  loop = 1
) {
  .args = c("-loop", loop, "-delay", delay, infiles)
  .args = c(.args, "-layers", "Optimize", outfile)
  message(paste(c("magick", .args), collapse = " "))
  system2("magick", .args)
}

.do = function(population, graph, snapshots, dest, ...) {
  .tbl = .add_clade_to_snapshots(population, graph, snapshots)
  .plt = .plot_snapshots(.tbl)
  .pngdir = file.path(dest, "png")
  dir.create(.pngdir, recursive = TRUE, mode = "0755")
  purrr::iwalk(.plt, \(.x, .y) {
    .outfile = file.path(.pngdir, sprintf("snapshot_%03d.png", .y))
    message(.outfile)
    ggplot2::ggsave(.outfile, .x, width = 1, height = 1, scale = 6, dpi = 72)
  })
  .infiles = file.path(.pngdir, "snapshot_*.png")
  magick_gif_animation(.infiles, sprintf("%s/%s.gif", dest, dest), delay = 8)
}

results |>
  dplyr::slice(3L) |>
  dplyr::mutate(dest = sprintf("C%s_L%s", .data$coord, .data$local)) |>
  purrr::pmap(.do)
