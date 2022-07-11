library(tidyverse)
library(tumopp)

.args = list(
  "-D2 -k100 -Cmoore -Lconst -R256 -N256 -o Cmoore_Lconst",
  "-D2 -k100 -Cmoore -Lstep -R256 -N256 -o Cmoore_Lstep",
  "-D2 -k100 -Cmoore -Llinear -R256 -N256 -o Cmoore_Llinear",
  "-D2 -k100 -Chex -Lconst -R256 -N256 -o Chex_Lconst",
  "-D2 -k100 -Chex -Lstep -R256 -N256 -o Chex_Lstep",
  "-D2 -k100 -Chex -Llinear -R256 -N256 -o Chex_Llinear"
)
results = tumopp(.args)
write_results(results)

.add_clade_to_snapshots = function(population, graph, snapshots, ...) {
  clade_info = tumopp:::sort_clades(population, graph, 4L)
  dplyr::left_join(snapshots, clade_info, by = "id")
}
# .tbl = (dplyr::slice(results, 1L) %>% purrr::pmap(.add_clade_to_snapshots))[[1L]]

.plot_snapshot = function(data, limit) {
  tumopp::plot_lattice2d(data, "clade", alpha = 1.0, limit = limit) +
    scale_colour_brewer(palette = "Spectral", na.value = "grey50", guide = FALSE) +
    theme_void()
}

.plot_snapshots = function(.tbl) {
  .lim = tumopp::max_abs_xyz(.tbl)
  tidyr::nest(.tbl, -time)$data %>%
    parallel::mclapply(.plot_snapshot, limit = .lim)
}

magick_gif_animation = function(infiles, outfile = "animation.gif", delay = 15, loop = 1) {
  args = c(
    "-loop", loop, "-delay", delay,
    infiles, "-layers", "Optimize", outfile
  )
  message(paste(args, collapse = " "))
  system2("magick", args)
}

.do = function(population, graph, snapshots, outdir, ...) {
  .tbl = .add_clade_to_snapshots(population, graph, snapshots)
  .plt = .plot_snapshots(.tbl)
  .pngdir = file.path(outdir, "png")
  dir.create(.pngdir, recursive = TRUE, mode = "0755")
  purrr::iwalk(.plt, ~ {
    .outfile = file.path(.pngdir, sprintf("snapshot_%03d.png", .y))
    message(.outfile)
    ggsave(.outfile, .x, width = 1, height = 1, scale = 6, dpi = 72)
  })
  .infiles = file.path(.pngdir, "snapshot_*.png")
  magick_gif_animation(.infiles, sprintf("%s/%s.gif", outdir, outdir), delay = 8)
}
dplyr::slice(results, 1L) %>% purrr::pmap(.do)
