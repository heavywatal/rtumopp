library(tidyverse)
library(wtl)
library(tumopp)
refresh('rtumopp')

# tumopp -D2 -k100 -Cmoore -Lconst -O4 -R256 -N256 0 0 -o Cmoore_Lconst
# tumopp -D2 -k100 -Cmoore -Lstep -O4 -R256 -N256 0 0 -o Cmoore_Lstep
# tumopp -D2 -k100 -Cmoore -Llinear -O4 -R256 -N256 0 0 -o Cmoore_Llinear
# tumopp -D2 -k100 -Chex -Lconst -O4 -R256 -N256 0 0 -o Chex_Lconst
# tumopp -D2 -k100 -Chex -Lstep -O4 -R256 -N256 0 0 -o Chex_Lstep
# tumopp -D2 -k100 -Chex -Llinear -O4 -R256 -N256 0 0 -o Chex_Llinear

(.args = wtl::command_args())
indir = .args$args[1]
if (!is.na(indir)) {
  setwd(indir)
}
setwd("Cmoore_Lconst")
results = tumopp::read_snapshots() %>% print()
.tbl = results$snapshots[[1L]]
.lim = tumopp::max_abs_xyz(.tbl)

plot_snapshot = function(data, time) {
  .N = nrow(data)
  tumopp::plot_lattice2d(data, "clade", alpha = 1.0, limit = .lim) +
    scale_color_brewer(palette = "Spectral", guide = FALSE) +
    # labs(title=sprintf('t = %.5f, N =%4d', time, .N))+
    theme_bw()+
    wtl::erase(axis.title, axis.text, axis.ticks)
}

.out = .tbl %>%
  tidyr::nest(-time) %>%
  dplyr::mutate(plt = purrr::pmap(., plot_snapshot)) %>%
  print()

dir.create("png", mode = "0755")
purrr::iwalk(.out$plt, ~{
  .outfile = sprintf("png/snapshot_%d.png", .y)
  message(.outfile)
  ggsave(.outfile, .x, width = 1, height = 1, scale = 6, dpi = 72)
})

# magick -loop 1 -delay 8 png/*.png snapshot-noloop.gif
# magick -loop 1 snapshot-noloop.gif -layers Optimize snapshot-noloop-opt.gif

if (FALSE) {
  # tiled-png for CSS sprite
  grob = cowplot::plot_grid(plotlist = .out$plt, nrow = 1)
  ggsave("earlysteps.png", grob, width = 54, height = 1, scale = 3, limitsize = FALSE)
}
