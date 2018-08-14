#' Utilities for writing files
#'
#' @description
#' `write_results` writes nested data.frame to files
#' @param results data.frame returned from tumopp()
#' @rdname write
#' @export
write_results = function(results) {
  for (i in seq_len(nrow(results))) {
    dplyr::slice(results, i) %>% .write_result()
  }
}

.write_result = function(results, outdir = NULL, force = FALSE) {
  if (is.null(outdir)) outdir = results$outdir
  message("outdir: ", outdir)
  stopifnot(dir.create(outdir, mode = "0755") || force)
  cols = names(results)
  dfs = c("population", "snapshots", "drivers", "distances")
  purrr::walk(dfs, ~{
    if (.x %in% cols) {
      content = results[[.x]][[1L]]
      if (.x %in% c("population", "snapshots") &&
          results$coord == "hex" &&
          !is.integer(content$x)) {
        content = revert_coord_hex(content)
      }
      outfile = file.path(outdir, paste0(.x, ".tsv.gz"))
      readr::write_tsv(content, outfile)
    }
  })
  conf = suppressWarnings(dplyr::select(results, -dplyr::one_of(c(dfs, "graph"))))
  readr::write_tsv(conf, file.path(outdir, "program_options.tsv.gz"))
  invisible(results)
}
