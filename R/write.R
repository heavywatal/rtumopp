#' Utilities for writing files
#'
#' @details
#' `write_results` writes nested data.frame to files
#' @param results data.frame returned from [tumopp()]
#' @rdname write
#' @export
write_results = function(results) {
  dplyr::rowwise(results) |> dplyr::group_walk(\(.x, .y) .write_result(.x))
}

.write_result = function(result, outdir = NULL, force = FALSE) {
  if (is.null(outdir)) {
    outdir = result$outdir
  }
  message("outdir: ", outdir)
  stopifnot(dir.create(outdir, mode = "0755") || force)
  cols = names(result)
  dfs = c("population", "snapshots", "drivers", "distances")
  purrr::walk(dfs, \(.x) {
    if (.x %in% cols) {
      content = result[[.x]][[1L]]
      if (
        .x %in%
          c("population", "snapshots") &&
          result$coord == "hex" &&
          !is.integer(content$x)
      ) {
        content = revert_coord_hex(content)
      }
      outfile = file.path(outdir, paste0(.x, ".tsv.gz"))
      readr::write_tsv(content, outfile)
    }
  })
  conf = suppressWarnings(dplyr::select(result, !dplyr::one_of(c(dfs, "graph"))))
  json = to_json(as.list(conf))
  cat(json, file = file.path(outdir, "config.json"))
  invisible(result)
}

to_json = function(x, ...) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    digits = I(15),
    pretty = 2L,
    always_decimal = TRUE,
    ...
  )
}
