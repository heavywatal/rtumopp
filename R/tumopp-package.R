#' @useDynLib tumopp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom rlang .data
#' @aliases NULL tumopp-package
#' @keywords internal
"_PACKAGE"

.onLoad = function(libname, pkgname) {
  igraph::igraph_options(
    add.vertex.names = FALSE,
    return.vs.es = FALSE
  )
}

.onUnload = function(libpath) {
  message("Unloading tumopp in ", libpath)
  library.dynam.unload("tumopp", libpath)
}
