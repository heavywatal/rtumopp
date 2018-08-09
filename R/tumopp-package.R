#' @useDynLib tumopp, .registration = TRUE
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @aliases NULL tumopp-package
#' @keywords internal
"_PACKAGE"

# to suppress NOTE
utils::globalVariables(c(".", "n"))

.onLoad = function(libname, pkgname) {
  igraph::igraph_options(return.vs.es = FALSE)
}

.onUnload = function(libpath) {
  message("Unloading tumopp")
  library.dynam.unload("tumopp", libpath)
}
