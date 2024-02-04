#' @useDynLib tumopp, .registration = TRUE
#' @importFrom rlang :=
#' @importFrom rlang .data
#' @aliases NULL tumopp-package
#' @keywords internal
"_PACKAGE"

.onUnload = function(libpath) {
  message("Unloading tumopp in ", libpath)
  library.dynam.unload("tumopp", libpath)
}
