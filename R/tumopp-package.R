#' @useDynLib tumopp, .registration = TRUE
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @aliases NULL tumopp-package
#' @keywords internal
"_PACKAGE"

# to suppress NOTE
utils::globalVariables(c(".", "n"))

.onUnload = function(libpath) {
  message("Unloading tumopp")
  library.dynam.unload("tumopp", libpath)
}
