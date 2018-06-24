#' R interface to tumopp, tumor growth simulator in C++
#' @aliases NULL tumopp-package
#' @useDynLib tumopp, .registration = TRUE
#' @importFrom magrittr %>%
#' @importFrom rlang .data
"_PACKAGE"

# to suppress NOTE
utils::globalVariables(c(".", "n"))

.onUnload = function(libpath) {
  message("Unloading tumopp")
  library.dynam.unload("tumopp", libpath)
}
