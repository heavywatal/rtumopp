#' R interface to tumopp, tumor growth simulator in C++
#' @aliases NULL tumopp-package
#' @useDynLib tumopp
#' @importFrom magrittr %>%
#' @importFrom rlang .data
"_PACKAGE"

# to suppress NOTE
utils::globalVariables(c(".", "n"))
