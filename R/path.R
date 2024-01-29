#' Get the path to the tumopp installation.
#'
#' @rdname path
#' @export
tumopp_path = function() {
  scan(system.file("path", package = "tumopp"), what = character(0), quiet = TRUE)
}

#' @rdname path
#' @export
tumopp_version = function() {
  system2(tumopp_path(), "--version", stdout=TRUE, stderr=FALSE)
}
