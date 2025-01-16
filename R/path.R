#' Get the path to the tumopp installation.
#'
#' @rdname path
#' @export
tumopp_path = function() {
  x = tumopp_path_config()
  if (!file.exists(x)) {
    x = tumopp_path_exec()
    if (!file.exists(x)) x = "tumopp"
  }
  x
}

tumopp_path_exec = function() {
  system.file("exec", "tumopp", package = "tumopp")
}

#' @rdname path
#' @export
tumopp_version = function() {
  system2(tumopp_path(), "--version", stdout = TRUE, stderr = FALSE)
}
