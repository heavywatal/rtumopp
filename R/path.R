#' Get information about the `tumopp` executable
#'
#' The simulator `tumopp` is a command-line tool written in C++.
#' The executable is built and installed along with this R package
#' unless the pre-installed one was found during the package installation.
#' @returns `tumopp_path()` returns the path to the `tumopp` executable.
#' @seealso [tumopp()] to run the simulator.
#' @rdname path
#' @export
#' @examples
#' \dontrun{
#' tumopp_path()
#'
#' tumopp_version()
#' }
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
  v_exe = tumopp("--version")
  v_pkg = utils::packageVersion("tumopp") |> as.character()
  if (extract_release(v_exe) != extract_release(v_pkg)) {
    warning(
      "C++ tumopp ",
      v_exe,
      " may be incompatible with R tumopp ",
      v_pkg,
      call. = FALSE
    )
  }
  v_exe
}

extract_release = function(x) {
  x |>
    stringr::str_remove("^v") |>
    stringr::str_extract("^\\d+\\.\\d+([-.]\\d+)?")
}
