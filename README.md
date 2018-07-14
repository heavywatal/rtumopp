# R interface to tumopp

This is an R interface to [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Installation

1.  Build and install [tumopp](https://github.com/heavywatal/tumopp) with Homebrew/Linuxbrew or CMake.

1.  Install [devtools](https://github.com/hadley/devtools) in R:
    `install.packages('devtools')`

1.  Execute `devtools::install_github('heavywatal/rtumopp')` in R.
    You may need `Sys.setenv(CMAKE_PREFIX_PATH='/prefix/to/tumopp')` to tell R the location of tumopp installation.
