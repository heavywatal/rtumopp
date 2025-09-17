# R interface to tumopp

[![R-CMD-check](https://github.com/heavywatal/rtumopp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heavywatal/rtumopp/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/heavywatal/rtumopp/graph/badge.svg?token=S8QayYTeJ0)](https://codecov.io/gh/heavywatal/rtumopp)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/rtumopp)](https://cran.r-project.org/package=rtumopp)

This is an R interface to [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Requirements

- Unix-like OS (macOS, Linux, etc.)
- C++17 compiler (clang++ >= Apple clang 11.0, g++ >= 9.1)
- [CMake](https://cmake.org/)

## Installation

```r
install.packages("pak")
pak::pak("heavywatal/rtumopp")
```

R packages are updated at random times.
Please try to check updates once in a while and use the latest versions.

```r
update.packages()
pak::pak("heavywatal/rigraphlite")
pak::pak("heavywatal/rtumopp")
```

## Basic usage

See ["Get started" page](https://heavywatal.github.io/rtumopp/articles/tumopp.html).

Available parameters are listed in
[the API document of C++ tumopp](https://heavywatal.github.io/tumopp/group__params.html).
