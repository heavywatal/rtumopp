# R interface to tumopp

[![R-CMD-check](https://github.com/heavywatal/rtumopp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heavywatal/rtumopp/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/heavywatal/rtumopp/branch/master/graph/badge.svg)](https://app.codecov.io/gh/heavywatal/rtumopp?branch=master)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/rtumopp)](https://cran.r-project.org/package=rtumopp)

This is an R interface to [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Requirements

- Unix-like OS (macOS, Linux, etc.)
- C++17 compiler (clang++ >= Apple LLVM 12, g++ >= 8)
- [CMake](https://cmake.org/) (>= 3.18.0)

## Installation

```r
install.packages("devtools")
devtools::install_github("heavywatal/rtumopp")
```

R packages are updated at random times.
Please try to check updates once in a while and use the latest versions.

```r
update.packages()
devtools::install_github("heavywatal/rigraphlite")
devtools::install_github("heavywatal/rtumopp")
```

## Basic usage

See ["Get started" page](https://heavywatal.github.io/rtumopp/articles/tumopp.html).

Available parameters are listed in
[the API document of C++ tumopp](https://heavywatal.github.io/tumopp/group__params.html).
