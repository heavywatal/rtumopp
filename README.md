# R interface to tumopp

[![Build Status](https://travis-ci.org/heavywatal/rtumopp.svg?branch=master)](https://travis-ci.org/heavywatal/rtumopp)

This is an R interface to [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Requirements

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.1.0)

## Installation

```r
# install.packages("devtools")
devtools::install_github("heavywatal/rtumopp")
```

## Basic usage

See ["Get started" page](http://heavywatal.github.io/rtumopp/articles/tumopp.html).

Available parameters are listed in
[the API document of C++ tumopp](https://heavywatal.github.io/tumopp/group__params.html).
