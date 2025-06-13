cmake_minimum_required(VERSION 3.21)
project($ENV{R_PACKAGE_NAME}
  VERSION 0.9.0
  LANGUAGES CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

set(tag "$ENV{TUMOPP_GIT_TAG}")
if(NOT tag)
  set(tag "v${PROJECT_VERSION}")
  find_package(tumopp ${PROJECT_VERSION} QUIET)
  cmake_print_variables(tumopp_FOUND tumopp_DIR)
endif()
if(tumopp_FOUND)
  set(exec_INSTALL_PREFIX "${tumopp_DIR}/../../..")
else()
  set(exec_INSTALL_PREFIX "$ENV{R_LIBRARY_DIR}/${PROJECT_NAME}") # load_all()
  set(FETCHCONTENT_SOURCE_DIR_TUMOPP "$ENV{FETCHCONTENT_SOURCE_DIR_TUMOPP}")
  if(FETCHCONTENT_SOURCE_DIR_TUMOPP AND NOT DEFINED ENV{TUMOPP_GIT_TAG})
    set(tag "HEAD")
  endif()
  cmake_print_variables(FETCHCONTENT_SOURCE_DIR_TUMOPP tag)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(
    tumopp
    GIT_REPOSITORY https://github.com/heavywatal/tumopp.git
    GIT_TAG ${tag}
    GIT_SHALLOW ON
    # FIND_PACKAGE_ARGS  # 3.24
  )
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_INSTALL_BINDIR "exec")
  set(BUILD_TESTING OFF CACHE BOOL "")
  FetchContent_MakeAvailable(tumopp)
  cmake_print_variables(tumopp_SOURCE_DIR tumopp_BINARY_DIR)
endif()

set(R_PACKAGE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(exec_path "${exec_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/tumopp")
cmake_print_variables(exec_path)
configure_file("config.cpp.in" "${R_PACKAGE_SOURCE}/src/config.cpp" @ONLY)
