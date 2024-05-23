cmake_minimum_required(VERSION 3.18)
project(rtumopp
  VERSION 0.8.7.1
  LANGUAGES CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

find_package(tumopp ${PROJECT_VERSION} QUIET)
cmake_print_variables(tumopp_FOUND tumopp_DIR)
if(tumopp_FOUND)
  get_target_property(INCLUDEDIR tumopp::tumopp INTERFACE_INCLUDE_DIRECTORIES)
  get_filename_component(CMAKE_INSTALL_PREFIX "${INCLUDEDIR}" DIRECTORY)
  cmake_print_variables(CMAKE_INSTALL_PREFIX)
  install(CODE "message(STATUS ${tumopp_DIR})") # to make target 'install'
else()
  if(DEFINED ENV{FETCHCONTENT_SOURCE_DIR_TUMOPP})
    set(FETCHCONTENT_SOURCE_DIR_TUMOPP $ENV{FETCHCONTENT_SOURCE_DIR_TUMOPP})
    set(tag "HEAD")
  else()
    set(tag "v${PROJECT_VERSION}")
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
  set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
  set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
  set(CMAKE_INSTALL_BINDIR "exec")
  set(CMAKE_INSTALL_LIBDIR "libs")
  set(BUILD_EXECUTABLE ON CACHE BOOL "")
  set(BUILD_TESTING OFF CACHE BOOL "")
  FetchContent_MakeAvailable(tumopp)
  cmake_print_variables(tumopp_SOURCE_DIR tumopp_BINARY_DIR)
endif()

set(BINDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
cmake_print_variables(BINDIR)
get_filename_component(R_PACKAGE_SOURCE ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
file(WRITE ${R_PACKAGE_SOURCE}/inst/path "${BINDIR}/tumopp")
