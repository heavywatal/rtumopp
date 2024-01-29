cmake_minimum_required(VERSION 3.15)
project(rtumopp
  VERSION 0.8.6
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
  get_filename_component(CMAKE_INSTALL_PREFIX ${INCLUDEDIR}/.. ABSOLUTE)
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
  set(BUILD_EXECUTABLE OFF CACHE BOOL "")
  set(BUILD_TESTING OFF CACHE BOOL "")
  FetchContent_MakeAvailable(tumopp)
  cmake_print_variables(tumopp_SOURCE_DIR tumopp_BINARY_DIR)
endif()

if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()
set(LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CPPFLAGS "-I${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
set(LDFLAGS "-L${LIBDIR} -Wl,-rpath,${LIBDIR}")
set(LDLIBS "-ltumopp")
cmake_print_variables(CPPFLAGS)
cmake_print_variables(LDFLAGS)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars.in
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars @ONLY
)
