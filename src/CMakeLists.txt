cmake_minimum_required(VERSION 3.15)
project(rtumopp
  VERSION 0.8.5
  LANGUAGES CXX)

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(tumopp ${PROJECT_VERSION} QUIET)
if(tumopp_FOUND)
  get_target_property(INCLUDEDIR tumopp::tumopp INTERFACE_INCLUDE_DIRECTORIES)
  get_filename_component(CMAKE_INSTALL_PREFIX ${INCLUDEDIR}/.. ABSOLUTE)
  install(CODE "message(STATUS ${tumopp_DIR})") # to make target 'install'
else()
  message(STATUS "FETCHCONTENT_SOURCE_DIR_TUMOPP: ${FETCHCONTENT_SOURCE_DIR_TUMOPP}")
  if(IS_DIRECTORY ${FETCHCONTENT_SOURCE_DIR_TUMOPP})
    set(tag "HEAD")
  else()
    set(tag "v${PROJECT_VERSION}")
  endif()
  message(STATUS "tag: ${tag}")
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)
  FetchContent_Declare(
    tumopp
    GIT_REPOSITORY https://github.com/heavywatal/tumopp.git
    GIT_TAG ${tag}
    GIT_SHALLOW ON
  )
  set(CMAKE_BUILD_TYPE Release)
  set(BUILD_EXECUTABLE OFF CACHE BOOL "")
  set(BUILD_TESTING OFF CACHE BOOL "")
  FetchContent_MakeAvailable(tumopp)
  message(STATUS "tumopp_SOURCE_DIR: ${tumopp_SOURCE_DIR}")
endif()

if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()
set(LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CPPFLAGS "-I${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
set(LDFLAGS "-L${LIBDIR} -Wl,-rpath,${LIBDIR}")
set(LDLIBS "-ltumopp")
message(STATUS "CPPFLAGS: ${CPPFLAGS}")
message(STATUS "LDFLAGS: ${LDFLAGS}")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars.in
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars @ONLY
)
