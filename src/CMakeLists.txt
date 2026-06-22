cmake_minimum_required(VERSION 3.28)
project($ENV{R_PACKAGE_NAME}
  VERSION 0.9.3
  LANGUAGES CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

if(DEFINED ENV{TUMOPP_GIT_TAG})
  FetchContent_Declare(tumopp
    GIT_REPOSITORY https://github.com/heavywatal/tumopp.git
    GIT_TAG "$ENV{TUMOPP_GIT_TAG}"
  )
else()
  FetchContent_Declare(tumopp
    URL https://github.com/heavywatal/tumopp/archive/refs/tags/v0.9.3.tar.gz
    URL_HASH SHA256=c87b9d2c440b0c5c951016650b59b051261007f13b05715219ae3cc63da7991c
    FIND_PACKAGE_ARGS ${PROJECT_VERSION}
  )
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_INSTALL_BINDIR "exec")
FetchContent_MakeAvailable(tumopp)

cmake_print_variables(tumopp_FOUND tumopp_DIR)
if(tumopp_FOUND)
  get_property(_cfg TARGET tumopp::tumopp-exe PROPERTY IMPORTED_CONFIGURATIONS)
  get_property(exec_path TARGET tumopp::tumopp-exe PROPERTY IMPORTED_LOCATION_${_cfg})
else()
  cmake_print_variables(tumopp_SOURCE_DIR tumopp_BINARY_DIR)
  set(exec_path "$ENV{R_LIBRARY_DIR}/${PROJECT_NAME}/exec/tumopp")
endif()

cmake_print_variables(exec_path)
configure_file("config.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/config.cpp" @ONLY)

cmake_print_variables(CMAKE_CURRENT_BINARY_DIR)
if("${CMAKE_CURRENT_BINARY_DIR}" MATCHES "/build$")
  add_library(compile_commands OBJECT config.cpp cpp11.cpp)
  target_include_directories(compile_commands PRIVATE "$ENV{R_INCLUDE_DIR}")
  string(REPLACE ":" ";" R_LIBS_USER "$ENV{R_LIBS_USER}")
  string(REPLACE ":" ";" R_LIBS_SITE "$ENV{R_LIBS_SITE}")
  set(R_HOME "$ENV{R_HOME}/library")
  foreach(r_lib IN LISTS R_LIBS_USER R_LIBS_SITE R_HOME)
    if(IS_DIRECTORY "${r_lib}/cpp11/include")
      target_include_directories(compile_commands PRIVATE "${r_lib}/cpp11/include")
    endif()
  endforeach()
endif()
