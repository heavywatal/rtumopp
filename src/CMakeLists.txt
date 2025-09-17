cmake_minimum_required(VERSION 3.28)
project($ENV{R_PACKAGE_NAME}
  VERSION 0.9.2
  LANGUAGES CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)
cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

set(git_tag "$ENV{TUMOPP_GIT_TAG}")
if(NOT git_tag)
  set(git_tag 2b1c8195e8de54b41171fcd6c25f30149a3ba886)
  set(find_package_args FIND_PACKAGE_ARGS ${PROJECT_VERSION})
endif()
FetchContent_Declare(tumopp
  GIT_REPOSITORY https://github.com/heavywatal/tumopp.git
  GIT_TAG ${git_tag}
  ${find_package_args}
)

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
