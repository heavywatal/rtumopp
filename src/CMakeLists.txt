cmake_minimum_required(VERSION 3.1)
project(tumoppr CXX)

find_package(tumopp 0.6.0 REQUIRED CONFIG)

include(GNUInstallDirs)
get_target_property(CXX_INCLUDEDIR tumopp::tumopp INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "INCLUDEDIR: ${CXX_INCLUDEDIR}")
if(${CXX_INCLUDEDIR} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()
get_filename_component(CXX_LIBDIR ${CXX_INCLUDEDIR}/../${CMAKE_INSTALL_LIBDIR} ABSOLUTE)
message(STATUS "LIBDIR: ${CXX_LIBDIR}")

message(STATUS "Configure Makevars")
set(CPPFLAGS "-I${CXX_INCLUDEDIR}")
set(LDFLAGS "-L${CXX_LIBDIR} -Wl,-rpath,${CXX_LIBDIR}")
set(LDLIBS "-ltumopp")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars.in
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars @ONLY
)
