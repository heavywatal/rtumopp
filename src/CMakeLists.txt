cmake_minimum_required(VERSION 3.1)
project(rtumopp
  VERSION 0.6.2
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
  find_package(Git)
  set(package tumopp)
  set(repo heavywatal/${package})
  if(${PROJECT_VERSION} VERSION_LESS 0.1.0)
    set(version HEAD)
    set(branch master)
  else()
    set(version ${PROJECT_VERSION})
    set(branch v${version})
  endif()
  set(subdir external/${package}-${version})
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/.git)
    execute_process(COMMAND
      ${GIT_EXECUTABLE} clone --recursive --depth=1 --branch=${branch} https://github.com/${repo}.git ${subdir}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  endif()
  message(STATUS "add_subdirectory(${subdir})")
  set(CMAKE_BUILD_TYPE Release)
  add_subdirectory(${subdir})
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
