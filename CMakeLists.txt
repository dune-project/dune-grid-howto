# set up project
project(dune-grid-howto C CXX)

# general stuff
cmake_minimum_required(VERSION 3.1)

message("\n\
        DUNE-grid-howto is deprecated. Changes in the DUNE core modules are\n\
        no longer reflected in this module, so it might break any time.\n\
        Release 2.7 was the last release, there will be no 2.8 release.\n")

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

option(DUNE_GRID_GRIDTYPE_SELECTOR "" ON)

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# The directory where the example grid files used for the examples are located
set_property(DIRECTORY APPEND PROPERTY
  COMPILE_DEFINITIONS "DUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/grids/\"")

# force that all tests are built
set(DUNE_BUILD_TESTS_ON_MAKE_ALL TRUE)

# enable the flags for all third party features found
dune_enable_all_packages()

dune_add_test(SOURCES adaptivefinitevolume.cc
  COMPILE_DEFINITIONS "GRIDDIM=2" "WORLDDIM=2" "YASPGRID")
add_dune_alberta_flags(adaptivefinitevolume WORLDDIM 2)

dune_add_test(SOURCES adaptiveintegration.cc
  COMPILE_DEFINITIONS "GRIDDIM=2" "WORLDDIM=2" "YASPGRID")

dune_add_test(SOURCES finiteelements.cc)
add_dune_alberta_flags(finiteelements WORLDDIM 2)

dune_add_test(SOURCES finitevolume.cc
  COMPILE_DEFINITIONS "GRIDDIM=2" "WORLDDIM=2" "YASPGRID")

dune_add_test(SOURCES integration.cc
  COMPILE_DEFINITIONS "GRIDDIM=2" "WORLDDIM=2" "YASPGRID")

dune_add_test(SOURCES gettingstarted.cc)

dune_add_test(SOURCES othergrids.cc)
add_dune_ug_flags(othergrids)

dune_add_test(SOURCES parfinitevolume.cc)

dune_add_test(SOURCES traversal.cc)

dune_add_test(SOURCES visualization.cc
  COMPILE_DEFINITIONS "GRIDDIM=2" "WORLDDIM=2" "YASPGRID")

add_subdirectory(doc)
add_subdirectory(grids)

# install headers, cc files and executables
install(FILES
  basicunitcube.hh
  elementdata.hh
  evolve.hh
  finitevolumeadapt.hh transportproblem.hh
  functors.hh unitcube_albertagrid.hh
  initialize.hh
  integrateentity.hh
  parfvdatahandle.hh
  parevolve.hh
  shapefunctions.hh
  transportproblem2.hh
  unitcube.hh
  unitcube_alugrid.hh
  unitcube_yaspgrid.hh
  vertexdata.hh
  vtkout.hh
  adaptivefinitevolume.cc
  adaptiveintegration.cc
  gettingstarted.cc
  integration.cc
  othergrids.cc
  finiteelements.cc
  finitevolume.cc
  parfinitevolume.cc
  traversal.cc
  visualization.cc
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)

install(TARGETS
  adaptivefinitevolume
  adaptiveintegration
  gettingstarted
  integration
  othergrids
  finiteelements
  finitevolume
  parfinitevolume
  traversal
  visualization
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
