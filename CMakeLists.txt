# set up project
project("dune-grid-howto" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

# guesse dune-common build directory
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("doc")
add_subdirectory("grids")

set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)

# The directory where the example grid files used for the examples are located
set_property(DIRECTORY APPEND PROPERTY
  COMPILE_DEFINITIONS "DUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/grids/\"")

set(examples_DATA  gettingstarted.cc traversal.cc integration.cc othergrids.cc
         adaptiveintegration.cc finitevolume.cc adaptivefinitevolume.cc
         parfinitevolume.cc visualization.cc finiteelements.cc)

set(TESTS
  adaptivefinitevolume
  adaptiveintegration
  finiteelements
  finitevolume
  gettingstarted
  integration
  othergrids
  parfinitevolume
  traversal
  visualization)

add_executable(adaptivefinitevolume adaptivefinitevolume.cc)
target_link_libraries(adaptivefinitevolume dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(adaptivefinitevolume)

add_executable(adaptiveintegration adaptiveintegration.cc)
target_link_libraries(adaptiveintegration dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(adaptiveintegration)

add_executable(finiteelements finiteelements.cc)
target_link_libraries(finiteelements dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(finiteelements)

add_executable(finitevolume finitevolume.cc)
target_link_libraries(finitevolume dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(finitevolume)

add_executable(integration integration.cc)
target_link_libraries(integration dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(integration)

add_executable(gettingstarted gettingstarted.cc)
target_link_libraries(gettingstarted dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(gettingstarted)

add_executable(othergrids othergrids.cc)
target_link_libraries(othergrids dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(othergrids)

add_executable(parfinitevolume parfinitevolume.cc)
target_link_libraries(parfinitevolume dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(parfinitevolume)

add_executable(traversal traversal.cc)
target_link_libraries(traversal dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(traversal)

add_executable(visualization visualization.cc)
target_link_libraries(visualization dunecommon dunegeometry dunegrid)
add_dune_mpi_flags(visualization)

set_property(TARGET
  adaptivefinitevolume
  adaptiveintegration
  finitevolume
  integration
  visualization APPEND
  PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" "WORLDDIM=${GRIDDIM}" "${GRIDTYPE}")

# We do not want want to build the tests during make all,
# but just build them on demand
SET(_test_target dune_grid_howto_tests)
add_custom_target(${_test_target})
dune_common_script_dir(SCRIPT_DIR)
configure_file(${SCRIPT_DIR}/BuildTests.cmake.in BuildTests.cmake @ONLY)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})

# install the header as done for the auto-tools
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
  unitcube_sgrid.hh
  unitcube_yaspgrid.hh
  vertexdata.hh
  vtkout.hh
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)

install(FILES ${examples_DATA} DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)
install(TARGETS ${TESTS} DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
