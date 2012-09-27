# set up project
project("dune-grid-howto" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)

# general stuff
cmake_minimum_required(VERSION 2.8)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("doc")
add_subdirectory("cmake/modules")

set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)

# TODO:
# set(examplesdir ${CMAKE_INSTALL_DOCDIR}/dune-grid-howto/examples)
#
# set(examples_HEADERS  elementdata.hh parfvdatahandle.hh unitcube_sgrid.hh
#         evolve.hh transportproblem2.hh
#         finitevolumeadapt.hh transportproblem.hh unitcube_yaspgrid.hh
#         functors.hh unitcube_albertagrid.hh vertexdata.hh
#         initialize.hh unitcube_alugrid.hh vtkout.hh integrateentity.hh
#         unitcube.hh parevolve.hh basicunitcube.hh
#         shapefunctions.hh)
#
# set(examples_PROGRAMS  gettingstarted traversal integration othergrids
#         adaptiveintegration finitevolume adaptivefinitevolume parfinitevolume
#         visualization finiteelements)
#
# set(examples_DATA  gettingstarted.cc traversal.cc integration.cc othergrids.cc
#         adaptiveintegration.cc finitevolume.cc adaptivefinitevolume.cc
#         parfinitevolume.cc visualization.cc finiteelements.cc)

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

set_property(TARGET adaptivefinitevolume adaptiveintegration finitevolume visualization APPEND
  PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" "WORLDDIM=${GRIDDIM}" "${GRIDTYPE}")

## TODO: take care for the following lines
##install(FILES "${examples_DATA}" DESTINATION ${examplesdir})
##install(TARGETS "${examples_PROGRAMS}" DESTINATION ${examplesdir})
##install(FILES "${examples_HEADERS}" DESTINATION ${examplesdir})

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
