# $Id$

SUBDIRS = doc grids

examplesdir=$(datadir)/doc/dune-grid-howto/examples

examples_HEADERS = elementdata.hh parfvdatahandle.hh unitcube_sgrid.hh \
	evolve.hh transportproblem2.hh unitcube_uggrid.hh \
	finitevolumeadapt.hh transportproblem.hh unitcube_yaspgrid.hh \
	functors.hh unitcube_albertagrid.hh vertexdata.hh \
	initialize.hh unitcube_alugrid.hh vtkout.hh integrateentity.hh \
	unitcube.hh parevolve.hh unitcube_onedgrid.hh

examples_PROGRAMS = gettingstarted traversal integration othergrids\
	adaptiveintegration finitevolume adaptivefinitevolume parfinitevolume\
	visualization

examples_DATA = gettingstarted.cc traversal.cc integration.cc othergrids.cc \
	adaptiveintegration.cc finitevolume.cc adaptivefinitevolume.cc \
	parfinitevolume.cc visualization.cc

gettingstarted_SOURCES = gettingstarted.cc
gettingstarted_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
gettingstarted_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

traversal_SOURCES = traversal.cc
traversal_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
traversal_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

integration_SOURCES = integration.cc
integration_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
integration_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

othergrids_SOURCES = othergrids.cc
othergrids_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
othergrids_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

adaptiveintegration_SOURCES = adaptiveintegration.cc
adaptiveintegration_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
adaptiveintegration_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

finitevolume_SOURCES = finitevolume.cc
finitevolume_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
finitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

adaptivefinitevolume_SOURCES = adaptivefinitevolume.cc
adaptivefinitevolume_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
adaptivefinitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

parfinitevolume_SOURCES = parfinitevolume.cc
parfinitevolume_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
parfinitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

visualization_SOURCES = visualization.cc
visualization_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
visualization_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid="$(DUNE_GRID_ROOT)" --with-dune-istl="$(DUNE_ISTL_ROOT)" --with-dune-disc="$(DUNE_DISC_ROOT)" CXX="$(CXX)" CC="$(CC)" --enable-parallel=@ENABLE_PARALLEL@ MPICC="$(MPICC)" --disable-documentation

include $(top_srcdir)/am/global-rules
