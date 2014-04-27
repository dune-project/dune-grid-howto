SUBDIRS = doc grids

GRIDTYPE=YASPGRID
GRIDDIM=2

examplesdir=$(datadir)/doc/dune-grid-howto/examples

examples_HEADERS = elementdata.hh parfvdatahandle.hh unitcube_sgrid.hh \
	evolve.hh transportproblem2.hh \
	finitevolumeadapt.hh transportproblem.hh unitcube_yaspgrid.hh \
	functors.hh unitcube_albertagrid.hh vertexdata.hh \
	initialize.hh unitcube_alugrid.hh vtkout.hh integrateentity.hh \
	unitcube.hh parevolve.hh basicunitcube.hh \
	shapefunctions.hh

examples_PROGRAMS = gettingstarted traversal integration othergrids\
	adaptiveintegration finitevolume adaptivefinitevolume parfinitevolume\
	visualization finiteelements

examples_DATA = gettingstarted.cc traversal.cc integration.cc othergrids.cc \
	adaptiveintegration.cc finitevolume.cc adaptivefinitevolume.cc \
	parfinitevolume.cc visualization.cc finiteelements.cc

# The directory where the example grid files used for the examples are located
AM_CPPFLAGS+=-DDUNE_GRID_HOWTO_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/grids/\"

#
## define the programs
#

gettingstarted_SOURCES = gettingstarted.cc
gettingstarted_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
gettingstarted_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

traversal_SOURCES = traversal.cc
traversal_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
traversal_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

integration_SOURCES = integration.cc
integration_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
integration_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

othergrids_SOURCES = othergrids.cc
othergrids_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
othergrids_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

adaptiveintegration_SOURCES = adaptiveintegration.cc
adaptiveintegration_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
adaptiveintegration_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

finitevolume_SOURCES = finitevolume.cc
finitevolume_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
finitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

adaptivefinitevolume_SOURCES = adaptivefinitevolume.cc
adaptivefinitevolume_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
adaptivefinitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

parfinitevolume_SOURCES = parfinitevolume.cc
parfinitevolume_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
parfinitevolume_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

visualization_SOURCES = visualization.cc
visualization_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
visualization_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

finiteelements_SOURCES = finiteelements.cc
finiteelements_CXXFLAGS = $(DUNEMPICPPFLAGS) $(ALL_PKG_CPPFLAGS)
finiteelements_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

# don't follow the full GNU-standard
# we need automake 1.9 or newer
AUTOMAKE_OPTIONS = foreign 1.9

# we need the module file to be able to build via dunecontrol
EXTRA_DIST= CMakeLists.txt dune.module

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid="$(DUNE_GRID_ROOT)" --with-dune-istl="$(DUNE_ISTL_ROOT)" --with-dune-disc="$(DUNE_DISC_ROOT)" CXX="$(CXX)" CC="$(CC)" --enable-parallel=@ENABLE_PARALLEL@ MPICC="$(MPICC)" --disable-documentation

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
