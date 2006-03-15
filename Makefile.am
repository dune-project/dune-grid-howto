# $Id$


noinst_PROGRAMS = gettingstarted traversal integration othergrids adaptiveintegration

gettingstarted_SOURCES = gettingstarted.cc
gettingstarted_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
gettingstarted_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

traversal_SOURCES = traversal.cc
traversal_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
traversal_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

integration_SOURCES = integration.cc
integration_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
integration_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

othergrids_SOURCES = othergrids.cc
othergrids_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
othergrids_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

adaptiveintegration_SOURCES = adaptiveintegration.cc
adaptiveintegration_CXXFLAGS = $(MPI_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
adaptiveintegration_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
