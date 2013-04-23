# $id: $

set(EPS_FILES  EPS/adaptiveintegration_alberta2d.eps EPS/adaptiveintegration_ug3d.eps
    EPS/alberta2d.eps EPS/alberta2d-view1.eps EPS/alberta2d-view2.eps EPS/alberta3d.eps
    EPS/alucube3d.eps EPS/alusimplex3d.eps EPS/conc.eps EPS/fem2d.eps EPS/fem3d.eps EPS/iso.eps EPS/ug2dquad-view1.eps
    EPS/ug2dquad-view2.eps EPS/ug2dtri-view1.eps EPS/ug2dtri-view2.eps EPS/ugcube2d.eps
    EPS/ugcube3d.eps EPS/ugsimplex2d.eps EPS/ugsimplex3d.eps EPS/yasp3d.eps)


set(MAINDOCS  grid-howto.pdf)

if(WML_FOUND)
 set(DOCSOURCE  grid-howto.tex navigation.html css.html)
 set(HTDOCS  navigation.html css.html)
if(TEX4HT_FOUND)
# HTDOCS += grid-howto.html
endif(TEX4HT_FOUND )
endif(WML_FOUND )


if( BUILD_DOCS)
if(CONVERT_FOUND)
 set(DOCFILES  ${HTDOCS} ${MAINDOCS} ${WMLDOCS})
endif(CONVERT_FOUND )
endif( BUILD_DOCS )

set(EXTRA_DIST  ${DOCFILES} grid-howto.tex ${EPS_FILES})

# setting like in dune-web
set(BASEDIR ../..)
set(CURDIR doc/grid-howto)

set(docdir ${CMAKE_INSTALL_DOCDIR}/dune-grid-howto)

set(EXTRAINSTALL  grid-howto.pdf grid-howto*.html *.png *.css)

# Install documentation
# We assume that ${DOCFILES} and ${DOCFILES_EXTRA} are targets
# TODO: care about the following line
### install(TARGETS "${DOCFILES}" ${DOCFILES_EXTRA} DESTINATION ${docdir})

set(SOURCES
  ../adaptivefinitevolume.cc
  ../adaptiveintegration.cc
  ../basicunitcube.hh
  ../elementdata.hh
  ../evolve.hh
  ../finiteelements.cc
  ../finitevolumeadapt.hh
  ../finitevolume.cc
  ../functors.hh
  ../gettingstarted.cc
  ../initialize.hh
  ../integrateentity.hh
  ../integration.cc
  ../othergrids.cc
  ../parevolve.hh
  ../parfinitevolume.cc
  ../parfvdatahandle.hh
  ../shapefunctions.hh
  ../transportproblem2.hh
  ../transportproblem.hh
  ../traversal.cc
  ../unitcube_albertagrid.hh
  ../unitcube_alugrid.hh
  ../unitcube.hh
  ../unitcube_sgrid.hh
  ../unitcube_yaspgrid.hh
  ../vertexdata.hh
  ../visualization.cc
  ../vtkout.hh
)

set(IMAGES
  EPS/adaptiveintegration_alberta2d.jpg
  EPS/adaptiveintegration_ug3d.jpg
  EPS/alberta2d.jpg
  EPS/alberta2d-view1.jpg
  EPS/alberta2d-view2.jpg
  EPS/alberta3d.jpg
  EPS/alucube3d.jpg
  EPS/alusimplex3d.jpg
  EPS/conc.jpg
  EPS/fem2d.jpg
  EPS/fem3d.jpg
  EPS/intersection.eps
  EPS/islocalref.eps
  EPS/iso.jpg
  EPS/partitionsingle.eps
  EPS/ug2dquad-view1.jpg
  EPS/ug2dquad-view2.jpg
  EPS/ug2dtri-view1.jpg
  EPS/ug2dtri-view2.jpg
  EPS/ugcube2d.jpg
  EPS/ugcube3d.jpg
  EPS/ugsimplex2d.jpg
  EPS/ugsimplex3d.jpg
  EPS/yasp3d.jpg
)

# create file config.tex containing the version of this module
set(PACKAGE_VERSION ${DUNE_MOD_VERSION})
configure_file(config.inc.in config.inc)

# Install documentation
# We assume that ${DOCFILES} and ${DOCFILES_EXTRA} are targets
dune_add_latex_document(grid-howto.tex
  FATHER_TARGET doc
  DEFAULT_PDF
  INPUTS grid-howto.tex config.inc ${SOURCES}
  IMAGES ${IMAGES})

