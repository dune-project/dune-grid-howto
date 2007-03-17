# all eps files which are not part of the distribution
EPS_FILES = adaptiveintegration_alberta2d.eps adaptiveintegration_ug3d.eps \
    alberta2d.eps alberta2d-view1.eps alberta2d-view2.eps alberta3d.eps \
    alucube3d.eps alusimplex3d.eps conc.eps iso.eps ug2dquad-view1.eps \
    ug2dquad-view2.eps ug2dtri-view1.eps ug2dtri-view2.eps ugcube2d.eps \
    ugcube3d.eps ugsimplex2d.eps ugsimplex3d.eps yasp3d.eps

EXTRA_DIST = adaptiveintegration_alberta2d.jpg adaptiveintegration_ug3d.jpg \
    alberta2d.jpg alberta2d-view1.jpg alberta2d-view2.jpg alberta3d.jpg \
    alucube3d.jpg alusimplex3d.jpg conc.eps intersection.eps islocalref.eps \
    iso.jpg partitionsingle.eps partitionsingle.svg ug2dquad-view1.jpg \
    ug2dquad-view2.jpg ug2dtri-view1.jpg ug2dtri-view2.jpg ugcube2d.jpg \
    ugcube3d.jpg ugsimplex2d.jpg ugsimplex3d.jpg yasp3d.jpg

if CONVERT
CLEANFILES = $(EPS_FILES)

all-am: $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	convert $< $@
endif

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
