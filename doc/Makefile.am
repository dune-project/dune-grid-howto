# $id: $

EPS_FILES = EPS/adaptiveintegration_alberta2d.eps EPS/adaptiveintegration_ug3d.eps \
    EPS/alberta2d.eps EPS/alberta2d-view1.eps EPS/alberta2d-view2.eps EPS/alberta3d.eps \
    EPS/alucube3d.eps EPS/alusimplex3d.eps EPS/conc.eps EPS/fem2d.eps EPS/fem3d.eps EPS/iso.eps EPS/ug2dquad-view1.eps \
    EPS/ug2dquad-view2.eps EPS/ug2dtri-view1.eps EPS/ug2dtri-view2.eps EPS/ugcube2d.eps \
    EPS/ugcube3d.eps EPS/ugsimplex2d.eps EPS/ugsimplex3d.eps EPS/yasp3d.eps


MAINDOCS = grid-howto.pdf

if WML
 DOCSOURCE = grid-howto.tex navigation.html css.html
 HTDOCS = navigation.html css.html
if TEX4HT
# HTDOCS += grid-howto.html
endif
endif


if BUILD_DOCS
if CONVERT
 DOCFILES = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif

EXTRA_DIST = CMakeLists.txt $(DOCFILES) grid-howto.tex $(EPS_FILES)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/grid-howto

docdir=$(datadir)/doc/dune-grid-howto

EXTRAINSTALL = grid-howto.pdf grid-howto*.html *.png *.css

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.ilg *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += grid-howto*.html grid-howto*.png navigation.html config.tex
CLEANFILES += $(DOCS)
if CONVERT
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	convert $< $@
grid-howto.tex: $(EPS_FILES)
endif

if WML
if TEX4HT
grid-howto.html:  $(DOCSOURCE) dune.cfg tex4ht.env
	$(MAKE) grid-howto.dvi
	$(MK4HT) htlatex grid-howto.tex "dune"; \
	$(MK4HT) htlatex grid-howto.tex "dune"; \
	$(T4HT) grid-howto.tex
	$(TEX) grid-howto.tex
endif
endif

# if dune-web is not used, we don't have any navigation to include
# create dummy file.
# depend on Makefile to make sure the file is recreated
# after enabling to dune-web
navigation.html: Makefile
if DUNEWEB
else
navigation.html:
	echo > navigation.html
endif

EXTRA_TEXINPUTS:=$(top_srcdir)

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
