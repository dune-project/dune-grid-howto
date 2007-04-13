# $id: $

EPS_FILES = EPS/adaptiveintegration_alberta2d.eps EPS/adaptiveintegration_ug3d.eps \
    EPS/alberta2d.eps EPS/alberta2d-view1.eps EPS/alberta2d-view2.eps EPS/alberta3d.eps \
    EPS/alucube3d.eps EPS/alusimplex3d.eps EPS/conc.eps EPS/iso.eps EPS/ug2dquad-view1.eps \
    EPS/ug2dquad-view2.eps EPS/ug2dtri-view1.eps EPS/ug2dtri-view2.eps EPS/ugcube2d.eps \
    EPS/ugcube3d.eps EPS/ugsimplex2d.eps EPS/ugsimplex3d.eps EPS/yasp3d.eps

DOCSOURCE = grid-howto.tex navigation.html css.html

MAINDOCS = grid-howto.dvi grid-howto.pdf grid-howto.ps

if WML
 WMLDOCS = index.html
endif

if TEX4HT
 HTDOCS =  navigation.html css.html grid-howto.html
endif

# EXTRA_DIST = $(DOCSOURCE) $(WMLDOCS) $(HTDOCS)
EXTRA_DIST = grid-howto.pdf

if BUILD_DOCS
if CONVERT
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/grid-howto

docdir=$(datadir)/doc/dune-grid-howto

doc_DATA = grid-howto.pdf

EXTRAINSTALL = grid-howto.pdf grid-howto*.html *.png *.css

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += grid-howto*.html grid-howto*.png
CLEANFILES += $(DOCS)
if CONVERT
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	convert $< $@
grid-howto.tex: $(EPS_FILES)
endif

grid-howto.html:  $(DOCSOURCE)
	$(MAKE) grid-howto.dvi
	$(MK4HT) htlatex grid-howto.tex "dune"; \
	$(T4HT) grid-howto.tex
	$(TEX) grid-howto.tex

dist-hook:
	cp $(srcdir)/Makefile.dist.am $(distdir)/Makefile.am
	cp $(srcdir)/Makefile.dist.in $(distdir)/Makefile.in

EXTRA_TEXINPUTS:=$(top_srcdir)

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
