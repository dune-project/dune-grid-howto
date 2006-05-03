# $id: $

SUBDIRS =

MAINDOCS = grid-howto.dvi grid-howto.pdf grid-howto.ps

if WML
 WMLDOCS = index.html
endif

if TEX4HT
 HTDOCS =  navigation.html css.html grid-howto.html
endif

if BUILD_DOCS
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif

dist_pkgdata_DATA = $(DOCS)

EXTRAINSTALL = grid-howto*.html *.png *.css

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += grid-howto*.html grid-howto*.png
CLEANFILES += $(DOCS)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/grid-howto

grid-howto.html: grid-howto.tex navigation.html css.html
	$(MAKE) grid-howto.dvi
	$(MK4HT) htlatex grid-howto.tex "dune"; \
	$(T4HT) grid-howto.tex
	$(TEX) grid-howto.tex

include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
