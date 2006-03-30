# $id: $

SUBDIRS =

DOCS = grid-howto.dvi grid-howto.pdf grid-howto.ps

dist_pkgdata_DATA = $(DOCS)

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += grid-howto*.html grid-howto*.png

# setting like in dune-web
#CURDIR=doc/

%.html: %.dvi
	mk4ht htlatex $*.tex "xhtml,2,info,png"
	t4ht $*.tex

include $(top_srcdir)/am/latex
