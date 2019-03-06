add_subdirectory(EPS)

# create file config.inc containing the version of this module
set(PACKAGE_VERSION ${DUNE_MOD_VERSION})
configure_file(config.inc.in config.inc)


dune_add_latex_document(
  SOURCE grid-howto.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR})
