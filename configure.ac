#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT(dune-grid-howto, 0.0, dune@hal.iwr.uni-heidelberg.de)
AM_INIT_AUTOMAKE(dune-grid-howto, 0.0, dune@hal.iwr.uni-heidelberg.de)
AC_CONFIG_SRCDIR([gettingstarted.cc])
AM_CONFIG_HEADER([config.h])

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
