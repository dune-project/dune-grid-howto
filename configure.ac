#! /bin/bash
# $Id$
# Process this file with autoconf to produce a configure script.
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([./gettingstarted.cc])
AC_CONFIG_HEADERS([config.h])

AC_CHECK_PROGS([TEX], [latex], [true])
AC_CHECK_PROGS([BIBTEX], [bibtex], [true])
AC_CHECK_PROGS([DVIPDF], [dvipdf], [true])
AC_CHECK_PROGS([DVIPS], [dvips], [true])
AC_CHECK_PROGS([WML], [wml], [true])
AM_CONDITIONAL([WML], [test "x$WML" != xtrue])
AC_CHECK_PROGS([TEX4HT], [tex4ht], [true])
AC_CHECK_PROGS([MK4HT], [mk4ht], [true])
AC_CHECK_PROGS([T4HT], [t4ht], [true])
AM_CONDITIONAL([TEX4HT], [test "x$TEX4HT" != xtrue])
AC_CHECK_PROGS([CONVERT], [convert], [false])
AM_CONDITIONAL([CONVERT], [test "x$CONVERT" != xfalse])

# special variable to include the documentation into the website
AC_ARG_WITH(duneweb,
    AC_HELP_STRING([--with-duneweb=PATH],[Only needed for website-generation, path to checked out version of dune-web]))

if test x$with_duneweb != x ; then
   # parameter is set. Check it
   AC_MSG_CHECKING([whether passed Dune-Web directory appears correct])
   WEBTESTFILE="$with_duneweb/layout/default.wml"
   if test -d "$with_duneweb" && test -e "$WEBTESTFILE" ; then
      AC_MSG_RESULT([ok])
      # normalize path
      with_duneweb=`(cd $with_duneweb && pwd)` ;
   else
      if test -d "$with_duneweb" ; then
        AC_MSG_ERROR([$WEBTESTFILE not found in Dune-web dir $with_duneweb!])
      else
        AC_MSG_ERROR([Dune-Web directory $with_duneweb not found!])
      fi
   fi
fi
AC_SUBST(DUNEWEBDIR, $with_duneweb)

# check all dune-module stuff
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')

# write output
AC_CONFIG_FILES([Makefile
        doc/Makefile
        doc/config.inc
        grids/Makefile
        dune-grid-howto.pc])
AC_OUTPUT

DUNE_SUMMARY_ALL
