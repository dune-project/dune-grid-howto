#! /bin/bash
# $Id$
# Process this file with autoconf to produce a configure script.
AC_INIT(dune-grid-howto, 0.0, dune@hal.iwr.uni-heidelberg.de)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([./gettingstarted.cc])
AM_CONFIG_HEADER([config.h])

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
DUNE_CHECK_ALL([dunecommon], [dunegrid])

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')

echo
echo Note: Most of the libraries checked for above are only used for the self-test
echo of Dune.  The library itself will build and the headers will work even if
echo Alberta, MPI, etc. cannot be found.  An exception to this are UG and AmiraMesh
echo which need to be found right now if you want to use them later.
echo

# write output
AC_CONFIG_FILES([Makefile doc/Makefile bin/wmlwrap])
AC_OUTPUT

chmod +x bin/wmlwrap

DUNE_SUMMARY_ALL
