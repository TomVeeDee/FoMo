AC_INIT(visual,0.99)
AC_PREREQ(2.59)
AM_INIT_AUTOMAKE(1.9)
AC_CONFIG_HEADERS(config.h)
# check if mpi is wanted
AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi],[build MPI version of POLLUX])],
mpi=$withval,mpi="no")
# search for the fortran compiler
if test "a$mpi" = "ayes" ; then
   CXXCOMP="[mpicxx]"
else
   CXXCOMP="[icc]"
fi
AC_PROG_CXX($CXXCOMP)
AC_LANG(C++)
PKG_PROG_PKG_CONFIG
#check if MPI is available
if test "a$mpi" = "ayes" ; then
   AC_CHECK_HEADER(mpi.h, [AC_DEFINE([HAVEMPI],"1",[Is MPI used?])],AC_MSG_ERROR([MPI header not found]))
   AC_SEARCH_LIBS(MPI_Init,[],[],AC_MSG_ERROR([MPI libraries not found]))
fi
# check if png is wanted
AC_ARG_WITH([png],
  [AS_HELP_STRING([--with-png],
    [support png files @<:@default=yes@:>@])],
  [png=$withval],
  [png=yes])
AS_CASE(["$png"],
  [yes], [PKG_CHECK_MODULES([libpng], [libpng], [HAVE_LIBPNG=1])],
  [no], [HAVE_LIBPNG=0],
  [PKG_CHECK_MODULES([libpng], [libpng], [HAVE_LIBPNG=1], [HAVE_LIBPNG=0])])
# use this line to conditionally include code
#AM_CONDITIONAL([USE_LIBPNG], [test "$png" != no -a "$HAVE_LIBPNG" -eq 1])
AS_IF([test "$png" = yes -a "$HAVE_LIBPNG" -eq 1], [AC_DEFINE([HAVEPNG], [1], [Define if using libpng.])
CXXFLAGS+=" $libpng_CFLAGS"
LDFLAGS+=" $libpng_LIBS"
])

# check if mpeg is wanted
AC_ARG_WITH([mpeg],
  [AS_HELP_STRING([--with-mpeg],
    [support mpeg files @<:@default=yes@:>@])],
  [mpeg=$withval],
  [mpeg=yes])
AS_CASE(["$mpeg"],
  [yes], [PKG_CHECK_MODULES([libavcodec], [libavcodec], [HAVE_LIBAVCODEC=1])],
  [no], [HAVE_LIBPNG=0],
  [PKG_CHECK_MODULES([libavcodec], [libavcodec], [HAVE_LIBAVCODEC=1], [HAVE_LIBAVCODEC=0])])
# use this line to conditionally include code
#AM_CONDITIONAL([USE_LIBAVCODEC], [test "$mpeg" != no -a "$HAVE_LIBAVCODEC" -eq 1])
AS_IF([test "$mpeg" = yes -a "$HAVE_LIBAVCODEC" -eq 1], [AC_DEFINE([HAVEMPEG], [1], [Define if using ffmpeg.])
CXXFLAGS+=" $libavcodec_CFLAGS"
LDFLAGS+=" $libavcodec_LIBS"
])

#if test "a$mpeg" = "ayes" ; then
#   AC_CHECK_HEADER(avcodec.h, [],AC_MSG_ERROR([MPEG header not found]))
#   AC_SEARCH_LIBS(avcodec_init,[avcodec],[],AC_MSG_ERROR([MPEG libraries not found]))
#fi

# write the Makefiles
AC_OUTPUT(Makefile src/Makefile)
