AC_INIT(visual,0.99)
AC_PREREQ(2.59)
AM_INIT_AUTOMAKE(1.9)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_PROG_CXX
AC_ENABLE_SHARED
AC_DISABLE_STATIC
LT_INIT
AC_PROG_LIBTOOL(libtool)
# check if mpi is wanted
AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi],[build MPI version of POLLUX])],
mpi=$withval,mpi="no")
# search for the fortran compiler
if test "a$mpi" = "ayes" ; then
   CXXCOMP="[mpicxx]"
else
   CXXCOMP="[g++]"
fi
AC_PROG_CXX($CXXCOMP)
AC_PROG_CXXCPP([cpp])
AC_LANG(C++)
AC_LANG_CPLUSPLUS
AX_CXX_COMPILE_STDCXX_0X
AX_BOOST_BASE([1.36])
AC_CHECK_HEADER(boost/progress.hpp,[],[AC_MSG_ERROR(Boost header not found)])
PKG_PROG_PKG_CONFIG
#check if MPI is available
if test "a$mpi" = "ayes" ; then
   AC_CHECK_HEADER(mpi.h, [AC_DEFINE([HAVEMPI],"1",[Is MPI used?])],AC_MSG_ERROR([MPI header not found]))
   AC_SEARCH_LIBS(MPI_Init,[],[],AC_MSG_ERROR([MPI libraries not found]))
fi

AC_CHECK_HEADER([gsl/gsl_const_mksa.h],[],AC_MSG_ERROR([GSL physical constants not found]))
AC_CHECK_HEADER([CGAL/Delaunay_triangulation_2.h],[],AC_MSG_ERROR([CGAL Delaunay headers not found]))
#AC_SEARCH_LIBS(a_function_what?,[CGAL],[],AC_MSG_ERROR([CGAL libraries not found]))
LDFLAGS+="-lCGAL -lgmp"
CXXFLAGS+=" -frounding-math"

# write the Makefiles
AC_OUTPUT(Makefile src/Makefile example/Makefile)
