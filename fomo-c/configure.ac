m4_define(FOMOVERSION,3.1)
AC_INIT(FoMo,FOMOVERSION,tom.vandoorsselaere@wis.kuleuven.be)
AC_SUBST(fomoversion,FOMOVERSION)
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
AX_CXX_COMPILE_STDCXX_11
AX_LIB_HDF5([serial])
AC_CHECK_PROG(H5CPP,h5c++,`which h5c++`)
if test -z "${H5CPP}"; then
    CXXCPP="[h5c++]"
    CXX="[h5c++]"
    CPP="[h5c++]"
else
    CXXCPP="[g++]"
fi
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([sstream])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([cmath])
AC_CHECK_HEADERS([cstdlib])
AC_CHECK_HEADERS([iomanip])
AC_CHECK_HEADERS([map])
AC_CHECK_HEADERS([cassert])

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

AC_ARG_WITH([tbb],
	    [AS_HELP_STRING([--with-tbb],
			    [link with Intel Thread Building Block (TBB) for parallelisation @<:@default=yes@:>@])],
	    [tbb=$withval],
	    [tbb=yes])
AS_CASE(["$tbb"],
	[yes], [PKG_CHECK_MODULES([tbb], [tbb], [HAVE_TBB=1])],
	[no], [HAVE_TBB=0],
	[PKG_CHECK_MODULES([tbb], [tbb], [HAVE_TBB=1], [HAVE_TBB=0])])
AS_IF([test "$tbb" = yes -a "$HAVE_TBB" -eq 1], [PKG_CHECK_MODULES([tbbmalloc], [tbbmalloc], [HAVE_TBB=1], [HAVE_TBB=0])])
AS_IF([test "$tbb" = yes -a "$HAVE_TBB" -eq 1], [AC_DEFINE([CGAL_LINKED_WITH_TBB], [1], [Define if using TBB.])
CXXFLAGS+=" $tbb_CFLAGS $tbbmalloc_CFLAGS"
LDFLAGS+=" $tbb_LIBS $tbbmalloc_LIBS"
])

# generate documentation by doxygen
DX_DOXYGEN_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_MAN_FEATURE(ON)
DX_PDF_FEATURE(ON)
DX_INIT_DOXYGEN($PACKAGE_NAME,docfiles/fomo-doxygen.cfg,doc)

# write the Makefiles
AC_OUTPUT(Makefile src/Makefile example/Makefile example_FLASH/Makefile)
