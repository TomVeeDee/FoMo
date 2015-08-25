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

# check for HDF5 for the hdf5 example
AX_LIB_HDF5([serial])
if test "x$with_hdf5" = "xyes"; then
   AC_CHECK_PROGS(HDF5_CXX,[h5c++ h5c++-64])
   if [test -z "$HDF5_CXX"]; then
	   AC_MSG_ERROR([h5c++ not found])
   else
	   HDF5_SHOW=$(eval $HDF5_CXX -show -shlib)
           HDF5_tmp_flags=$(eval $HDF5_CXX -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

           HDF5_tmp_inst=$(eval $HDF5_CXX -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

           HDF5_CPPFLAGS="-I${HDF5_tmp_inst}/include"

           dnl Now sort the flags out based upon their prefixes
           for arg in $HDF5_SHOW $HDF5_tmp_flags ; do
             case "$arg" in
               -I*) echo $HDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                     || HDF5_CPPFLAGS="$arg $HDF5_CPPFLAGS"
                 ;;
               -L*) echo $HDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                     || HDF5_LDFLAGS="$arg $HDF5_LDFLAGS"
                ;;
               -l*) echo $HDF5_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                     || HDF5_LIBS="$arg $HDF5_LIBS"
                 ;;
             esac
           done

           HDF5_LIBS="$HDF5_LIBS -lhdf5"
   fi
fi
AM_CONDITIONAL(HAVE_HDF5, [test -n "$HDF5_CXX"])

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
AC_OUTPUT(Makefile src/Makefile example/Makefile example/example_FLASH/Makefile)
