m4_define(FOMOVERSION,3.3)
AC_INIT(FoMo,FOMOVERSION,tom.vandoorsselaere@kuleuven.be)
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
AC_CHECK_HEADERS([climits])
AC_CHECK_HEADERS([bitset chrono])

# Checks for OpenCL
# Copied from CGAL code above
AC_ARG_WITH([opencl],
	    [AS_HELP_STRING([--with-opencl=path/yes/no],
			    [OpenCL for GPU-accelerated rendering (found in non-default location path) @<:@default=yes@:>@])],
	    [opencl=$withval],
	    [opencl=yes])
link_opencl=false
if test "a$opencl" != "ano" ; then
	CLDIRS="/usr/lib64/nvidia /usr/lib /usr/lib64 /usr/lib/nvidia "
	if test "a$opencl" != "ayes" ; then
		CPPFLAGS+=" -I$opencl/include"
		LDFLAGS+=" -L$opencl/lib -Wl,-rpath -Wl,$opencl/lib"
		CLDIRS+=" $opencl/lib $opencl/lib64"
	fi
	AC_CHECK_HEADERS([CL/cl.hpp],[],AC_MSG_ERROR([OpenCL header not found]))
	# headers have been found successfully, now check the linking
	tmp_LDFLAGS=$LDFLAGS
	for CLDIR in $CLDIRS; do
		LDFLAGS+=" -L$CLDIR -Wl,-rpath -Wl,$CLDIR"
		AC_SEARCH_LIBS(clCompileProgram,[OpenCL nvidia-opencl],[
						 link_opencl=true
						 ],[
						    unset ac_cv_search_clCompileProgram
						    LDFLAGS=$tmp_LDFLAGS])
		if $link_opencl; then 
			break
		fi
	done
	if ! $link_opencl; then
		AC_MSG_ERROR([OpenCL libraries not found])
	fi
fi

# checks for building gui, depends on QT5
AC_ARG_WITH([qt5],
	    [AS_HELP_STRING([--with-qt5=yes/no],
			    [Enable Qt5 for GUI (add correct qmake to PATH) @<:@default=yes@:>@])],
	    [qt5=$withval],
	    [qt5=yes])
have_qt=no
if test "a$qt5" = "ayes"; then
	# fail if the OpenCL libraries have not been found
	if ! $link_opencl; then
		AC_MSG_ERROR([OpenCL needed for GUI])
	fi
	AX_HAVE_QT
	if test "a$have_qt" = "ano" ; then
		AC_MSG_ERROR([Qt5 not found])
	fi
fi
AC_DEFINE(HAVE_QT5,,[Qt5 libraries are available])
AM_CONDITIONAL(HAVE_QT5,[test "a$qt5" = "ayes" -a "a$have_qt" = "ayes"])

AX_BOOST_BASE([1.54])
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
AC_CHECK_HEADER(boost/progress.hpp,[],[AC_MSG_ERROR(Boost progress header not found)])
AC_CHECK_HEADER(boost/geometry.hpp,[],[AC_MSG_ERROR(Boost geometry header not found)])
AC_CHECK_HEADER(boost/geometry/index/rtree.hpp,[],[AC_MSG_ERROR(Boost rtree header not found)],[#include <iostream>])
AC_CHECK_HEADER(boost/regex.hpp,[],[AC_MSG_ERROR(Boost regex header not found)])
AC_CHECK_HEADER(boost/algorithm/string/replace.hpp,[],[AC_MSG_ERROR(Boost replace header not found)])
AC_CHECK_HEADER(boost/filesystem/path.hpp,[],[AC_MSG_ERROR(Boost path header not found)])
AC_CHECK_HEADER(boost/iostreams/filtering_streambuf.hpp,[],[AC_MSG_ERROR(Boost filtering_streambuf header not found)])
AC_CHECK_HEADER(boost/iostreams/copy.hpp,[],[AC_MSG_ERROR(Boost copy header not found)])
AC_CHECK_HEADER(boost/iostreams/filter/gzip.hpp,[],[AC_MSG_ERROR(Boost zlib header not found)])

PKG_PROG_PKG_CONFIG



#check if MPI is available
if test "a$mpi" = "ayes" ; then
   AC_CHECK_HEADER(mpi.h, [AC_DEFINE([HAVEMPI],"1",[Is MPI used?])],AC_MSG_ERROR([MPI header not found]))
   AC_SEARCH_LIBS(MPI_Init,[],[],AC_MSG_ERROR([MPI libraries not found]))
fi

AC_CHECK_HEADER([gsl/gsl_const_mksa.h],[],AC_MSG_ERROR([GSL physical constants not found]))
AC_ARG_WITH([cgal],
	    [AS_HELP_STRING([--with-cgal=yes/no],
			    [link with CGAL for Delaunay triangulation for interpolation @<:@default=no@:>@])],
	    [cgal=$withval],
	    [cgal=no])
if test "a$cgal" = "ayes" ; then
	AC_CHECK_HEADERS([CGAL/Delaunay_triangulation_2.h],[],AC_MSG_ERROR([CGAL Delaunay headers not found]))
	#AC_SEARCH_LIBS(a_function_what?,[CGAL],[],AC_MSG_ERROR([CGAL libraries not found]))
	LDFLAGS+=" -lCGAL -lgmp"
	CXXFLAGS+=" -frounding-math"
fi

AC_ARG_WITH([tbb],
	    [AS_HELP_STRING([--with-tbb=yes/no/PATH],
			    [link with Intel Thread Building Block (TBB) for parallelisation of CGAL triangulation @<:@default=no@:>@])],
	    [tbb=$withval],
	    [tbb=no])
HAVE_TBB=0
AS_IF([test "a$tbb" != ano -a "$HAVE_TBB" -eq 0 -a "a$cgal" = "ayes"], 
      # pkg_config did not find the TBB yet, it is probably in a special dir
      [
       AS_IF([test "a$tbb" = ayes],tbbdirs="/usr /usr/local /opt /opt/local",tbbdirs="$withval")
       found=false
       if ! $found; then
         tbb_CPPFLAGS=
         for tbbdir in $tbbdirs; do
             AC_MSG_CHECKING([for tbb/tbb.h in $tbbdir])
             if test -f "$tbbdir/include/tbb/tbb.h"; then
                tbb_CPPFLAGS="-I$tbbdir/include"
                tbb_LDFLAGS="-L$tbbdir/lib"
                tbb_LIBS="-ltbb -ltbbmalloc"
                found=true
		HAVE_TBB=1
                AC_MSG_RESULT([yes])
                break
            else
                AC_MSG_RESULT([no])
            fi
         done
       # as a last resort, check if the TBB may be found with pkg-config
         if ! $found; then 
	       PKG_CHECK_MODULES([tbb],[tbb tbbmalloc], [HAVE_TBB=1])
	 fi
       fi
       AC_MSG_CHECKING([whether compiling and linking against tbb works])
       echo "Trying link with tbb_LDFLAGS=$tbb_LDFLAGS;" \
        "tbb_LIBS=$tbb_LIBS; tbb_CPPFLAGS=$tbb_CPPFLAGS" >&AS_MESSAGE_LOG_FD

       save_LIBS="$LIBS"
       save_LDFLAGS="$LDFLAGS"
       save_CPPFLAGS="$CPPFLAGS"
       LDFLAGS="$LDFLAGS $tbb_LDFLAGS"
       LIBS="$tbb_LIBS $LIBS"
       CPPFLAGS="$tbb_CPPFLAGS $CPPFLAGS"

       AC_LINK_IFELSE(
           [AC_LANG_PROGRAM([#include <tbb/tbb.h>],
                           [tbb::concurrent_hash_map<int,int>()])],
           [ 
               AC_MSG_RESULT([yes])
           ], [
               AC_MSG_RESULT([no])
           ])
    CPPFLAGS="$save_CPPFLAGS"
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS" 
      ]
      )
AS_IF([test "a$tbb" != ano -a "$HAVE_TBB" -eq 1], [AC_DEFINE([CGAL_LINKED_WITH_TBB], [1], [Define if using TBB.])
CXXFLAGS+=" $tbb_CFLAGS"
LDFLAGS+=" $tbb_LIBS"
LDFLAGS="$LDFLAGS $tbb_LDFLAGS"
LIBS="$tbb_LIBS $LIBS"
CPPFLAGS="$tbb_CPPFLAGS $CPPFLAGS"
])

# generate documentation by doxygen
DX_DOXYGEN_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_MAN_FEATURE(ON)
DX_PDF_FEATURE(ON)
DX_INIT_DOXYGEN($PACKAGE_NAME,docfiles/fomo-doxygen.cfg,doc)

# write the Makefiles
AC_OUTPUT(Makefile src/Makefile example/Makefile example/example_FLASH/Makefile example/example_mpi_amrvac/Makefile gui/fomo-gui/Makefile)
