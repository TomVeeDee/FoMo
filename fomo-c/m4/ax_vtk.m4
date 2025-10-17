AC_DEFUN([AX_VTK], [
    AC_MSG_CHECKING([for VTK])

    # Allow user to specify the VTK installation directory (optional)
    AC_ARG_WITH([vtk-dir],
        [AS_HELP_STRING([--with-vtk-dir=DIR], [Specify VTK install directory if not in a standard path])],
        [vtk_dir=$withval],
        [vtk_dir=""])

    # Check for VTK include directory
    AC_MSG_CHECKING([for VTK include directory])
    VTK_INCLUDE_DIR=""

    if test -d "/usr/include/vtk" ; then
        VTK_INCLUDE_DIR="/usr/include/vtk"
    elif ls -d /usr/include/vtk-* >/dev/null 2>&1; then
        VTK_INCLUDE_DIR=$(ls -d /usr/include/vtk-* | head -n1)
    elif test -d "/usr/local/include/vtk" ; then
        VTK_INCLUDE_DIR="/usr/local/include/vtk"
    elif ls -d /usr/local/include/vtk-* >/dev/null 2>&1; then
        VTK_INCLUDE_DIR=$(ls -d /usr/local/include/vtk-* | head -n1)
    elif test -n "$vtk_dir" && test -d "$vtk_dir/include"; then
        VTK_INCLUDE_DIR="$vtk_dir/include"
    fi

    if test -n "$VTK_INCLUDE_DIR"; then
        AC_MSG_RESULT([$VTK_INCLUDE_DIR])
    else
        AC_MSG_ERROR([VTK not found in standard locations. Please specify with --with-vtk-dir=DIR])
    fi

    # Check for VTK library directory
    AC_MSG_CHECKING([for VTK library directory])
    if test -n "$vtk_dir" && test -d "$vtk_dir/lib64"; then
        VTK_LIB_DIR="$vtk_dir/lib64"
    elif test -n "$vtk_dir" && test -d "$vtk_dir/lib"; then
        VTK_LIB_DIR="$vtk_dir/lib"
    elif test -d "/usr/lib64"; then
        VTK_LIB_DIR="/usr/lib64"
    elif test -d "/usr/local/lib64"; then
        VTK_LIB_DIR="/usr/local/lib64"
    elif test -d "/usr/lib"; then
        VTK_LIB_DIR="/usr/lib"
    elif test -d "/usr/local/lib"; then
        VTK_LIB_DIR="/usr/local/lib"
    else
        AC_MSG_ERROR([VTK libraries not found. Please specify with --with-vtk-dir=DIR])
    fi
    AC_MSG_RESULT([$VTK_LIB_DIR])

    # Export variables for Makefile substitution
    AC_SUBST([VTK_INCLUDE_DIR])
    AC_SUBST([VTK_LIB_DIR])

    AC_MSG_NOTICE([VTK found: includes in $VTK_INCLUDE_DIR, libs in $VTK_LIB_DIR])
])