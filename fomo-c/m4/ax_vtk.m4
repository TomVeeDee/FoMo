AC_DEFUN([AX_VTK], [
    AC_MSG_CHECKING([for VTK])

    # Allow user to specify the VTK installation directory
    AC_ARG_WITH([vtk-dir],
        [AS_HELP_STRING([--with-vtk-dir=DIR],
                        [Specify VTK install directory (default: /usr/local)])],
        [vtk_dir=$withval],
        [vtk_dir=""])

    # Check for VTK include directory
    AC_MSG_CHECKING([for VTK include directory])
    VTK_INCLUDE_DIR=""
    for d in \
        /usr/include/vtk* \
        /usr/local/include/vtk* \
        if test -d "$d"; then
            VTK_INCLUDE_DIR="$d"
            break
        fi
    done
    if test -d "$VTK_INCLUDE_DIR"; then
        AC_MSG_RESULT([$VTK_INCLUDE_DIR])
    elif test -d $vtk_dir/include/vtk"; then
        VTK_INCLUDE_DIR="$vtk_dir/include/vtk"
        AC_MSG_RESULT([$VTK_INCLUDE_DIR])
    else
        AC_MSG_ERROR([Cannot find VTK include directory in $VTK_INCLUDE_DIR. Please check your VTK installation or use --with-vtk-dir.])
    fi

    # Check for VTK library directory
    AC_MSG_CHECKING([for VTK library directory])
    VTK_LIB_DIR="$vtk_dir/lib64"
    if test -d "$VTK_LIB_DIR"; then
        AC_MSG_RESULT([$VTK_LIB_DIR])
    else
        # fallback if lib64 not found
        VTK_LIB_DIR="$vtk_dir/lib"
        if test -d "$VTK_LIB_DIR"; then
            AC_MSG_RESULT([$VTK_LIB_DIR])
        else
            AC_MSG_ERROR([Cannot find VTK library directory in $vtk_dir.])
        fi
    fi

    # Export variables for Makefile substitution
    AC_SUBST([VTK_INCLUDE_DIR])
    AC_SUBST([VTK_LIB_DIR])

    AC_MSG_NOTICE([VTK found in $vtk_dir])
])