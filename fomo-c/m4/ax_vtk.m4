AC_DEFUN([AX_VTK], [
    AC_MSG_CHECKING([for VTK])

    # Allow user to specify the VTK installation directory
    AC_ARG_WITH([vtk-dir],
        [AS_HELP_STRING([--with-vtk-dir=DIR], [Specify VTK install directory])],
        [vtk_dir=$with_val],
        [vtk_dir=$HOME/Darias/code/VTK-install/usr/local])

    # Check for VTK include directory
    AC_MSG_CHECKING([for VTK include directory])
    VTK_INCLUDE_DIR=$HOME/Darias/code/VTK-install/usr/local/include/vtk-9.3
    AC_MSG_NOTICE([VTK_INCLUDE_DIR is: $VTK_INCLUDE_DIR]) 
    if test -d "$VTK_INCLUDE_DIR"; then
        AC_MSG_RESULT([found $VTK_INCLUDE_DIR])
    else
        AC_MSG_ERROR([Cannot find VTK include directory. Please check the VTK installation.])
    fi

    # Check for VTK library directory
    AC_MSG_CHECKING([for VTK library directory])
    VTK_LIB_DIR=$vtk_dir/lib64
    if test -d "$VTK_LIB_DIR"; then
        AC_MSG_RESULT([found $VTK_LIB_DIR])
    else
        AC_MSG_ERROR([Cannot find VTK library directory. Please check the VTK installation.])
    fi

    # Set the flags for compilation
    AC_SUBST([VTK_INCLUDE_DIR])
    AC_SUBST([VTK_LIB_DIR])

    AC_MSG_RESULT([found])
])
