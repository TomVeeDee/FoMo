export fomoversion=@fomoversion@
ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src example
EXTRA_DIST = idl docfiles
@DX_RULES@
AM_CPPFLAGS = -I$(VTK_INCLUDE_DIR)
LDADD = -L$(VTK_LIB_DIR) -lvtkCommonCore-9.3 -lvtkIOCore-9.3 -lvtksys-9.3 
