if HAVE_HDF5
  noinst_PROGRAMS = example_hdf5
  LDADD = -L$(top_builddir)/src/.libs/ -lFoMo $(HDF5_LDFLAGS) $(HDF5_LIBS)
  AM_CPPFLAGS = -I$(top_srcdir)/src $(HDF5_CPPFLAGS)
  example_hdf5_SOURCES = example_hdf5.cpp
  EXTRA_DIST = example.hdf5
endif
