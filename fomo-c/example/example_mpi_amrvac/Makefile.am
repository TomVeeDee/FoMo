  noinst_PROGRAMS = example_mpi_amrvac render_all_datfiles
  LDADD = -L$(top_builddir)/src/.libs/ -lFoMo 
  AM_CPPFLAGS = -I$(top_srcdir)/src
  example_mpi_amrvac_SOURCES = example_mpi_amrvac.cpp read_amrvac_files.cpp
  render_all_datfiles_LDADD=-L$(top_builddir)/src/.libs/ -lFoMo -lboost_system -lboost_filesystem -lboost_program_options -lboost_regex
  render_all_datfiles_SOURCES=render_all_datfiles.cpp read_amrvac_files.cpp
