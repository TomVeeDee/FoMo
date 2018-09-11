lib_LTLIBRARIES = libFoMo.la
libFoMo_la_LDFLAGS = -shared -release @fomoversion@ -lboost_iostreams
libFoMo_ladir=$(includedir)
libFoMo_la_HEADERS=FoMo.h
libFoMo_la_SOURCES= \
	$(libFoMo_la_HEADERS) \
	FoMo-internal.h \
	../config.h \
	fomo-CGAL.cpp \
	fomo-CGAL2D.cpp \
	fomo-object.cpp \
	fomo-datacube.cpp \
	fomo-operations.cpp \
	fomo-goftcube.cpp \
	fomo-rendercube.cpp \
	fomo-CHIANTI.cpp \
	fomo-io.cpp sun_coronal.cpp \
	fomo-nearestneighbour.cpp \
	fomo-gpu-nearestneighbour.cpp \
	gpu-regulargrid-kernel.cpp \
	fomo-gpu-regulargrid.cpp \
	fomo-projection.cpp

