#include "../config.h"
#include <iostream>
#include <complex>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <boost/call_traits.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

extern void build_polyhedron();
extern void writegrid(double xl, double xr, int dimx, double yl, double yr, int dimy, double zl, double zr, int dimz);
extern void polyhedron_prog_inside_test();
