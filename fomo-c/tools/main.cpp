#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <boost/call_traits.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#include <boost/call_traits.hpp>

double xl, xr, yl, yr, zl, zr;
int dimr, dimx, dimy, dimz;
std::string mystr, line;
double x,y,z;
double i,j,k;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> 					 Polyhedron;
typedef K::Point_3 				      				 Point_3, Point;
typedef std::vector<Point_3>	 					 Points;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;

void writegrid(double xl, double xr, int dimx, double yl, double yr, int dimy, double zl, double zr, int dimz) {

	std::ofstream outfile ("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/grids/pointgridcompletecyl.dat");
	std::cout << "Writing grid..." << std::endl;
	
	for (int i=0; i< dimx; i++) {
		for (int j=0; j<dimy; j++) {
			for (int k=0; k<dimz; k++) {
				outfile << xl + (xr - xl) * i/(dimx-1) << "\t" << yl + (yr - yl) *j/(dimy-1) << "\t" << zl + (zr - zl) *k/(dimz-1) << std::endl;
			}
		}
	}
	outfile.close();
}

double Convertstringtofloat (std::string a) {
double fl;
std::stringstream convert(a); // stringstream used for the conversion initialized with the contents of Text
	if ( !(convert >> fl) ) return 0;
	else  return fl;
}

Polyhedron build_polyhedron(int t)

// Build a polyhedron, which is the convex hull of all points that lie
// in the loop interior at the given time step.

{
	Point_3 gridpoint;
	Points initial, withgridpoint;	
	Polyhedron poly, polywithgridpoint;			// to store the points of the convex hull (vertices, edges, faces and adjacency relations)
	std::stringstream ss;
	ss << "/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoordright";
	ss << std::setfill('0') << std::setw(3) << t;
	ss << ".dat";
	std::string filename= ss.str();
	std::ifstream myfile(filename.c_str());

	std::cout << "Building polyhedron " << t << std::endl;

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			std::getline (myfile,line);
		  x = Convertstringtofloat(line);
			std::getline (myfile,line);
		  y = Convertstringtofloat(line);
			std::getline (myfile,line);
		  z = Convertstringtofloat(line);	
			initial.push_back(Point_3(x,y,z));
			}
	}

	else std::cout << "Unable to open file";

	myfile.close();
	CGAL::convex_hull_3(initial.begin(), initial.end(), poly);

	std::cout << "There are " << poly.size_of_vertices() << " vertices on the convex hull." << std::endl;
	std::stringstream sp;
	sp << "/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/hullpolyhedra/secondrun/largeba5104completecyl/hulltr";
	sp << std::setfill('0') << std::setw(3) << t;
	sp << ".dat";
	std::string filenamep=sp.str();
	std::ofstream writepoly(filenamep.c_str());

	if (writepoly.is_open()) {
 		// Write polyhedron in Object File Format (OFF).
		CGAL::set_ascii_mode( writepoly);
		writepoly << "OFF" << std::endl << poly.size_of_vertices() << ' '<< poly.size_of_facets() << " 0" << std::endl;
		std::copy( poly.points_begin(), poly.points_end(), std::ostream_iterator<Point_3>( writepoly, "\n"));
		for ( Facet_iterator i = poly.facets_begin(); i != poly.facets_end(); ++i) {
				Halfedge_facet_circulator j = i->facet_begin();
				// Facets in polyhedral surfaces are at least triangles.
				CGAL_assertion( CGAL::circulator_size(j) >= 3);
				writepoly << CGAL::circulator_size(j) << ' ';
				do {
						writepoly << ' ' << std::distance(poly.vertices_begin(), j->vertex());
						} while ( ++j != i->facet_begin());
						writepoly << std::endl;
		}
	}
	else {std::cout << "Unable to open file";}
	return (poly);
}

void polyhedron_prog_inside_test(Polyhedron P, int t)

// Fast algorithm to decide for each point of the grid whether it lies inside or
// outside the loop, which is constructed as the convex hull of the interior
// points. It returns a locfile with '1' if the points lies within the convex hull,
// and '0' otherwise.

{
  std::ifstream points_file("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/grids/pointgridcompletecyl.dat");
	std::stringstream ss;
	ss << "/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/locfiles/largeba5103completecyl/locfilet";
	ss << std::setfill('0') << std::setw(3) << t;
	ss << ".dat";
	std::string filename=ss.str();
	std::ofstream locfile(filename.c_str());

	std::cout << "Checking grid points for convex hull at time " << t << std::endl;

  std::vector<Point> points;
  std::copy(std::istream_iterator<Point>(points_file),
            std::istream_iterator<Point>(),
            std::back_inserter(points));

  int nb_points=points.size();
  
  std::vector<bool> ray_res(nb_points);
  std::vector<bool> grid_res(nb_points);
  
  CGAL::Timer timer;
  timer.start();
  CGAL::Point_inside_polyhedron_3<Polyhedron,K> inside_with_ray(P);
  timer.stop();
  std::cerr <<"Using ray"<< std::endl;
  std::cerr << "  Preprocessing took " << timer.time() << " sec." << std::endl;
  timer.reset();  
  int n_inside = 0;  
  
  timer.start();
  for(int k=0;k<nb_points;++k){
    ray_res[k]=inside_with_ray(points[k]);
    if(ray_res[k]){
      ++n_inside;
    }
		locfile << ray_res[k] << std::endl;
  }
  timer.stop();
  std::cerr << "  " << n_inside << " points inside " << std::endl;
  std::cerr << "  " << points.size() - n_inside << " points outside "  << std::endl;
  std::cerr << " Queries took " << timer.time() << " sec." << std::endl; 

	locfile.close();
}

int main() {

int dimt = 1;
std::ifstream myfile ("/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/griddata.dat");

// Read the grid parameters from IDL to write out Eulerian grid points.

std::getline (myfile,mystr);
std::stringstream(mystr) >> dimr;
std::getline (myfile,mystr);
std::stringstream(mystr) >> dimz;
std::getline (myfile,mystr);
std::stringstream(mystr) >> xr;
std::getline (myfile,mystr);
std::stringstream(mystr) >> zr;
xl = -xr;
dimx = dimr;
yl = -xr;		// (Must be adapted whether we take half of or the complete loop)
yr = xr;
dimy = dimr;
zl = 0;

// writegrid(xl, xr, dimx, yl, yr, dimy, zl, zr, dimz);

for(int t = 0; t < dimt; t++){
	Polyhedron poly = build_polyhedron(t);
//	polyhedron_prog_inside_test(poly, t);
	}
return 0;
}
