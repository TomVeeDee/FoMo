// CGAL and vector environment
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <vector>
// Input and output environment
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 				      				 Point_3;
typedef CGAL::Polyhedron_3<K>            Polyhedron_3;
typedef std::vector<Point_3>	 					 Points;
typedef Polyhedron_3::Facet_iterator Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator Halfedge_facet_circulator;

std::string line;
double x,y,z;
int counter = 1;
double i,j,k;

double Convertstringtofloat (std::string a) {

double fl;
std::stringstream convert(a); // stringstream used for the conversion initialized with the contents of Text
	if ( !(convert >> fl) ) return 0;
	else  return fl;
}


void build_polyhedron()
{
	Point_3 gridpoint;
	Points initial, withgridpoint;
	Polyhedron_3 poly, polywithgridpoint;			// to store the points of the convex hull in (vertices, edges, faces and adjacency relations)
	std::ifstream myfile ("/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/advectedcoord000.dat");

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
			counter++;	
			initial.push_back(Point_3(x,y,z));

			if(counter==100 ||	counter==1000 || counter==10000 || counter == 100000 || counter == 200000 || counter == 500000) {
				std::cout << counter << std::endl;
			}
		}
	}

	else std::cout << "Unable to open file";		// Important to catch wrong input while programming

	myfile.close();
	CGAL::convex_hull_3(initial.begin(), initial.end(), poly);

	std::cout << "There are " << poly.size_of_vertices() << " vertices on the convex hull." << std::endl;
	std::ofstream writepoly ("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/hullpolyhedra/hullt00.off");

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
	else std::cout << "Unable to open file";
	
	return 0;
}
