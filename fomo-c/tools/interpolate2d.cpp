#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>               Traits;
typedef K::FT                                         Coord_type;
typedef K::Point_2                                    Point, Point_2;


double Convertstringtofloat (std::string a) {

double fl;
std::stringstream convert(a); // stringstream used for the conversion initialized with the contents of Text
	if ( !(convert >> fl) ) return 0;
	else  return fl;
}


void interpolate(int t, int z) {

std::string line;
double x,y,rhoin,tein,brin,btin,bzin;

Delaunay_triangulation T;
std::map<Point, Coord_type, K::Less_xy_2> rho, te, br, bt, bz;
typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2 > >
                                            Value_access;

	// Files for input and output

	CGAL::Timer timerrr;
	timerrr.start();

	std::stringstream sin, sout;
		sin << "/users/cpa/sgijsen/idl_general/test_advectioncube/eigft";
		sin << std::setfill('0') << std::setw(3) << t;
		sin << "z";
		sin << std::setfill('0') << std::setw(3) << z;
		sin << ".dat";

		sout << "/users/cpa/sgijsen/idl_general/test_advectioncube/eigftcube";
		sout << std::setfill('0') << std::setw(3) << t;
		sout << "z";
		sout << std::setfill('0') << std::setw(3) << z;
		sout << ".dat";

	std::string filenamein = sin.str();
	std::string filenameout = sout.str();
	std::ifstream myfile(filenamein.c_str());
	std::ofstream outfile(filenameout.c_str());

	// Read in data from input file; triangulate points and make pairs

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			std::getline (myfile,line);
				x = Convertstringtofloat(line);
			std::getline (myfile,line);
				y = Convertstringtofloat(line);
			//std::getline (myfile,line);
			//  vxin = Convertstringtofloat(line);
			//std::getline (myfile,line);
			//  vyin = Convertstringtofloat(line);
			std::getline (myfile,line);
			  rhoin = Convertstringtofloat(line);
			std::getline (myfile,line);
			  tein = Convertstringtofloat(line);			
			std::getline (myfile,line);
			  brin = Convertstringtofloat(line);
			std::getline (myfile,line);
			  btin = Convertstringtofloat(line);
			std::getline (myfile,line);
			  bzin = Convertstringtofloat(line);
      K::Point_2 p(x,y);
      T.insert(p);

				//vx.insert(std::make_pair(p,vxin));
				//vy.insert(std::make_pair(p,vyin));
				rho.insert(std::make_pair(p,rhoin));
				te.insert(std::make_pair(p,tein));
				br.insert(std::make_pair(p,brin));
				bt.insert(std::make_pair(p,btin));
				bz.insert(std::make_pair(p,bzin));
		}
	}

	else std::cout << "Unable to open file";		// Catch wrong input while programming   (16/07/14)


  // Linear interpolation using natural neighbor weights
	// See CGAL documentation on 2D interpolation

	double dimx=180, dimy=180, xl=-2*6.6, xr=2*6.6, yl=-2*6.6, yr=2*6.6;

	for (int i=0; i< dimx; i++) {
		for (int j=0; j<dimy; j++) {
  		K::Point_2 p(xl + (xr - xl) * i/(dimx-1), yl + (yr - yl) * j/(dimy-1));
  		std::vector< std::pair< Point, Coord_type > > coords;
  		Coord_type norm =
    		CGAL::natural_neighbor_coordinates_2
    		(T, p,std::back_inserter(coords)).second;

  		//Coord_type vxres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
			//			               norm,
			//				       Value_access(vx));
  		//Coord_type vyres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
			//			               norm,
			//				       Value_access(vy));
  		Coord_type rhores =  CGAL::linear_interpolation(coords.begin(), coords.end(),
						               norm,
							       Value_access(rho));
  		Coord_type teres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
						               norm,
							       Value_access(te));
  		Coord_type brres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
						               norm,
							       Value_access(br));
  		Coord_type btres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
						               norm,
							       Value_access(bt));
  		Coord_type bzres =  CGAL::linear_interpolation(coords.begin(), coords.end(),
						               norm,
							       Value_access(bz));

  		 outfile << xl + (xr - xl) * i/(dimx-1) << "\t" << yl + (yr - yl) * j/(dimy-1) << "\t" << rhores << "\t" << teres << "\t" << brres << "\t" << btres << "\t" << bzres << std::endl;

		}
	}
	outfile.close();
	// std::cout << "Intermediate lap " << z << ": " << timerrr.time() << " sec." << std::endl; 
}


int main() {

	CGAL::Timer timer;
	timer.start();
	int dimz=210, dimt=16;
	for (int t=0; t<dimt; t++){
		 for (int z=0; z<dimz; z++) {
			interpolate(t,z);
		}	
std::cout << " Queries took " << timer.time() << " sec." << std::endl;
	}

return 0;

}

