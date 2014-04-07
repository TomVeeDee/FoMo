#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>

#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


#include <CGAL/Timer.h>
#include <boost/call_traits.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;


void polyhedron_prog_inside_test()
{
  std::ifstream points_file("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/grids/pointgrid.dat");
	std::ofstream locfile("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/locfiles/locfilet0.dat");

  std::vector<Point> points;
  std::copy(std::istream_iterator<Point>(points_file),
            std::istream_iterator<Point>(),
            std::back_inserter(points));


  int nb_points=points.size();
  
  std::vector<bool> ray_res(nb_points);
  std::vector<bool> grid_res(nb_points);
  
  Polyhedron polyhedron;
  std::ifstream polyhedron_file("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/hullpolyhedra/hullt39.off");
  polyhedron_file >> polyhedron;
  std::cerr << "|V| = " << polyhedron.size_of_vertices() << std::endl;
  
  CGAL::Timer timer;
  timer.start();
  CGAL::Point_inside_polyhedron_3<Polyhedron,K> inside_with_ray(polyhedron);
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
	
  
  return 0;
}
