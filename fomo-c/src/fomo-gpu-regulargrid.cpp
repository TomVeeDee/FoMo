#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <set>
#include <chrono>

#include <CL/cl.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <functional>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, unsigned> value;

typedef struct __attribute__ ((packed)) RegularPoint {
	cl_float peak;
	cl_float fwhm;
	cl_float vx;
	cl_float vy;
	cl_float vz;
} RegularPoint;

typedef struct RegularGrid {
	float minx;
	float maxx;
	int x_pixel;
	float miny;
	float maxy;
	int y_pixel;
	float minz;
	float maxz;
	int z_pixel;
	RegularPoint* points; // Index: y*x_pixel*z_pixel + x*z_pixel + z
} RegularGrid;

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
const double pi=M_PI; //pi

static std::chrono::time_point<std::chrono::high_resolution_clock> time_now() {
	return std::chrono::high_resolution_clock::now();
}

RegularGrid* constructRegularGrid(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel)
{
//
// results is an array of at least dimension (x2-x1+1)*(y2-y1+1)*lambda_pixel and must be initialized to zero
// 
// determine contributions per pixel

// Allocates RegularGrid and RegularGrid->points on the heap, must be deleted by caller.
	
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	
	// Start timing
	std::chrono::time_point<std::chrono::high_resolution_clock> start = time_now();
	
	FoMo::tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();
	int dim=goftcube.readdim();
	
	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) std::cout << "Preparing coordinates ... " << std::flush;
	// Read the physical variables
	FoMo::tphysvar peakvec=goftcube.readvar(0);//Peak intensity 
	FoMo::tphysvar fwhmvec=goftcube.readvar(1);// line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);  
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);
	
	// initialisations for boost nearest neighbour
	point boostpoint, targetpoint;
	value boostpair;
	std::vector<value> input_values(ng),returned_values;
	box maxdistancebox;
	
	for (int i = 0; i < ng; i++) {
		// build r-tree from gridpoints, this part is not parallel
		boostpoint = point(grid[0][i], grid[1][i], (dim == 2 ? 0 : grid[2][i]));
		boostpair = std::make_pair(boostpoint, i);
		input_values.at(i) = boostpair;
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	if (commrank==0) std::cout << "Building R-tree... " << std::flush;
	// take an rtree with the quadratic packing algorithm, it takes (slightly) more time to build, but queries are faster for large renderings
	bgi::rtree< value, bgi::quadratic<16> > rtree(input_values.begin(), input_values.end());
	
	// compute the bounds of the input data points, so that we can equidistantly distribute the target pixels
	double minz = *(min_element(grid[2].begin(), grid[2].end()));
	double maxz = *(max_element(grid[2].begin(), grid[2].end()));
	double minx = *(min_element(grid[0].begin(), grid[0].end()));
	double maxx = *(max_element(grid[0].begin(), grid[0].end()));
	double miny = *(min_element(grid[1].begin(), grid[1].end()));
	double maxy = *(max_element(grid[1].begin(), grid[1].end()));
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	if (commrank==0) std::cout << "Building frame: " << std::flush;
	double x,y,z;
	
	// Initialize grid
	RegularPoint *points = new RegularPoint[y_pixel*x_pixel*z_pixel];
	RegularPoint* regular_point;
	
	// maxdistance is the furthest distance between a grid point and a simulation point at which the emission is interpolated
	// it is computed as the half diagonal of the rectangle around this ray, with the sides equal to the x and y distance between rays
	// i.e. it needs to be closer to this ray than to any other ray
	double maxdistance; 
	maxdistance = std::sqrt(std::pow((maxx-minx)/(x_pixel-1),2)+std::pow((maxy-miny)/(y_pixel-1),2))/2.;
	// However, it is better to just take the minimum of the pixel size in either direction, because it is then used in the maxdistancebox
	maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/2.;
	// If the viewing is along one of the axis, the previous value does not work very well, and the rendering almost always shows dark stripes: make the value 6 times larger!
	maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/.3;
	if ((maxx-minx)/std::pow(ng,1./3.)>maxdistance || (maxy-miny)/std::pow(ng,1./3.)>maxdistance) std::cout << std::endl << "Warning: maximum distance to interpolated point set to " << maxdistance << "Mm. If it is too small, you have too many interpolating rays and you will have dark stripes in the image plane. Reduce x-resolution or y-resolution." << std::endl;

	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	double deltaz = (maxz - minz);
	if (z_pixel != 1) deltaz /= (z_pixel - 1);
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(2) shared (rtree) private (x, y, z, returned_values, targetpoint, maxdistancebox, regular_point)
#endif
	for (int i = 0; i < y_pixel; i++)
		for (int j = 0; j < x_pixel; j++)
		{
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel - 1)*(maxx - minx) + minx;
			y = double(i)/(y_pixel - 1)*(maxy - miny) + miny;
						
			std::vector<double> p;
			
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k=0; k<z_pixel; k++) {
				z = double(k)*deltaz+minz;
				//p = {x, y, z};
				
				// initialise nearestindex to point -1
				int nearestindex=-1;
				
				// look for nearest point to targetpoint
				targetpoint=point(x, y, z);
				returned_values.clear();
				// the second condition ensures the point is not further away than 
				// - half the x-resolution in the x-direction
				// - half the y-resolution in the y-direction
				// - the maximum of both the previous numbers in the z-direction (sort of improvising a convex hull approach)
				//maxdistancebox=box(point(p.at(0)-(maxx-minx)/(x_pixel-1)/2.,p.at(1)-(maxy-miny)/(y_pixel-1)/2.,p.at(2)-maxdistance),point(p.at(0)+(maxx-minx)/(x_pixel-1)/2.,p.at(1)+(maxy-miny)/(y_pixel-1)/2.,p.at(2)+maxdistance));
				// it seems the expression above is the culprit for simulations with very stretched grids producing striped emissions, let's make the box of size maxdistance
				maxdistancebox=box(point(x-maxdistance,y-maxdistance,z-maxdistance),point(x+maxdistance,y+maxdistance,z+maxdistance));
				rtree.query(bgi::nearest(targetpoint, 1) && bgi::within(maxdistancebox), std::back_inserter(returned_values));
				regular_point = &(points[i*x_pixel*z_pixel + j*z_pixel + k]);
				if (returned_values.size() >= 1) {
					nearestindex = returned_values.at(0).second;
					regular_point->peak = 1e8*peakvec.at(nearestindex);
					regular_point->fwhm = fwhmvec.at(nearestindex);
					regular_point->vx = vx.at(nearestindex);
					regular_point->vy = vy.at(nearestindex);
					regular_point->vz = vz.at(nearestindex);
				} else {
					// All values must be initialized to deal with NaN and such
					regular_point->peak = 0;
					regularPoint.fwhm = 0;
					regularPoint.vx = 0;
					regularPoint.vy = 0;
					regularPoint.vz = 0;
				}

				// Print progress
				++show_progress;
			}
		}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	if (commrank==0) std::cout << "Constructing grid object ..." << std::flush;
	RegularGrid* regular_grid = new RegularGrid();
	regular_grid->minx = minx;
	regular_grid->maxx = maxx;
	regular_grid->x_pixel = x_pixel;
	regular_grid->miny = miny;
	regular_grid->maxy = maxy;
	regular_grid->y_pixel = y_pixel;
	regular_grid->minz = minz;
	regular_grid->maxz = maxz;
	regular_grid->z_pixel = z_pixel;
	regular_grid->points = points;
	
	if (commrank==0) std::cout << "Time spent in constructRegularGrid (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	return regular_grid;
	
}

namespace FoMo
{
	FoMo::RenderCube RenderWithGPURegularGrid(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		
		// Start timing
		std::chrono::time_point<std::chrono::high_resolution_clock> start = time_now();
		
		// Start of pre-processing
		
		// Initialize OpenCL
		if (commrank==0) std::cout << "Initializing OpenCL ... " << std::flush;
		// Find platforms
		std::vector<cl::Platform> cl_platforms;
		cl::Platform::get(&cl_platforms);
		if(cl_platforms.size() == 0) {
			std::cerr << "Error: No OpenCL platforms found!" << std::endl;
			exit(1);
		}
		// Make context
		cl_int err;
		cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(cl_platforms[0])(), 0};
		cl::Context cl_context(CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err);
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not create OpenCL context!" << std::endl;
			exit(1);
		}
		// Find devices
		std::vector<cl::Device> cl_devices;
		cl_devices = cl_context.getInfo<CL_CONTEXT_DEVICES>();
		if(cl_devices.size() == 0) {
			std::cerr << "Error: No OpenCL devices found!" << std::endl;
			exit(1);
		}
		// Load the program
		std::filebuf file;
		if(file.open("src/gpu-regulargrid.cl", std::ios_base::in | std::ios_base::binary) == NULL) {
			std::cerr << "Error: Could not load OpenCL program source!" << std::endl;
			exit(1);
		}
		std::streampos size = file.pubseekoff(0, std::ios_base::end);
		if(size < 0) {
			std::cerr << "Error: Could not load OpenCL program source!" << std::endl;
			exit(1);
		}
		std::vector<char> prog(size);
		file.pubseekoff(0, std::ios_base::beg);
		std::streamsize read = file.sgetn(prog.data(), prog.size());
		if((size_t) read != prog.size()) {
			std::cerr << "Error: Could not load OpenCL program source!" << std::endl;
			exit(1);
		}
		prog.push_back('\0');
		// Compile the program
		cl::Program::Sources cl_program_source(1, std::make_pair(prog.data(), prog.size()));
		cl::Program cl_program(cl_context, cl_program_source);
		std::string build_options = "-cl-nv-verbose -D MAX_LAMBDA_PIXEL=" + std::to_string(lambda_pixel) + " -D MAX_DEPTH=" + std::to_string(int(floor(log2(ng))));
		err = cl_program.build(cl_devices, build_options.c_str());
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not compile OpenCL program!" << std::endl;
			std::cerr << "OpenCL build log:\n" << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl_devices[0]) << std::endl;
			exit(1);
		}
		// Allocate buffers
		int input_amount = ng;
		int pixels = x_pixel*y_pixel;
		int output_amount = std::min(chunk_size, pixels)*lambda_pixel;
		cl::Buffer cl_buffer_coords(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, 3*sizeof(float)*(input_amount + 1));
		cl::Buffer cl_buffer_nodes(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(uint8_t)*(input_amount + 1)); // Tree metadata: index of axis of split for non-leaves, 3 for leaves
		cl::Buffer cl_buffer_data_in(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, data_in_per_point*sizeof(float)*(input_amount + 1)); // Peak, fwhm, losvel
		cl::Buffer cl_buffer_parameters0(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(parameters_struct)); // Various constant parameters      
		cl::Buffer cl_buffer_parameters1(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(parameters_struct)); // Various constant parameters
		cl::Buffer cl_buffer_parameters[] = {cl_buffer_parameters0, cl_buffer_parameters1};
		cl::Buffer cl_buffer_data_out0(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, data_out_per_point*sizeof(float)*output_amount);
		cl::Buffer cl_buffer_data_out1(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, data_out_per_point*sizeof(float)*output_amount);
		cl::Buffer cl_buffer_data_out[] = {cl_buffer_data_out0, cl_buffer_data_out1};
		//cl::Buffer cl_buffer_debug(cl_context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, 20*sizeof(double));
		// Create queue
		cl::CommandQueue cl_queue0(cl_context, cl_devices[0], 0, &err); // Used for first kernel and general stuff
		cl::CommandQueue cl_queue1(cl_context, cl_devices[0], 0, &err); // Used for second kernel
		cl::CommandQueue queues[] = {cl_queue0, cl_queue1};
		assert(err == CL_SUCCESS);
		// Map buffers
		float *coords = (float*) queues[0].enqueueMapBuffer(cl_buffer_coords, CL_FALSE, CL_MAP_WRITE, 0, 3*sizeof(float)*(input_amount + 1));
		uint8_t *nodes = (uint8_t*) queues[0].enqueueMapBuffer(cl_buffer_nodes, CL_FALSE, CL_MAP_WRITE, 0, sizeof(uint8_t)*(input_amount + 1));
		float *data_in = (float*) queues[0].enqueueMapBuffer(cl_buffer_data_in, CL_FALSE, CL_MAP_WRITE, 0, data_in_per_point*sizeof(float)*(input_amount + 1));
		parameters_struct *parameters0 = (parameters_struct*) queues[0].enqueueMapBuffer(cl_buffer_parameters[0], CL_FALSE, CL_MAP_WRITE, 0, sizeof(parameters_struct));
		parameters_struct *parameters1 = (parameters_struct*) queues[1].enqueueMapBuffer(cl_buffer_parameters[0], CL_FALSE, CL_MAP_WRITE, 0, sizeof(parameters_struct));
		parameters_struct *parameters[] = {parameters0, parameters1};
		float *data_out0 = (float*) queues[0].enqueueMapBuffer(cl_buffer_data_out[0], CL_FALSE, CL_MAP_READ, 0, data_out_per_point*sizeof(float)*output_amount);
		float *data_out1 = (float*) queues[1].enqueueMapBuffer(cl_buffer_data_out[0], CL_FALSE, CL_MAP_READ, 0, data_out_per_point*sizeof(float)*output_amount);
		float *data_out[] = {data_out0, data_out1};
		// Create kernels
		cl::Kernel cl_kernel0(cl_program, "calculate_ray", &err);
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not create kernel: " << err << std::endl;
			exit(1);
		}
		cl::Kernel cl_kernel1(cl_program, "calculate_ray", &err);
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not create kernel: " << err << std::endl;
			exit(1);
		}
		cl::Kernel kernels[] = {cl_kernel0, cl_kernel1};
		for(int i = 0; i < 2; i++) {
			kernels[i].setArg(0, cl_buffer_coords);
			kernels[i].setArg(1, cl_buffer_nodes);
			kernels[i].setArg(2, cl_buffer_data_in);
			kernels[i].setArg(3, cl_buffer_parameters[i]);
			kernels[i].setArg(4, cl_buffer_data_out[i]);
			//kernels[i].setArg(5, cl_buffer_debug);
		}
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		RegularGrid* regular_grid = constructRegularGrid(goftcube, x_pixel, y_pixel, z_pixel);
		std::cout << (regular_grid->points)[0].vx << std::endl;
		
		FoMo::RenderCube rendercube(goftcube);
		/*for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit) {
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit)
			{
				rendercube=nearestneighbourinterpolation(goftcube,*lit,*bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
				rendercube.setangles(*lit,*bit);
				std::stringstream ss;
				// if outfile is "", then this should not be executed.
				ss << outfile;
				ss << "l";
				ss << std::setfill('0') << std::setw(3) << std::round(*lit/pi*180.);
				ss << "b";
				ss << std::setfill('0') << std::setw(3) << std::round(*bit/pi*180.);
				ss << ".txt";
				rendercube.writegoftcube(ss.str());
				ss.str("");
			}
		}*/
		return rendercube;
	}
}

