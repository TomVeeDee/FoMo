#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <cstdlib>
#include <sstream>
#include <fstream>
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

#define GPU_REGULAR_GRID_DEBUG 0
#define GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE 200

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, unsigned> value;

typedef struct __attribute__ ((packed)) Parameters {
	cl_float rxx;
	cl_float rxy;
	cl_float rxz;
	cl_float ryx;
	cl_float ryy;
	cl_float ryz;
	cl_float rzx;
	cl_float rzy;
	cl_float rzz;
	cl_float x_offset;
	cl_float y_offset;
} Parameters;

typedef struct RegularGrid {
	// Offset indicates the coordinates of the new origin (the middle of the grid) in the old coordinate system
	float x_offset;
	float minx;
	float y_offset;
	float miny;
	float z_offset;
	float minz;
	cl_float8* points; // Index: y*x_pixel*z_pixel + x*z_pixel + z. Vector format: {peak, fwhm, vx, vy, vz, undefined, undefined, undefined}
} RegularGrid;

const double pi=M_PI; //pi

inline static std::chrono::time_point<std::chrono::high_resolution_clock> time_now() {
	return std::chrono::high_resolution_clock::now();
}

inline int sgn(float val) {
    return (0 < val) - (val < 0);
}

void enqueueKernel2(cl::CommandQueue queue, cl::Kernel kernel, int offset, int size) {
	int err = queue.enqueueNDRangeKernel(kernel, cl::NDRange(offset), cl::NDRange(size), cl::NullRange);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue the kernel for execution: " << err << std::endl;
        exit(1);
    }
}

void enqueueRead2(cl::CommandQueue queue, cl::Buffer cl_buffer_data_out, int size, float *data_out) {
	int err = queue.enqueueReadBuffer(cl_buffer_data_out, CL_TRUE, 0, size, data_out);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue buffer read: " << err << std::endl;
        exit(1);
    }
}

inline void rotateAroundZ(float* in, float angle, float* out) {
	
	// in/out are 3-component vectors
	// angle is in radians
	
	float cosa = cos(angle);
	float sina = sin(angle);
	out[0] = cosa*in[0] - sina*in[1];
	out[1] = sina*in[0] + cosa*in[1];
	out[2] = in[2];
	
}

inline void rotateAroundY(float* in, float angle, float* out) {
	
	// in/out are 3-component vectors
	// angle is in radians
	
	float cosa = cos(angle);
	float sina = sin(angle);
	out[0] = cosa*in[0] + sina*in[2];
	out[1] = in[1];
	out[2] = -sina*in[0] + cosa*in[2];
	
}

RegularGrid* constructRegularGrid(FoMo::GoftCube goftcube, cl_float8* points, const int x_pixel, const int y_pixel, const int z_pixel)
{
//
// results is an array of at least dimension (x2-x1+1)*(y2-y1+1)*lambda_pixel and must be initialized to zero
// 
// determine contributions per pixel

// Allocates RegularGrid on the heap, must be deleted by caller. The points array must be allocated by the caller already and a pointer to the start of it is passed as an argument.
	
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
	
	if (commrank==0) std::cout << "Building regular grid: " << std::flush;
	double x,y,z;
	
	// maxdistance is the furthest distance between a grid point and a simulation point at which the emission is interpolated
	// it is computed as the half diagonal of the rectangle around this ray, with the sides equal to the x and y distance between rays
	// i.e. it needs to be closer to this ray than to any other ray
	double maxdistance;
	//maxdistance = std::sqrt(std::pow((maxx-minx)/(x_pixel-1),2)+std::pow((maxy-miny)/(y_pixel-1),2))/2.;
	// However, it is better to just take the minimum of the pixel size in either direction, because it is then used in the maxdistancebox
	//maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/2.;
	// If the viewing is along one of the axis, the previous value does not work very well, and the rendering almost always shows dark stripes: make the value 6 times larger!
	maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/.3;
	if ((maxx-minx)/std::pow(ng,1./3.)>maxdistance || (maxy-miny)/std::pow(ng,1./3.)>maxdistance) std::cout << std::endl << "Warning: maximum distance to interpolated point set to " << maxdistance << "Mm. If it is too small, you have too many interpolating rays and you will have dark stripes in the image plane. Reduce x-resolution or y-resolution." << std::endl;

	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	double deltaz = (maxz - minz);
	if (z_pixel != 1) deltaz /= (z_pixel - 1);
	
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic) collapse(2) shared (rtree, points) private (x, y, z, returned_values, targetpoint, maxdistancebox/*, regular_point*/)
	#endif
	for (int i = 0; i < y_pixel; i++) {
		for (int j = 0; j < x_pixel; j++)
		{
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel - 1)*(maxx - minx) + minx;
			y = double(i)/(y_pixel - 1)*(maxy - miny) + miny;
						
			std::vector<double> p;
			
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k = 0; k < z_pixel; k++) {
				//std::cout << i << " " << j << " " << k << std::endl << std::flush;
				z = double(k)*deltaz+minz;
				
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
				if (returned_values.size() >= 1) {
					nearestindex = returned_values.at(0).second;
					points[i*x_pixel*z_pixel + j*z_pixel + k] = {peakvec.at(nearestindex), fwhmvec.at(nearestindex), vx.at(nearestindex), vy.at(nearestindex), vz.at(nearestindex), 0, 0, 0};
				} else {
					// All values must be initialized to deal with NaN and such
					points[i*x_pixel*z_pixel + j*z_pixel + k] = {0, 1, 0, 0, 0, 0, 0, 0};
				}

				// Print progress
				++show_progress;
			}
		}
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	std::cout << minx << " " << maxx << " " << miny << " " << maxy << " " << minz << " " << maxz << " " << std::endl;
	
	if (commrank==0) std::cout << "Constructing grid object ..." << std::flush;
	RegularGrid* regular_grid = new RegularGrid();
	regular_grid->x_offset = (minx + maxx)/2.0;
	regular_grid->minx = (minx - maxx)/2.0*float(x_pixel)/float(x_pixel - 1);
	regular_grid->y_offset = (miny + maxy)/2.0;
	regular_grid->miny = (miny - maxy)/2.0*float(y_pixel)/float(y_pixel - 1);
	regular_grid->z_offset = (minz + maxz)/2.0;
	regular_grid->minz = (dim < 3 || z_pixel == 1 ? -0.5 : (minz - maxz)/2.0*float(z_pixel)/float(z_pixel - 1)); // If dimension is smaller than 3, set Z-grid-size to 1 Mm.
	regular_grid->points = points;
	if (dim < 3 || z_pixel == 1) {
		std::cout << "Assuming that this is a 2D simulations: setting thickness of simulation to 1 Mm." << std::endl << std::flush;
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	if (commrank==0) std::cout << "Time spent in constructRegularGrid (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	return regular_grid;
	
}

namespace FoMo
{
	FoMo::RenderCube RenderWithGPURegularGrid(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
		std::vector<double> lvec, std::vector<double> bvec, std::string outfile) {
	
		int commrank;
		#ifdef HAVEMPI
			MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
		#else
			commrank = 0;
		#endif
		
		// Start timing
		std::chrono::time_point<std::chrono::high_resolution_clock> start = time_now();
		
		// Constants
		const int data_out_per_point = 1;
		const int chunk_size = 1024*2; // Amount of jobs submitted to the GPU simultaneously
		
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
		cl::Program::Sources cl_program_source(1, std::make_pair(prog.data(), prog.size()));
		cl::Program cl_program(cl_context, cl_program_source);
		// Allocate buffers
		int input_size = sizeof(cl_float8)*x_pixel*y_pixel*z_pixel;
		int pixels = x_pixel*y_pixel;
		int output_amount = std::min(chunk_size, pixels)*lambda_pixel;
		cl::Buffer cl_buffer_points(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, input_size);
		cl::Buffer cl_buffer_lambdaval(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(float)*lambda_pixel);
		cl::Buffer cl_buffer_parameters(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(Parameters)); // Various constant parameters
		cl::Buffer cl_buffer_data_out0(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, data_out_per_point*sizeof(float)*output_amount);
		cl::Buffer cl_buffer_data_out1(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, data_out_per_point*sizeof(float)*output_amount);
		cl::Buffer cl_buffer_data_out[] = {cl_buffer_data_out0, cl_buffer_data_out1};
		#if (GPU_REGULAR_GRID_DEBUG == 1)
			cl::Buffer cl_buffer_debug(cl_context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(float));
		#endif
		// Create queue
		cl::CommandQueue cl_queue0(cl_context, cl_devices[0], 0, &err); // Used for first kernel and general stuff
		cl::CommandQueue cl_queue1(cl_context, cl_devices[0], 0, &err); // Used for second kernel
		cl::CommandQueue queues[] = {cl_queue0, cl_queue1};
		assert(err == CL_SUCCESS);
		// Map buffers
		cl_float8 *points = (cl_float8*) queues[0].enqueueMapBuffer(cl_buffer_points, CL_FALSE, CL_MAP_WRITE, 0, input_size);
		float *lambdaval = (float*) queues[0].enqueueMapBuffer(cl_buffer_lambdaval, CL_FALSE, CL_MAP_WRITE, 0, sizeof(float)*lambda_pixel);
		Parameters *parameters = (Parameters*) queues[0].enqueueMapBuffer(cl_buffer_parameters, CL_FALSE, CL_MAP_WRITE, 0, sizeof(Parameters));
		float *data_out0 = (float*) queues[0].enqueueMapBuffer(cl_buffer_data_out[0], CL_FALSE, CL_MAP_READ, 0, data_out_per_point*sizeof(float)*output_amount);
		float *data_out1 = (float*) queues[1].enqueueMapBuffer(cl_buffer_data_out[1], CL_FALSE, CL_MAP_READ, 0, data_out_per_point*sizeof(float)*output_amount);
		float *data_out[] = {data_out0, data_out1};
		#if (GPU_REGULAR_GRID_DEBUG == 1)
			float *debug_buffer = (float*) queues[0].enqueueMapBuffer(cl_buffer_debug, CL_FALSE, CL_MAP_READ, 0, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(float));
		#endif
		queues[0].finish();
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		// Construct regular grid
		if (commrank==0) std::cout << "Constructing regular grid ... " << std::flush;
		RegularGrid *regular_grid = constructRegularGrid(goftcube, points, x_pixel, y_pixel, z_pixel);
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		// Compile program and create kernels
		if (commrank==0) std::cout << "Compiling OpenCL program and creating kernels ... " << std::flush;
		// Build program
		float g[] = {-2*regular_grid->minx/x_pixel, -2*regular_grid->miny/y_pixel, -2*regular_grid->minz/z_pixel};
		float pixel_width = g[0];
		float pixel_height = g[1];
		float lambda0 = float(goftcube.readlambda0());
		float ox = x_pixel/2.0;
		float oy = y_pixel/2.0;
		std::string build_options = "-cl-nv-verbose -D GPU_REGULAR_GRID_DEBUG=" + std::to_string(GPU_REGULAR_GRID_DEBUG) + " -D X_PIXEL=" + std::to_string(x_pixel)
			+ " -D PIXEL_WIDTH=" + std::to_string(pixel_width) + " -D PIXEL_HEIGHT=" + std::to_string(pixel_height)
			+ " -D LAMBDA_PIXEL=" + std::to_string(lambda_pixel) + " -D LAMBDA0=" + std::to_string(lambda0)
			+ " -D MINX=" + std::to_string(regular_grid->minx) + " -D MAXX=" + std::to_string(-regular_grid->minx)
			+ " -D MINY=" + std::to_string(regular_grid->miny) + " -D MAXY=" + std::to_string(-regular_grid->miny)
			+ " -D MINZ=" + std::to_string(regular_grid->minz) + " -D MAXZ=" + std::to_string(-regular_grid->minz)
			+ " -D GX=" + std::to_string(g[0]) + " -D GSX=" + std::to_string(x_pixel)
			+ " -D GY=" + std::to_string(g[1]) + " -D GSY=" + std::to_string(y_pixel)
			+ " -D GZ=" + std::to_string(g[2]) + " -D GSZ=" + std::to_string(z_pixel)
			+ " -D OX=" + std::to_string(ox) + " -D OY=" + std::to_string(oy);
		err = cl_program.build(cl_devices, build_options.c_str());
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not compile OpenCL program!" << std::endl;
			std::cerr << "OpenCL build log:\n" << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl_devices[0]) << std::endl;
			exit(1);
		}
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
		// Set kernel arguments
		for(int i = 0; i < 2; i++) {
			kernels[i].setArg(0, cl_buffer_points);
			kernels[i].setArg(1, cl_buffer_lambdaval);
			kernels[i].setArg(2, cl_buffer_parameters);
			kernels[i].setArg(3, cl_buffer_data_out[i]);
			#if (GPU_REGULAR_GRID_DEBUG == 1)
				kernels[i].setArg(4, cl_buffer_debug);
			#endif
		}
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		// Other initialization
		if (commrank==0) std::cout << "Other initialization ... " << std::flush;
		FoMo::RenderCube rendercube(goftcube);
		float lambda_width_in_A = lambda_width*lambda0/GSL_CONST_MKSA_SPEED_OF_LIGHT;
		for(int i = 0; i < lambda_pixel; i++)
			lambdaval[i] = float(i)/(lambda_pixel - 1)*lambda_width_in_A - lambda_width_in_A/2.0;
		queues[0].enqueueWriteBuffer(cl_buffer_points, CL_FALSE, 0, input_size, points);
		queues[0].enqueueWriteBuffer(cl_buffer_lambdaval, CL_FALSE, 0, sizeof(float)*lambda_pixel, lambdaval);
		queues[0].finish();
		FoMo::tgrid newgrid;
		FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel), yvec(x_pixel*y_pixel*lambda_pixel);
		newgrid.push_back(xvec);
		newgrid.push_back(yvec);
		if (lambda_pixel > 1) {
			FoMo::tcoord lambdavec(x_pixel*y_pixel*lambda_pixel);
			newgrid.push_back(lambdavec);
		}
		FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel, 0);
		FoMo::tvars newdata;
		newdata.push_back(intens);
		float xs[x_pixel];
		for(int x = 0; x < x_pixel; x++)
			xs[x] = (x + 0.5 - ox)*pixel_width + parameters->x_offset;
		float ys[y_pixel];
		for(int y = 0; y < y_pixel; y++)
			ys[y] = (y + 0.5 - oy)*pixel_height + parameters->y_offset;
		float lambdas[lambda_pixel];
		for(int l = 0; l < lambda_pixel; l++)
			lambdas[l] = lambdaval[l] + lambda0;
		float* intensity = new float[y_pixel*x_pixel*lambda_pixel];
		float temp[3];
		float inx[] = {1, 0, 0};
		float iny[] = {0, 1, 0};
		float inz[] = {0, 0, 1};
		float global_offset_vector[] = {regular_grid->x_offset, regular_grid->y_offset, regular_grid->z_offset};
		float rx[3];
		float ry[3];
		float rz[3];
		float local_offset_vector[3];
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		// End of pre-processing
		
		// Render frames
		int frame_index = 0;
		double frame_time = 0;
		for (std::vector<double>::iterator lit = lvec.begin(); lit != lvec.end(); ++lit) {
			for (std::vector<double>::iterator bit = bvec.begin(); bit != bvec.end(); ++bit) {
				
				double l = *lit;
				double b = *bit;
				
				for(int frame_counter = 0; frame_counter < 5; frame_counter++) {
				
				if (commrank==0) std::cout << "Rendering frame " << frame_index++ << " at " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
				
				// Calculate frame parameters
				if (commrank==0) std::cout << "Calculating frame parameters ... " << std::flush;
				rotateAroundY(inx, b, temp);
				rotateAroundZ(temp, -l, rx);
				rotateAroundZ(iny, -l, ry);
				rotateAroundY(inz, b, temp);
				rotateAroundZ(temp, -l, rz);
				rotateAroundZ(global_offset_vector, l, temp);
				rotateAroundY(temp, -b, local_offset_vector);
				parameters->rxx = rx[0]; parameters->rxy = rx[1]; parameters->rxz = rx[2];
				parameters->ryx = ry[0]; parameters->ryy = ry[1]; parameters->ryz = ry[2];
				parameters->rzx = rz[0]; parameters->rzy = rz[1]; parameters->rzz = rz[2];
				parameters->x_offset = local_offset_vector[0]; parameters->y_offset = local_offset_vector[1];
				queues[0].enqueueWriteBuffer(cl_buffer_parameters, CL_FALSE, 0, sizeof(Parameters), parameters);
				queues[0].finish();
				if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
				
				if (commrank==0) std::cout << "Building frame and extracting output ... " << std::flush;
				std::chrono::time_point<std::chrono::high_resolution_clock> frame_start = time_now();
				
				// Ping-pong in between two kernels until work is finished
				int index = 0;
				enqueueKernel2(queues[index], kernels[index], 0, std::min(chunk_size, pixels));
				for(int offset = chunk_size; offset < pixels; offset += chunk_size) {
					enqueueKernel2(queues[1 - index], kernels[1 - index], offset, std::min(chunk_size, pixels - offset)); // Queue other kernel execution
					enqueueRead2(queues[index], cl_buffer_data_out[index], chunk_size*lambda_pixel*data_out_per_point*sizeof(float), data_out[index]); // Wait for this kernel to finish execution
					// Extract output from this kernel, other kernel is already queued for execution so GPU is not waiting
					for(int i = 0; i < chunk_size*lambda_pixel; i++) {
						int output_index = (offset - chunk_size)*lambda_pixel + i;
						intensity[output_index] = data_out[index][i];
					}
					index = 1 - index;
				}
				enqueueRead2(queues[index], cl_buffer_data_out[index], (pixels%chunk_size)*lambda_pixel*data_out_per_point*sizeof(float), data_out[index]); // Wait for the last kernel to finish execution
				// Extract output from this kernel, other kernel is already queued for execution so GPU is not waiting
				for(int i = 0; i < (pixels%chunk_size)*lambda_pixel; i++) {
					int output_index = pixels/chunk_size*chunk_size*lambda_pixel + i;
						intensity[output_index] = data_out[index][i];
				}
				
				if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
				if (commrank==0) std::cout << "Time spent building frame and extracting output: " << std::chrono::duration<double>(time_now() - frame_start).count() << std::endl << std::flush;
				frame_time += std::chrono::duration<double>(time_now() - frame_start).count();
				
				#if (GPU_REGULAR_GRID_DEBUG == 1)
					enqueueRead2(queues[0], cl_buffer_debug, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(float), debug_buffer); // Wait for the last kernel to finish execution
					for(int i = 0; i < GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE; i++) {
						std::cout << i << "\t" << debug_buffer[i] << std::endl << std::flush;
					}
				#endif
				
				if (commrank==0) std::cout << "Constructing RenderCube ... " << std::flush;
				index = 0;
				for(int y = 0; y < y_pixel; y++) {
					for(int x = 0; x < x_pixel; x++) {
						for(int l = 0; l < lambda_pixel; l++) {
							newgrid.at(0).at(index) = xs[x];
							newgrid.at(1).at(index) = ys[y];
							newgrid.at(2).at(index) = lambdas[l];
							newdata.at(0).at(index) = 1e8*intensity[index];
							index++;
						}
					}
				}
				rendercube.setdata(newgrid, newdata);
				rendercube.setrendermethod("GPURegularGrid");
				rendercube.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
				if (lambda_pixel == 1) {
					rendercube.setobservationtype(FoMo::Imaging);
				} else {
					rendercube.setobservationtype(FoMo::Spectroscopic);
				}
				rendercube.setangles(l, b);
				if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
				
				}
				
				std::cout << "Frame time: " << frame_time/5 << std::endl << std::flush;
				
				if (commrank==0) std::cout << "Writing frame to file ... " << std::flush;
				std::stringstream ss;
				// if outfile is "", then this should not be executed.
				ss << outfile;
				ss << "l";
				ss << std::setfill('0') << std::setw(3) << std::round(l/pi*180.);
				ss << "b";
				ss << std::setfill('0') << std::setw(3) << std::round(b/pi*180.);
				ss << ".txt";
				rendercube.writegoftcube(ss.str());
				ss.str("");
				if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
				
			}
		}
				
		if (commrank==0) std::cout << "Freeing grid and data ... " << std::flush;
		newgrid.clear();
		newdata.clear();
		if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
		delete regular_grid;
		delete[] intensity;
		
		// Only returns last rendercube!
		return rendercube;
		
	}
}

