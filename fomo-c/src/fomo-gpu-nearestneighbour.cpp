#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <gsl/gsl_const_mksa.h>
#include <cmath>
#include <algorithm>
#include <chrono>

#include <CL/cl.hpp>

struct __attribute__ ((packed)) parameters_struct {
	cl_int ng;
    cl_float minx;
    cl_float dx;
    cl_float miny;
    cl_float dy;
    cl_float minz;
    cl_float dz;
    cl_int z_pixel;
    cl_int lambda_pixel;
    cl_float lambda0;
    cl_float lambda_width_in_A;
    cl_float speedoflight;
};

void quick_sort(int *temp_indices, const float *temp_coords, const int axis, int low, int high) {
	
	// Quick-sorts the section of temp_indices in [low, high[ based on the axis'th component of the corresponding coordinates in temp_coords
	
	if (high - low > 1) {
		
		// Pick random pivot in range
		float pivot = temp_coords[3*temp_indices[rand()%(high - low) + low] + axis];
		
		// Partition (Hoare)
		int i = low - 1;
		int j = high;
		int pivot_index;
		while (true) {
			do {
				i = i + 1;
			} while temp_coords[3*temp_indices[i] + axis] < pivot;
			do {
				j = j - 1;
			} while temp_coords[3*temp_indices[j] + axis] > pivot;
			if (i >= j) {
				pivot_index = j;
				break;
			}
			int temp = temp_indices[i];
			temp_indices[i] = temp_indices[j];
			temp_indices[j] = temp;
		}
		
		// Recurse
		quick_sort(temp_indices, temp_coords, axis, low, pivot_index);
		quick_sort(temp_indices, temp_coords, axis, pivot_index + 1, high);
		
	}
	
}

void construct_kd_tree(int *temp_indices, int *out_indices, const float *temp_coords, uint8_t *nodes, int low, int high, int current_node, float dx, float dy, float dz) {
	
	// Recursively constructs a KD-tree
	// temp_indices: array of size ng, contains current order of points
	// out_indices: array of size ng + 1, contains output order of points. Seemed to be too difficult to store in temp_indices, would require a lot of extra movement.
	// temp_coords: array of size 3*ng, contains coordinate values in order x0, y0, z0, x1, y1, ...
	// nodes: array of size ng + 1, stores the index of the axis of the split of this node, 3 for leaves
	// low, high: mark range of points to be considered in temp_indices (high is non-inclusive)
	// current_node: index of current node being constructed in out_indices and nodes
	// dx, dy, dz: range across three dimensions, can be used for splitting heuristic
	// The constructed tree is balanced, with any excess elements moved to the left and is stored using Eytzinger's method in out_indices and nodes.
	// Node indexing starts at 1 so that children of i are at 2*i and 2*i + 1. The output arrays out_indices and nodes have an unused element at the start.
	
	if (high - low > 1) {
		
		// Splitting heuristic: longest axis
		int axis;
		if (dx > dy && dx > dz) {
			axis = 0;
		} else if (dy > dz) {
			axis = 1;
		} else {
			axis = 2;
		}
		
		// Sort the current range and use median as node
		quick_sort(temp_indices, temp_coords, axis, low, high);
		int median_index = (high - low)/2 + low;
		
		// Store node
		out_indices[current_node] = temp_indices[median_index];
		nodes[current_node] = axis;
		
		// Recurse
		construct_kd_tree(temp_indices, out_indices, temp_coords, nodes, low, median_index, 2*current_node, );
		
		
	} else {
		// Base case
		out_indices[current_node] = temp_indices[low];
		nodes[current_node] = 3;
	}
	
}

FoMo::RenderCube gpunearestneighbourinterpolation(FoMo::GoftCube goftcube, const double l, const double b, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel,
	const double lambda_width) {

	// results is an array of at least dimension (x2-x1+1)*(y2-y1+1)*lambda_pixel and must be initialized to zero
	// determine contributions per pixel
	// Simple algorithm to start with: rotate data points on CPU, then let GPU iterate through all data points to find nearest for each sampling point.

	// Determine commrank (whether or not to print info)
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
#else
	commrank = 0;
#endif

	// Start timing
	auto start = std::chrono::high_resolution_clock::now();
	
	// Constants
	const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
	const int data_in_per_point = 3;
	const int data_out_per_point = 4;
	
	// Read data
	FoMo::tgrid grid = goftcube.readgrid();
	int ng = goftcube.readngrid();
	int dim = goftcube.readdim();
	FoMo::tphysvar peakvec=goftcube.readvar(0); // Peak intensity
	FoMo::tphysvar fwhmvec=goftcube.readvar(1); // line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);
	
	// Initialize OpenCL
	if (commrank==0) std::cout << "Initializing OpenCL ..." << std::flush;
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
	if(file.open("src/gpu-nearestneighbour.cl", std::ios_base::in | std::ios_base::binary) == NULL) {
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
	err = cl_program.build(cl_devices, "-cl-nv-verbose");
	if(err != CL_SUCCESS) {
		std::cerr << "Error: Could not compile OpenCL program!" << std::endl;
		std::cerr << "OpenCL build log:\n" << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl_devices[0]) << std::endl;
		exit(1);
	}
	// Allocate buffers
	unsigned int input_amount = ng;
	unsigned int output_amount = x_pixel*y_pixel*lambda_pixel;
	cl::Buffer cl_buffer_coords(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, 3*sizeof(float)*(input_amount + 1));
	cl::Buffer cl_buffer_nodes(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(uint8_t)*(input_amount + 1)); // tree metadata: index of axis of split for non-leaves, 3 for leaves
	cl::Buffer cl_buffer_data_in(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, data_in_per_point*sizeof(float)*input_amount); // peak, fwhm, losvel
	cl::Buffer cl_buffer_parameters(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(parameters_struct)); // various constant parameters
	cl::Buffer cl_buffer_data_out(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, data_out_per_point*sizeof(float)*output_amount);
	// Create queue
	cl::CommandQueue cl_queue(cl_context, cl_devices[0], 0, &err);
	assert(err == CL_SUCCESS);
	// Map buffers
	float *coords = (float*) cl_queue.enqueueMapBuffer(cl_buffer_coords, CL_TRUE, CL_MAP_WRITE, 0, 3*sizeof(float)*input_amount);
	float *nodes = (float*) cl_queue.enqueueMapBuffer(cl_buffer_nodes, CL_TRUE, CL_MAP_WRITE, 0, 3*sizeof(float)*input_amount);
	float *data_in = (float*) cl_queue.enqueueMapBuffer(cl_buffer_data_in, CL_TRUE, CL_MAP_WRITE, 0, data_in_per_point*sizeof(float)*input_amount);
	parameters_struct *parameters = (parameters_struct*) cl_queue.enqueueMapBuffer(cl_buffer_parameters, CL_TRUE, CL_MAP_WRITE, 0, sizeof(parameters_struct));
	float *data_out = (float*) cl_queue.enqueueMapBuffer(cl_buffer_data_out, CL_TRUE, CL_MAP_READ, 0, data_out_per_point*sizeof(float)*output_amount);
	// Create kernel
	cl::Kernel cl_kernel(cl_program, "calculate_ray", &err);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not create kernel: " << err << std::endl;
        exit(1);
    }
	cl_kernel.setArg(0, cl_buffer_coords);
	cl_kernel.setArg(1, cl_buffer_nodes);
	cl_kernel.setArg(2, cl_buffer_data_in);
	cl_kernel.setArg(3, cl_buffer_parameters);
	cl_kernel.setArg(4, cl_buffer_data_out);
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) std::cout << "Rotating coordinates to POS reference ..." << std::flush;
	double sinl = sin(l);
	double sinb = sin(b);
	double cosl = cos(l);
	double cosb = cos(b);
	double unit[] = {sinb*cosl, -sinb*sinl, cosb}; // Define the unit vector along the line-of-sight
	float temp_coords[3*input_amount];
	// Rotate every point ad store in temporary array, store extra data in OpenCL buffer
	for (int i = 0; i < ng; i++) {
		double x = grid[0][i];
		double y = grid[1][i];
		double z;
		if (dim == 2)
			z = 0;
		else
			z = grid[2][i];
		temp_coords[3*i] = x*cosb*cosl - y*cosb*sinl - z*sinb;
		temp_coords[3*i + 1] = x*sinl + y*cosl;
		temp_coords[3*i + 2] = x*sinb*cosl - y*sinb*sinl + z*cosb;
		data_in[data_in_per_point*i] = peakvec.at(i);
		data_in[data_in_per_point*i + 1] = fwhmvec.at(i);
		data_in[data_in_per_point*i + 2] = vx[i]*unit[0] + vy[i]*unit[1] + vz[i]*unit[2]; // velocity along line of sight for position [i]/[ng]
	}
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;
	
	// Construct tree
	int temp_indices[input_amount]; // Stores the order of coordinates without having to constantly move 3 floats per point
	int out_indices[input_amount]; // Stores the order of coordinates without having to constantly move 3 floats per point
	for(int i = 0; i < input_amount; i++) {
		temp_indices[i] = i;
	}
		
	// Compute the bounds of the input data points, so that we can equidistantly distribute the target pixels
	if (commrank==0) std::cout << "Calculating bounds of data points ..." << std::flush;
	float minx = coords[0];
	float maxx = coords[0];
	float miny = coords[1];
	float maxy = coords[1];
	float minz = coords[2];
	float maxz = coords[2];
	for (int i = 0; i < ng; i++) {
		minx = std::min(minx, coords[3*i]);
		maxx = std::max(maxx, coords[3*i]);
		miny = std::min(miny, coords[3*i + 1]);
		maxy = std::max(maxy, coords[3*i + 1]);
		minz = std::min(minz, coords[3*i + 2]);
		maxz = std::max(maxz, coords[3*i + 2]);
	}
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;
	
	// Storing constant parameters
	if (commrank==0) std::cout << "Storing constant parameters ..." << std::flush;
	parameters->ng = ng;
	parameters->minx = minx;
	parameters->dx = (maxx - minx)/(x_pixel - 1);
	parameters->miny = miny;
	parameters->dy = (maxy - miny)/(y_pixel - 1);
	parameters->minz = minz;
	parameters->dz = (maxz - minz);
	if (z_pixel != 1) parameters->dz /= (z_pixel - 1);
	parameters->z_pixel = z_pixel;
	parameters->lambda_pixel = lambda_pixel;
	parameters->lambda0 = goftcube.readlambda0(); // lambda0=AIA bandpass for AIA imaging
	parameters->lambda_width_in_A = lambda_width*goftcube.readlambda0()/speedoflight;
	parameters->speedoflight = speedoflight;
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;
	
	// Building frame
	
	if (commrank==0) std::cout << "Building frame ..." << std::flush;
	
	// Write input
	cl_queue.enqueueWriteBuffer(cl_buffer_coords, CL_FALSE, 0, 3*sizeof(float)*input_amount, coords);
	cl_queue.enqueueWriteBuffer(cl_buffer_nodes, CL_FALSE, 0, 3*sizeof(float)*input_amount, nodes);
	cl_queue.enqueueWriteBuffer(cl_buffer_data_in, CL_FALSE, 0, data_in_per_point*sizeof(float)*input_amount, data_in);
	cl_queue.enqueueWriteBuffer(cl_buffer_parameters, CL_FALSE, 0, sizeof(parameters_struct), parameters);
	
	// Run kernel
	err = cl_queue.enqueueNDRangeKernel(cl_kernel, cl::NullRange, cl::NDRange(x_pixel, y_pixel), cl::NullRange);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue the kernel for execution: " << err << std::endl;
        exit(1);
    }
	
	// Read output
	err = cl_queue.enqueueReadBuffer(cl_buffer_data_out, CL_TRUE, 0, data_out_per_point*sizeof(float)*output_amount, data_out);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue buffer read: " << err << std::endl;
        exit(1);
    }
	
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;
	
	// Extracting output
	if (commrank==0) std::cout << "Extracting output ..." << std::flush;
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel), yvec(x_pixel*y_pixel*lambda_pixel);
	FoMo::tvars newdata;
	FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel, 0);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	// Output buffer seems to be completely zero
	for(unsigned int i = 0; i < 20; i++) {
		// Check if buffer is non-empty
		//if (i < 3) std::cout << "Data in at index " << (3*ng - 3 + i) << " is: " << data_in[(3*ng - 3 + i)] << std::endl << std::flush;
		std::cout << "Data out at index " << i << " is: " << data_out[i] << std::endl << std::flush;
	}
	if (lambda_pixel > 1) {
		FoMo::tcoord lambdavec(x_pixel*y_pixel*lambda_pixel);
		newgrid.push_back(lambdavec);
	}
	for(unsigned int i = 0; i < output_amount; i++) {
		for(int j = 0; j < data_out_per_point - 1; j++) {
			newgrid.at(j).at(i) = data_out[data_out_per_point*i + j];
		}
		intens.at(i) = data_out[data_out_per_point*i + (data_out_per_point - 1)];
	}
	double pathlength;
	// this does not work if only one z_pixel is given (e.g. for a 2D simulation), or the maxz and minz are equal (face-on on 2D simulation)
	// assume that the thickness of the slab is 1Mm. 
	if ((maxz==minz) || (z_pixel==1)) {
		pathlength = 1.;
		std::cout << "Assuming that this is a 2D simulation: setting thickness of simulation to " << pathlength << "Mm." << std::endl << std::flush;
	} else {
		pathlength = (maxz-minz)/(z_pixel-1);
	}
	intens=FoMo::operator*(pathlength*1e8,intens); // assume that the coordinates are given in Mm, and convert to cm
	newdata.push_back(intens);
	FoMo::RenderCube rendercube(goftcube);
	rendercube.setdata(newgrid, newdata);
	rendercube.setrendermethod("GPUNearestNeighbour");
	rendercube.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
	if (lambda_pixel == 1) {
		rendercube.setobservationtype(FoMo::Imaging);
	} else {
		rendercube.setobservationtype(FoMo::Spectroscopic);
	}
	if (commrank==0) std::cout << " Done!" << std::endl << std::flush;
	
	// End timing
	auto end = std::chrono::high_resolution_clock::now();
	if (commrank==0) std::cout << "Time spent in nearestneighbourinterpolation (milliseconds): "
	<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl << std::flush;
	
	return rendercube;
}

namespace FoMo
{
	FoMo::RenderCube RenderWithGPUNearestNeighbour(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile) {
		
		// Optional pre-processing (building tree, ...)
		FoMo::RenderCube rendercube(goftcube);
		
		// Iterate through viewing angles
		const double pi=M_PI; // pi
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit)
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit) {
				rendercube = gpunearestneighbourinterpolation(goftcube,*lit,*bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
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
		
		return rendercube;
	}
}

