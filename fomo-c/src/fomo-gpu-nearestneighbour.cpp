#include <cstdlib>
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <cstring>
#include <string>
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

#define LEAF 3
#define WENT_DOWN_LEFT 0
#define WENT_DOWN_RIGHT 1
#define WENT_DOWN_BOTH 2

struct __attribute__ ((packed)) parameters_struct {
	cl_int offset;
	cl_int ng;
    cl_float minx;
    cl_float dx;
    cl_int x_pixel;
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

struct comparator {
	const float *temp_coords;
	const int axis;
	bool operator()(int a, int b) {
		return temp_coords[3*a + axis] < temp_coords[3*b + axis];
	}
};

void construct_kd_tree_rec(int *temp_indices, int *out_indices, const float *temp_coords, uint8_t *nodes, int low, int high, int current_node, int depth,
	float minx, float maxx, float miny, float maxy, float minz, float maxz, int ng, int max_depth) {
	
	// Recursively constructs a KD-tree
	// temp_indices: array of size ng, contains current order of points
	// out_indices: array of size ng + 1, contains output order of points. Seemed to be too difficult to store in temp_indices, would require a lot of extra movement.
	// low, high: mark range of points to be considered in temp_indices (high is non-inclusive)
	// current_node: index of current node being constructed in out_indices and nodes
	// depth: current depth, root is at 0
	// max_depth: depth of deepest leaf (pre-calculated, these leaves don't have to be in the tree yet)
	
	//std::cout << current_node << std::endl;
	
	int amount = high - low;
	if (amount > 1) {
		
		// Splitting heuristic: longest axis
		int axis; 
		float dx = maxx - minx;
		float dy = maxy - miny;
		float dz = maxz - minz;
		if (dx > dy && dx > dz) {
			axis = 0;
		} else if (dy > dz) {
			axis = 1;
		} else {
			axis = 2;
		}
		
		// Sort the current range and split so that all nodes end up in a contiguous address space
		//std::cout << "Test3" << std::endl;
		//quick_sort(temp_indices, temp_coords, axis, low, high);
		struct comparator comparator_obj = {temp_coords, axis};
		std::sort(&(temp_indices[low]), &(temp_indices[high]), comparator_obj);
		//std::cout << "Test4" << std::endl;
		int depth_dist = max_depth - depth - 1;
		int pivot_index = low + (1 << depth_dist) - 1 + std::max(0, std::min(1 << depth_dist, ng + 1 - (current_node << (depth_dist + 1))));
		
		// Store node
		out_indices[current_node] = temp_indices[pivot_index];
		nodes[current_node] = axis;
		
		float split = temp_coords[3*temp_indices[pivot_index] + axis];
		
		// Recurse
		if (axis == 0) {
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, low, pivot_index, 2*current_node, depth + 1, minx, split, miny, maxy, minz, maxz, ng, max_depth);
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, pivot_index + 1, high, 2*current_node + 1, depth + 1, split, maxx, miny, maxy, minz, maxz, ng, max_depth);
		} else if (axis == 1) {
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, low, pivot_index, 2*current_node, depth + 1, minx, maxx, miny, split, minz, maxz, ng, max_depth);
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, pivot_index + 1, high, 2*current_node + 1, depth + 1, minx, maxx, split, maxy, minz, maxz, ng, max_depth);
		} else {
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, low, pivot_index, 2*current_node, depth + 1, minx, maxx, miny, maxy, minz, split, ng, max_depth);
			construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, pivot_index + 1, high, 2*current_node + 1, depth + 1, minx, maxx, miny, maxy, split, maxz, ng, max_depth);
		}
		
	} else if (amount == 1) {
		// Base case, leaf
		out_indices[current_node] = temp_indices[low];
		nodes[current_node] = LEAF;
	} else {
		// Base case, no node
	}
	
}

void construct_kd_tree(int *out_indices, const float *temp_coords, uint8_t *nodes, int ng, float minx, float maxx, float miny, float maxy, float minz, float maxz) {
	
	// Recursively constructs a KD-tree
	// out_indices: array of size ng + 1, contains output order of points
	// temp_coords: array of size 3*ng, contains coordinate values in order x0, y0, z0, x1, y1, ...
	// nodes: array of size ng + 1, stores the index of the axis of the split of this node, 3 for leaves
	// ng: number of data points
	// minx, maxx, miny, maxy, minz, maxz: range across three dimensions, can be used for splitting heuristic
	// The constructed tree is balanced, with any excess elements moved to the left and is stored using Eytzinger's method in out_indices and nodes.
	// Node indexing starts at 1 so that children of i are at 2*i and 2*i + 1. The output arrays out_indices and nodes have an unused element at the start.
	
	
	int *temp_indices = new int[ng]; // Stores the order of coordinates without having to constantly move 3 floats per point
	for(int i = 0; i < ng; i++) {
		temp_indices[i] = i;
	}
	//std::cout << "Test1" << std::endl;
	construct_kd_tree_rec(temp_indices, out_indices, temp_coords, nodes, 0, ng, 1, 0, minx, maxx, miny, maxy, minz, maxz, ng, int(log2(ng)));
	//std::cout << "Test2" << std::endl;
	delete[] temp_indices;
	
}

// Test functions

float rand_float_in_range(float min, float max) {
	return (static_cast <float> (rand()) / static_cast <float> (RAND_MAX))*(max - min) + min;
}

inline float dist2(const float *coords, const float *target) {
	float sum = 0;
	for(int i = 0; i < 3; i++) {
		float diff = coords[i] - target[i];
		sum += diff*diff;
	}
	return sum;
}

int nearest_neighbour_brute_force(const float *coords, int ng, float x, float y, float z) {
	float target[] = {x, y, z};
	int min_index = -1;
	float min_dist2 = 0;
	for(int i = 1; i < ng + 1; i++) {
		float temp = dist2(&(coords[3*i]), target);
		if (temp < min_dist2 || i == 1) {
			min_index = i;
			min_dist2 = temp;
		}
	}
	return min_index;
}

int nearest_neighbour_kd(const float *coords, const uint8_t* nodes, int ng, float x, float y, float z) {
	
	// Initialization
	float target[] = {x, y, z};
	int max_depth = int(floor(log2(ng)));
	uint8_t stack[max_depth]; // Keeps track of path
	float diffs[max_depth]; // Stores non-absolute distances to splitting planes
	int stack_pointer = 0;
	int current_node = 1;
	int single_child_node = 0;
	if (ng%2 == 0) single_child_node = ng/2;
	int best_index = -1;
	float best_dist2 = 0;
	
	// Continue going down and up the tree until we go up the root, indicating that the search has finished
	while (current_node > 0) {
		
		// Descend until we reach a leaf
		while (true) {
			
			// Stop at leaf
			int axis = nodes[current_node];
			if (axis == LEAF)
				break;
				
			// Continue descending and keep track of path on stack
			diffs[stack_pointer] = target[axis] - coords[3*current_node + axis];
			if (diffs[stack_pointer] < 0) {
				// Left child
				current_node *= 2;
				stack[stack_pointer++] = (current_node == single_child_node ? WENT_DOWN_BOTH : WENT_DOWN_LEFT);
			} else {
				// Right child
				stack[stack_pointer] = WENT_DOWN_RIGHT;
				if (current_node != single_child_node) {
					current_node *= 2;
					current_node++;
					stack_pointer++;
				} else {
					break;
				}
			}
			
		}
		
		// Ascend until we reach the root of the tree or go down a new branch
		while (current_node > 0) {
			
			// Check if current node is best candidate
			float temp_dist = dist2(&(coords[3*current_node]), target);
			if (temp_dist < best_dist2 || best_index == -1) {
				best_index = current_node;
				best_dist2 = temp_dist;
			}
			
			// Consider branching out
			int axis = nodes[current_node];
			if (axis != LEAF && stack[stack_pointer] != WENT_DOWN_BOTH && diffs[stack_pointer]*diffs[stack_pointer] < best_dist2) {
				// Node is not a leaf, has an unexplored branch and splitting plane is within best distance, so branch out
				if (stack[stack_pointer] == WENT_DOWN_RIGHT) {
					// Left child
					current_node *= 2;
				} else {
					// Right child
					current_node *= 2;
					current_node++;
				}
				stack[stack_pointer++] = WENT_DOWN_BOTH;
				break;
			} else {
				// Continue going up
				current_node /= 2;
				stack_pointer--;
			}
			
		}
		
	}
	
	return best_index;
	
}

void test_nearest_neighbour_kd(const float *coords, const uint8_t* nodes, int ng, float minx, float maxx, float miny, float maxy, float minz, float maxz, int amount) {
	
	// Tests nearest-neighbour search through the kd-tree versus brute-force
	
	bool success = true;
	srand(10);
	for(int i = 0; i < amount; i++) {
		float x = rand_float_in_range(minx, maxx);
		float y = rand_float_in_range(miny, maxy);
		float z = rand_float_in_range(minz, maxz);
		float target[] = {x, y, z};
		int res1 = nearest_neighbour_brute_force(coords, ng, x, y, z);
		int res2 = nearest_neighbour_kd(coords, nodes, ng, x, y, z);
		float dist21 = dist2(&(coords[3*res1]), target);
		float dist22 = dist2(&(coords[3*res2]), target);
		if (dist21 != dist22) {
			success = false;
			std::cout << "Test " << i << " failed!" << std::endl;
			std::cout << "Input: " << x << " " << y << " " << z << std::endl;
			std::cout << "Res1: " << res1 << ", coordinates " << coords[3*res1] << " " << coords[3*res1 + 1] << " " << coords[3*res1 + 2] << std::endl;
			std::cout << "Distance squared: " << dist21 << std::endl;
			std::cout << "Res2: " << res2 << ", coordinates " << coords[3*res2] << " " << coords[3*res2 + 1] << " " << coords[3*res2 + 2] << std::endl;
			std::cout << "Distance squared: " << dist22 << std::endl;
			break;
		} else if (i%1000 == 0) {
			std::cout << "Finished " << i << " tests succesfully!" << std::endl << std::flush;
		}
	}
	if (success) {
		std::cout << "All " << amount << " tests succeeded!" << std::endl;
	} else {
		std::cout << "Not all " << amount << " tests succeeded!" << std::endl;
	}
	
}

// Rendering functions

static std::chrono::time_point<std::chrono::high_resolution_clock> time_now() {
	return std::chrono::high_resolution_clock::now();
}

void enqueueKernel(cl::CommandQueue queue, cl::Kernel kernel, int offset, int size) {
	int err = queue.enqueueNDRangeKernel(kernel, cl::NDRange(offset), cl::NDRange(size), cl::NullRange);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue the kernel for execution: " << err << std::endl;
        exit(1);
    }
}

void enqueueRead(cl::CommandQueue queue, cl::Buffer cl_buffer_data_out, int size, float *data_out) {
	int err = queue.enqueueReadBuffer(cl_buffer_data_out, CL_TRUE, 0, size, data_out);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue buffer read: " << err << std::endl;
        exit(1);
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
	std::chrono::time_point<std::chrono::high_resolution_clock> start = time_now();
	
	// Constants
	const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
	const int data_in_per_point = 3;
	const int data_out_per_point = 4;
	const int chunk_size = 1024; // Amount of jobs submitted to the GPU simultaneously
	
	// Read data
	FoMo::tgrid grid = goftcube.readgrid();
	int ng = goftcube.readngrid();
	//ng = 392751; // Setting ng to something small seems to cause memory access bugs in the kernel? Investigate?
	std::cout << "ng = " << ng << std::endl;
	int dim = goftcube.readdim();
	FoMo::tphysvar peakvec=goftcube.readvar(0); // Peak intensity
	FoMo::tphysvar fwhmvec=goftcube.readvar(1); // line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);
	
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
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) std::cout << "Rotating coordinates to POS reference ... " << std::flush;
	double sinl = sin(l);
	double sinb = sin(b);
	double cosl = cos(l);
	double cosb = cos(b);
	float *temp_coords = new float[3*input_amount];
	// Rotate every point and store in temporary array
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
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
		
	// Compute the bounds of the input data points, so that we can equidistantly distribute the target pixels
	if (commrank==0) std::cout << "Calculating bounds of data points ... " << std::flush;
	float minx = temp_coords[0];
	float maxx = temp_coords[0];
	float miny = temp_coords[1];
	float maxy = temp_coords[1];
	float minz = temp_coords[2];
	float maxz = temp_coords[2];
	for (int i = 0; i < ng; i++) {
		minx = std::min(minx, temp_coords[3*i]);
		maxx = std::max(maxx, temp_coords[3*i]);
		miny = std::min(miny, temp_coords[3*i + 1]);
		maxy = std::max(maxy, temp_coords[3*i + 1]);
		minz = std::min(minz, temp_coords[3*i + 2]);
		maxz = std::max(maxz, temp_coords[3*i + 2]);
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// We need some OpenCL buffers for the next step, wait for allocation to finish first
	if (commrank==0) std::cout << "Waiting for OpenCL to finish buffer allocation ... " << std::flush;
	queues[0].finish();
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// Construct tree
	if (commrank==0) std::cout << "Constructing KD-tree ... " << std::flush;
	int *out_indices = new int[ng + 1]; // Stores the order of coordinates without having to constantly move 3 floats per point
	construct_kd_tree(out_indices, temp_coords, nodes, ng, minx, maxx, miny, maxy, minz, maxz);
	if (commrank==0) std::cout << "Done! Time spent since start (milliseconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// Store coordinates and extra data
	if (commrank==0) std::cout << "Storing coordinates and extra data ... " << std::flush;
	double unit[] = {sinb*cosl, -sinb*sinl, cosb}; // Define the unit vector along the line-of-sight
	for(int i = 1; i < ng + 1; i++) {
		int original_index = out_indices[i];
		for(int j = 0; j < 3; j++) coords[3*i + j] = temp_coords[3*original_index + j];
		data_in[data_in_per_point*i] = peakvec.at(original_index);
		data_in[data_in_per_point*i + 1] = fwhmvec.at(original_index);
		data_in[data_in_per_point*i + 2] = vx[original_index]*unit[0] + vy[original_index]*unit[1] + vz[original_index]*unit[2]; // velocity along line of sight for position [i]/[ng]
	}
	delete[] temp_coords;
	delete[] out_indices;
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	//test_nearest_neighbour_kd(coords, nodes, ng, minx, maxx, miny, maxy, minz, maxz, 100000);
	//exit(0);
	
	// Storing constant parameters
	queues[1].finish();
	if (commrank==0) std::cout << "Storing constant parameters ... " << std::flush;
	for(int i = 0; i < 2; i++) {
		parameters[i]->ng = ng;
		parameters[i]->minx = minx;
		parameters[i]->dx = (maxx - minx)/(x_pixel - 1);
		parameters[i]->x_pixel = x_pixel;
		parameters[i]->miny = miny;
		parameters[i]->dy = (maxy - miny)/(y_pixel - 1);
		parameters[i]->minz = minz;
		parameters[i]->dz = (maxz - minz);
		if (z_pixel != 1) parameters[i]->dz /= (z_pixel - 1);
		parameters[i]->z_pixel = z_pixel;
		parameters[i]->lambda_pixel = lambda_pixel;
		parameters[i]->lambda0 = goftcube.readlambda0(); // lambda0=AIA bandpass for AIA imaging
		parameters[i]->lambda_width_in_A = lambda_width*goftcube.readlambda0()/speedoflight;
		parameters[i]->speedoflight = speedoflight;
	}
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// Building frame and extracting output
	
	if (commrank==0) std::cout << "Building frame and extracting output ... " << std::flush;
	
	// Write input
	queues[0].enqueueWriteBuffer(cl_buffer_coords, CL_FALSE, 0, 3*sizeof(float)*(input_amount + 1), coords);
	queues[0].enqueueWriteBuffer(cl_buffer_nodes, CL_FALSE, 0, sizeof(uint8_t)*(input_amount + 1), nodes);
	queues[0].enqueueWriteBuffer(cl_buffer_data_in, CL_FALSE, 0, data_in_per_point*sizeof(float)*(input_amount + 1), data_in);
	
	// Initialize output data extraction
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel), yvec(x_pixel*y_pixel*lambda_pixel);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	if (lambda_pixel > 1) {
		FoMo::tcoord lambdavec(x_pixel*y_pixel*lambda_pixel);
		newgrid.push_back(lambdavec);
	}
	FoMo::tvars newdata;
	FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel, 0);
	double pathlength;
	// this does not work if only one z_pixel is given (e.g. for a 2D simulation), or the maxz and minz are equal (face-on on 2D simulation)
	// assume that the thickness of the slab is 1Mm. 
	if ((maxz == minz) || (z_pixel == 1)) {
		pathlength = 1.;
		std::cout << "Assuming that this is a 2D simulation: setting thickness of simulation to " << pathlength << "Mm." << std::endl << std::flush;
	} else {
		pathlength = (maxz - minz)/(z_pixel - 1);
	}
	pathlength *= 1e8; // assume that the coordinates are given in Mm, and convert to cm
	
	// Ping-pong in between two kernels until work is finished
	queues[0].finish(); // The general queue needs to be finished because the second queue won't wait for it to copy the input
	int index = 0;
	parameters[index]->offset = 0;
	queues[index].enqueueWriteBuffer(cl_buffer_parameters[index], CL_FALSE, 0, sizeof(parameters_struct), parameters[index]);
	enqueueKernel(queues[index], kernels[index], 0, std::min(chunk_size, pixels));
	for(int offset = chunk_size; offset < pixels; offset += chunk_size) {
		parameters[1 - index]->offset = offset;
		queues[1 - index].enqueueWriteBuffer(cl_buffer_parameters[1 - index], CL_FALSE, 0, sizeof(parameters_struct), parameters[1 - index]);
		enqueueKernel(queues[1 - index], kernels[1 - index], offset, std::min(chunk_size, pixels - offset)); // Queue other kernel execution
		enqueueRead(queues[index], cl_buffer_data_out[index], chunk_size*lambda_pixel*data_out_per_point*sizeof(float), data_out[index]); // Wait for this kernel to finish execution
		// Extract output this kernel, other kernel is already queued for execution so GPU is not waiting
		for(int i = 0; i < chunk_size*lambda_pixel; i++) {
			int output_index = (offset - chunk_size)*lambda_pixel + i;
			for(int j = 0; j < data_out_per_point - 1; j++) {
				newgrid.at(j).at(output_index) = data_out[index][data_out_per_point*i + j];
			}
			intens.at(output_index) = data_out[index][data_out_per_point*i + (data_out_per_point - 1)]*pathlength;
		}
		index = 1 - index;
	}
	enqueueRead(queues[index], cl_buffer_data_out[index], (pixels%chunk_size)*lambda_pixel*data_out_per_point*sizeof(float), data_out[index]); // Wait for the last kernel to finish execution
	// Extract output this kernel, other kernel is already queued for execution so GPU is not waiting
	for(int i = 0; i < (pixels%chunk_size)*lambda_pixel; i++) {
		int output_index = pixels/chunk_size*chunk_size*lambda_pixel + i;
		for(int j = 0; j < data_out_per_point - 1; j++) {
			newgrid.at(j).at(output_index) = data_out[index][data_out_per_point*i + j];
		}
		intens.at(output_index) = data_out[index][data_out_per_point*i + (data_out_per_point - 1)]*pathlength;
	}
	
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// Constructing RenderCube
	if (commrank==0) std::cout << "Constructing RenderCube ... " << std::flush;
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
	if (commrank==0) std::cout << "Done! Time spent since start (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
	// End timing
	if (commrank==0) std::cout << "Time spent in gpunearestneighbourinterpolation (seconds): " << std::chrono::duration<double>(time_now() - start).count() << std::endl << std::flush;
	
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
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit) {
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit) {
				rendercube = gpunearestneighbourinterpolation(goftcube, *lit, *bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
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
		}
		
		return rendercube;
	}
}

