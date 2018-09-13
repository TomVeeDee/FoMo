#include <config.h>
#include "FoMo.h"
#include "FoMo-internal.h"

#ifdef HAVE_CL_CL_HPP

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <set>
#include <chrono>
#include <limits>
#include <functional>

#include <CL/cl.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#define GPU_REGULAR_GRID_DEBUG 0
#define GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE 200

// Wrapper class

FoMo::RegularGridRendererWrapper::RegularGridRendererWrapper(FoMo::GoftCube *goftCube) {
	renderer = new RegularGridRenderer(goftCube);
}

FoMo::RegularGridRendererWrapper::~RegularGridRendererWrapper() {
	delete renderer;
}

void FoMo::RegularGridRendererWrapper::readBounds(float &minx, float &maxx, float &miny, float &maxy, float &minz, float &maxz) {
	renderer->readBounds(minx, maxx, miny, maxy, minz, maxz);
}

void FoMo::RegularGridRendererWrapper::constructRegularGrid(const int gridx, const int gridy, const int gridz, const float max_distance_x, const float max_distance_y,
	const float max_distance_z) {
	renderer->constructRegularGrid(gridx, gridy, gridz, max_distance_x, max_distance_y, max_distance_z);
}

void FoMo::RegularGridRendererWrapper::setRenderingSettings(const int x_pixel, const int y_pixel, const int lambda_pixel, const float lambda_width,
	const RegularGridRendererDisplayMode displayMode, const float max_intensity) {
	renderer->setRenderingSettings(x_pixel, y_pixel, lambda_pixel, lambda_width, displayMode, max_intensity);
}

void FoMo::RegularGridRendererWrapper::renderToBuffer(const float l, const float b, const float view_width, const float view_height, unsigned char *data) {
	renderer->renderToBuffer(l, b, view_width, view_height, data);
}

void FoMo::RegularGridRendererWrapper::renderToCube(const float l, const float b, const float view_width, const float view_height, std::string fileName, FoMo::RenderCube *renderCubePointer) {
	renderer->renderToCube(l, b, view_width, view_height, fileName, renderCubePointer);
}

// Public methods

/**
 * @brief Constructor for RegularGridRenderer.
 * 
 * Constructs a RegularGridRenderer with an internal R-tree for faster regular grid construction. The bounds on the input points are also already calculated.
 * 
 * @param goftCube The GoftCube pointer should remain valid throughout the lifetime of the renderer.
 */
FoMo::RegularGridRenderer::RegularGridRenderer(FoMo::GoftCube *goftCube) {
	
	// Constructs R-tree and calculates bounds of grid
	
	#ifdef HAVEMPI
		MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
	#else
		commrank = 0;
	#endif
	
	// Initialize OpenCL
	start_timing();
	// Find platforms
	std::vector<cl::Platform> cl_platforms;
	cl::Platform::get(&cl_platforms);
	if(cl_platforms.size() == 0) {
		std::cerr << "Error: No OpenCL platforms found!" << std::endl;
		exit(1);
	}
	// Make context
	cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(cl_platforms[0])(), 0};
	cl_context = cl::Context(CL_DEVICE_TYPE_ALL, cprops, NULL, NULL, &err);
	if(err != CL_SUCCESS) {
		std::cerr << "Error: Could not create OpenCL context!" << std::endl;
		exit(1);
	}
	// Find devices
	cl_devices = cl_context.getInfo<CL_CONTEXT_DEVICES>();
	if(cl_devices.size() == 0) {
		std::cerr << "Error: No OpenCL devices found!" << std::endl;
		exit(1);
	}
	// Allocate constant buffers
	cl_buffer_parameters = cl::Buffer(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(Parameters)); // Various constant parameters
	#if (GPU_REGULAR_GRID_DEBUG == 1)
		cl_buffer_debug = cl::Buffer(cl_context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(cl_float)); // Used to communicate data for kernel debugging
	#endif
	// Create queues
	queues[0] = cl::CommandQueue(cl_context, cl_devices[0], 0, &err); // Used for first kernel and general stuff
	queues[1] = cl::CommandQueue(cl_context, cl_devices[0], 0, &err); // Used for second kernel
	assert(err == CL_SUCCESS);
	// Map buffers
	parameters = (Parameters*) queues[0].enqueueMapBuffer(cl_buffer_parameters, CL_FALSE, CL_MAP_WRITE, 0, sizeof(Parameters));
	#if (GPU_REGULAR_GRID_DEBUG == 1)
		debug_buffer = (cl_float*) queues[0].enqueueMapBuffer(cl_buffer_debug, CL_FALSE, CL_MAP_READ, 0, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(cl_float));
	#endif
	finish_timing("Finished initializing OpenCL in ");
	
	// Initialization
	this->goftCube = goftCube;
	FoMo::tgrid grid = goftCube->readgrid();
	int ng = goftCube->readngrid();
	int dim = goftCube->readdim();
	
	// Preparing coordinates
	start_timing();
	std::vector<value> input_values(ng);
	// Build R-tree from gridpoints, this part is not parallel
	for (int i = 0; i < ng; i++) {
		//value boostpair = ;
		input_values[i] = std::make_pair(point(grid[0][i], grid[1][i], (dim == 2 ? 0 : grid[2][i])), i);
	}
	finish_timing("Finished preparing coordinates in ");
	
	// Build an rtree with the quadratic packing algorithm, it takes (slightly) more time to build, but queries are faster for large renderings
	rtree = boost::geometry::index::rtree<value, boost::geometry::index::quadratic<16>>(input_values.begin(), input_values.end());
	finish_timing("Finished building R-tree in ");
	
	// Compute bounds
	minx = *(min_element(grid[0].begin(), grid[0].end()));
	maxx = *(max_element(grid[0].begin(), grid[0].end()));
	miny = *(min_element(grid[1].begin(), grid[1].end()));
	maxy = *(max_element(grid[1].begin(), grid[1].end()));
	minz = *(min_element(grid[2].begin(), grid[2].end()));
	maxz = *(max_element(grid[2].begin(), grid[2].end()));
	finish_timing("Finished computing bounds in ");
	
}

/**
 * @brief Destructs this RegularGridRenderer and any necessary internal objects.
 */
FoMo::RegularGridRenderer::~RegularGridRenderer() {
	
}

/**
 * @brief Stores the bounds of the input data in the given parameters.
 */
void FoMo::RegularGridRenderer::readBounds(float &minx, float &maxx, float &miny, float &maxy, float &minz, float &maxz) {
	minx = this->minx;
	maxx = this->maxx;
	miny = this->miny;
	maxy = this->maxy;
	minz = this->minz;
	maxz = this->maxz;
}

/**
 * @brief Constructs the regular grid used for rendering based on the currently stored GoftCube.
 * 
 * The dimensions of the grid should be picked based on the resolution of the input data.
 * If the regular grid is too fine, some cells will not have an associated data point and will have zero emission.
 * 
 * @param gridx The amount of grid cells along the X-axis.
 * @param gridy The amount of grid cells along the Y-axis.
 * @param gridz The amount of grid cells along the Z-axis.
 */
void FoMo::RegularGridRenderer::constructRegularGrid(const int gridx, const int gridy, const int gridz, const float max_distance_x, const float max_distance_y, const float max_distance_z) {

	hasRegularGrid = true;
	hasRenderingSettings = false;
	start_timing();
	
	// Store grid parameters
	this->gridx = gridx;
	this->gridy = gridy;
	this->gridz = gridz;
	grid_mid_x = (minx + maxx)/2.0;
	grid_size_x = (maxx - minx)*float(gridx)/float(gridx - 1);
	grid_mid_y = (miny + maxy)/2.0;
	grid_size_y = (maxy - miny)*float(gridy)/float(gridy - 1);
	grid_mid_z = (minz + maxz)/2.0;
	grid_size_z = (goftCube->readdim() < 3 || gridz == 1 ? 1 : (maxz - minz)*float(gridz)/float(gridz - 1));
	if (goftCube->readdim() < 3 || gridz == 1) {
		std::cout << "Assuming that this is a 2D simulations: setting thickness of simulation to 1 Mm." << std::endl << std::flush;
	}
	
	// Allocate OpenCL buffers
	int input_size = gridx*gridy*gridz;
	cl_buffer_points = cl::Buffer(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float8)*input_size);
	cl_buffer_emissivity = cl::Buffer(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float)*input_size);
	cl_float8 *points = (cl_float8*) queues[0].enqueueMapBuffer(cl_buffer_points, CL_FALSE, CL_MAP_WRITE, 0, sizeof(cl_float8)*input_size);
	cl_float *emissivity = (cl_float*) queues[0].enqueueMapBuffer(cl_buffer_points, CL_FALSE, CL_MAP_WRITE, 0, sizeof(cl_float)*input_size);
	queues[0].finish();
	finish_timing("Finished allocating grid-dependent OpenCL buffers in ");
	
	// Iterate through cells
	double x, y, z;
	double deltaz = (maxz - minz);
	if (gridz != 1) deltaz /= (gridz - 1);
	point targetPoint;
	std::vector<value> returnedValues;
	box maxdistancebox;
	// Read the physical variables
	FoMo::tphysvar peakvec = goftCube->readvar(0); // Peak intensity 
	FoMo::tphysvar fwhmvec = goftCube->readvar(1); // Line width, = 1 for AIA imaging
	FoMo::tphysvar vx = goftCube->readvar(2);  
	FoMo::tphysvar vy = goftCube->readvar(3);
	FoMo::tphysvar vz = goftCube->readvar(4);
	int index;
	int counter = 0;
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic) collapse(2) shared (rtree, points, emissivity, counter) private (x, y, z, returnedValues, targetPoint, maxdistancebox, index)
	#endif
	for (int i = 0; i < gridy; i++) {
		for (int j = 0; j < gridx; j++) {
			
			y = double(i)/(gridy - 1)*(maxy - miny) + miny;
			x = double(j)/(gridx - 1)*(maxx - minx) + minx;
			
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k = 0; k < gridz; k++) {
				
				z = double(k)*deltaz + minz;
				
				// Search nearest point
				targetPoint = point(x, y, z);
				returnedValues.clear();
				// The second condition ensures the point is not further away than 
				// - half the x-resolution in the x-direction
				// - half the y-resolution in the y-direction
				// - the maximum of both the previous numbers in the z-direction (sort of improvising a convex hull approach)
				maxdistancebox = box(point(x - max_distance_x, y - max_distance_y, z - max_distance_z), point(x + max_distance_x, y + max_distance_y, z + max_distance_z));
				rtree.query(boost::geometry::index::nearest(targetPoint, 1) && boost::geometry::index::within(maxdistancebox), std::back_inserter(returnedValues));
				
				index = i*gridx*gridz + j*gridz + k;
				if (returnedValues.size() >= 1) {
					counter++;
					int nearestindex = returnedValues.at(0).second;
					points[index] = {float(1e8)*peakvec[nearestindex], fwhmvec[nearestindex], vx[nearestindex], vy[nearestindex], vz[nearestindex], 0, 0, 0};
					// Convert peak to emissivity, first constant is to convert from cm to Mm, last constant is sqrt(pi/(4*ln(2)))
					emissivity[index] = float(1e8)*peakvec[nearestindex]*fwhmvec[nearestindex]*1.064467019;
				} else {
					// All values must be initialized to deal with NaN and such
					points[index] = {0, 1, 0, 0, 0, 0, 0, 0};
					emissivity[index] = 0;
				}
				
			}
			
		}
	}
	
	std::cout << "Found matching data point within range for " << counter << " out of " << gridx*gridy*gridz << " grid points." << std::endl;
	
	finish_timing("Finished constructing regular grid in ");
	
	queues[0].enqueueWriteBuffer(cl_buffer_points, CL_FALSE, 0, sizeof(cl_float8)*input_size, points);
	queues[0].enqueueWriteBuffer(cl_buffer_emissivity, CL_FALSE, 0, sizeof(cl_float)*input_size, emissivity);
	finish_timing("Finished enqueuing grid-dependent OpenCL buffers write in ");
	
}

/**
 * @brief Set the rendering settings for this RegularGridRenderer.
 * 
 * Can only be called after constructRegularGrid has been called at least once.
 * @param x_pixel The amount of pixels along the X-axis.
 * @param y_pixel The amount of pixels along the Y-axis.
 * @param lambda_pixel The amount of wavelengths that should be sampled.
 * @param lambda_width The width of the wavelength window (in m/s).
 * @param max_intensity The maximal intensity that should be used for the IntegratedIntensity display mode.
 * @param displayMode The display mode that should be used for rendering.
 */
void FoMo::RegularGridRenderer::setRenderingSettings(const int x_pixel, const int y_pixel, const int lambda_pixel, const float lambda_width, DisplayMode displayMode,
	const float max_intensity) {
	
	// Stores and processes the given settings
	
	if (!hasRegularGrid) {
		std::cout << "Cannot call setRenderingSettings before calling constructRegularGrid! Returning." << std::endl << std::flush;
		return;
	}
	hasRenderingSettings = true;
	start_timing();
	
	// Store settings
	this->x_pixel = x_pixel;
	this->y_pixel = y_pixel;
	this->lambda_pixel = lambda_pixel;
	this->view_width = view_width;
	this->view_height = view_height;
	this->lambda_width = lambda_width;
	this->max_intensity = max_intensity;
	
	// Process settings
	// Allocate buffers
	ox = x_pixel/2.0;
	oy = y_pixel/2.0;
	int pixels = x_pixel*y_pixel;
	int output_amount = std::min(chunk_size, pixels)*lambda_pixel;
	cl_buffer_lambdaval = cl::Buffer(cl_context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float)*lambda_pixel);
	cl_buffer_bytes_out[0] = cl::Buffer(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, bytes_per_pixel*sizeof(cl_uchar)*output_amount);
	cl_buffer_bytes_out[1] = cl::Buffer(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, bytes_per_pixel*sizeof(cl_uchar)*output_amount);
	cl_buffer_floats_out[0] = cl::Buffer(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float)*output_amount);
	cl_buffer_floats_out[1] = cl::Buffer(cl_context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float)*output_amount);
	// Map buffers
	lambdaval = (cl_float*) queues[0].enqueueMapBuffer(cl_buffer_lambdaval, CL_FALSE, CL_MAP_WRITE, 0, sizeof(cl_float)*lambda_pixel);
	bytes_out[0] = (cl_uchar*) queues[0].enqueueMapBuffer(cl_buffer_bytes_out[0], CL_FALSE, CL_MAP_READ, 0, bytes_per_pixel*sizeof(cl_uchar)*output_amount);
	bytes_out[1] = (cl_uchar*) queues[1].enqueueMapBuffer(cl_buffer_bytes_out[1], CL_FALSE, CL_MAP_READ, 0, bytes_per_pixel*sizeof(cl_uchar)*output_amount);
	floats_out[0] = (cl_float*) queues[0].enqueueMapBuffer(cl_buffer_floats_out[0], CL_FALSE, CL_MAP_READ, 0, sizeof(cl_float)*output_amount);
	floats_out[1] = (cl_float*) queues[1].enqueueMapBuffer(cl_buffer_floats_out[1], CL_FALSE, CL_MAP_READ, 0, sizeof(cl_float)*output_amount);
	finish_timing("Finished setting and processing rendering settings in ");
	
	// Pre-compute and map lambda values
	float lambda_width_in_A = lambda_width*goftCube->readlambda0()/GSL_CONST_MKSA_SPEED_OF_LIGHT;
	for(int i = 0; i < lambda_pixel; i++)
		lambdaval[i] = float(i)/(lambda_pixel - 1)*lambda_width_in_A - lambda_width_in_A/2.0;
	queues[0].enqueueWriteBuffer(cl_buffer_lambdaval, CL_FALSE, 0, sizeof(float)*lambda_pixel, lambdaval);
	finish_timing("Finished pre-computing and mapping lambda values in ");
	
	// Set display mode
	setDisplayMode(displayMode);
	
}

/**
 * @brief Renders a frame from the given angle and stores it in the given buffer.
 * 
 * Can only be called after setRenderingSettings and when this RegularGridRenderer is in a screen display mode (so not AllIntensities).
 * @param l The l-angle.
 * @param b The b-angle.
 * @param view_width The width of the view in the space of the data points (in Mm).
 * @param view_height The height of the view in the space of the data points (in Mm).
 * @param data The output buffer. The caller must pass a buffer of appropriate size (right now, 1 byte is used per pixel). Indexing is done row by row, so index = y*x_pixel + x.
 */
void FoMo::RegularGridRenderer::renderToBuffer(const float l, const float b, const float view_width, const float view_height, unsigned char *data) {
	
	if (!hasRenderingSettings) {
		std::cout << "Cannot call render before calling setRenderingSettings! Returning." << std::endl << std::flush;
		return;
	}
	if (displayMode == DisplayMode::AllIntensities) {
		std::cout << "Cannot call render while in the AllIntensities display mode! Returning." << std::endl << std::flush;
		return;
	}
	render(l, b, view_width, view_height, data, NULL);
	
}

/**
 * @brief Renders a frame from the given angle and stores it in a RenderCube.
 * 
 * Output can be stored in a file, returned as a RenderCube or both.
 * Can only be called after setRenderingSettings.
 * Renders using the AllIntensities display mode, display mode gets set back to normal before the method returns.
 * @param l The l-angle.
 * @param b The b-angle.
 * @param view_width The width of the view in the space of the data points (in Mm).
 * @param view_height The height of the view in the space of the data points (in Mm).
 * @param fileName The name of the output file. Data will be stored in text format. If the fileName is empty, no file will be written (only useful if renderCubePointer is given).
 * @param renderCubePointer Optional parameter, if set to a non-null value the generated RenderCube will be stored in here as well.
 */
void FoMo::RegularGridRenderer::renderToCube(const float l, const float b, const float view_width, const float view_height, std::string fileName, FoMo::RenderCube *renderCubePointer) {
	
	// Check if output is necessary
	if (fileName.empty() && renderCubePointer == NULL)
		return;
	
	if (!hasRenderingSettings) {
		std::cout << "Cannot call renderToFile before calling setRenderingSettings! Returning." << std::endl << std::flush;
		return;
	}
	
	// Initialization
	DisplayMode tmpDisplayMode = displayMode;
	if (tmpDisplayMode != DisplayMode::AllIntensities)
		setDisplayMode(DisplayMode::AllIntensities);
	start_timing();
	float *data = new float[x_pixel*y_pixel*lambda_pixel];
	finish_timing("Finished allocating frame buffer in ");
	
	// Render
	render(l, b, view_width, view_height, NULL, data);
	finish_timing("Finished rendering frame in ");
	
	// Initialize RenderCube
	FoMo::RenderCube renderCube(*goftCube);
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
	finish_timing("Finished initializing RenderCube in ");
	
	// Extract data
	float global_offset_vector[] = {grid_mid_x, grid_mid_y, grid_mid_z};
	float temp[3];
	float local_offset_vector[3];
	rotateAroundZ(global_offset_vector, l, temp);
	rotateAroundY(temp, -b, local_offset_vector);
	float xs[x_pixel];
	for(int x = 0; x < x_pixel; x++)
		xs[x] = (x + 0.5 - ox)*(view_width/x_pixel) + local_offset_vector[0];
	float ys[y_pixel];
	for(int y = 0; y < y_pixel; y++)
		ys[y] = (y + 0.5 - oy)*(view_height/y_pixel) + local_offset_vector[1];
	float lambdas[lambda_pixel];
	float lambda0 = goftCube->readlambda0();
	for(int l = 0; l < lambda_pixel; l++)
		lambdas[l] = lambdaval[l] + lambda0;
	int index = 0;
	for(int y = 0; y < y_pixel; y++) {
		for(int x = 0; x < x_pixel; x++) {
			for(int l = 0; l < lambda_pixel; l++) {
				newgrid[0][index] = xs[x];
				newgrid[1][index] = ys[y];
				newgrid[2][index] = lambdas[l];
				newdata[0][index] = data[index];
				index++;
			}
		}
	}
	finish_timing("Finished extracting data in ");
	
	// Construct RenderCube
	renderCube.setdata(newgrid, newdata);
	renderCube.setrendermethod("GPURegularGrid");
	renderCube.setresolution(x_pixel, y_pixel, gridz, lambda_pixel, lambda_width);
	if (lambda_pixel == 1) {
		renderCube.setobservationtype(FoMo::Imaging);
	} else {
		renderCube.setobservationtype(FoMo::Spectroscopic);
	}
	renderCube.setangles(l, b);
	finish_timing("Finished constructing RenderCube in ");
	
	// Write to file if necessary
	if (!fileName.empty()) {
		renderCube.writegoftcube(fileName);
		finish_timing("Finished writing RenderCube to file in ");
	}
	
	// Copy RenderCube if necessary
	if (renderCubePointer != NULL) {
		*renderCubePointer = renderCube;
		finish_timing("Finished copying RenderCube in ");
	}
	
	// Free memory and switch back to old display mode if necessary
	delete[] data;
	if (tmpDisplayMode != DisplayMode::AllIntensities)
		setDisplayMode(tmpDisplayMode);
}

// Internal methods

void FoMo::RegularGridRenderer::setDisplayMode(DisplayMode displayMode) {
	
	// Stores and processes the given display mode
	
	start_timing();
	
	// Store display mode
	this->displayMode = displayMode;
	
	// Process display mode
	// Load program. It seems to be necessary to construct a new program object for every build to avoid error code CL_INVALID_OPERATION.
	cl::Program cl_program(cl_context, readKernelSource());
	// Compile program
	std::string build_options = "-cl-nv-verbose -D DEBUG=" + std::to_string(GPU_REGULAR_GRID_DEBUG)
		+ " -D ALL_INTENSITIES=" + std::to_string(displayMode == DisplayMode::AllIntensities)
		+ " -D INTEGRATED_INTENSITY=" + std::to_string(displayMode == DisplayMode::IntegratedIntensity)
		+ " -D MAX_INTENSITY=" + float_to_string(max_intensity)
		+ " -D X_PIXEL=" + std::to_string(x_pixel)
		+ " -D LAMBDA_PIXEL=" + std::to_string(lambda_pixel) + " -D LAMBDA0=" + float_to_string(float(goftCube->readlambda0()))
		+ " -D MINX=" + float_to_string(-grid_size_x/2) + " -D MAXX=" + float_to_string(grid_size_x/2)
		+ " -D MINY=" + float_to_string(-grid_size_y/2) + " -D MAXY=" + float_to_string(grid_size_y/2)
		+ " -D MINZ=" + float_to_string(-grid_size_z/2) + " -D MAXZ=" + float_to_string(grid_size_z/2)
		+ " -D GX=" + float_to_string(grid_size_x/gridx) + " -D GSX=" + std::to_string(gridx)
		+ " -D GY=" + float_to_string(grid_size_y/gridy) + " -D GSY=" + std::to_string(gridy)
		+ " -D GZ=" + float_to_string(grid_size_z/gridz) + " -D GSZ=" + std::to_string(gridz)
		+ " -D OX=" + float_to_string(ox) + " -D OY=" + float_to_string(oy);
	err = cl_program.build(cl_devices, build_options.c_str());
	if(err != CL_SUCCESS) {
		std::cerr << "Error: Could not compile OpenCL program! Error code: " << err << std::endl;
		std::cerr << "OpenCL build log:\n" << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl_devices[0]) << std::endl;
		exit(1);
	}
	finish_timing("Finished compiling OpenCL program in ");
	// Create kernels
	kernels[0] = cl::Kernel(cl_program, "calculate_ray", &err);
	if(err != CL_SUCCESS) {
		std::cerr << "Error: Could not create kernel: " << err << std::endl;
		exit(1);
	}
	kernels[1] = cl::Kernel(cl_program, "calculate_ray", &err);
	if(err != CL_SUCCESS) {
		std::cerr << "Error: Could not create kernel: " << err << std::endl;
		exit(1);
	}
	finish_timing("Finished creating kernels in ");
	
	// Set kernel arguments
	for(int i = 0; i < 2; i++) {
		if (displayMode == DisplayMode::IntegratedIntensity)
			kernels[i].setArg(0, cl_buffer_emissivity);
		else
			kernels[i].setArg(0, cl_buffer_points);
		kernels[i].setArg(1, cl_buffer_lambdaval);
		kernels[i].setArg(2, cl_buffer_parameters);
		if (displayMode == DisplayMode::AllIntensities)
			kernels[i].setArg(3, cl_buffer_floats_out[i]);
		else
			kernels[i].setArg(3, cl_buffer_bytes_out[i]);
		#if (GPU_REGULAR_GRID_DEBUG == 1)
			kernels[i].setArg(4, cl_buffer_debug);
		#endif
	}
	finish_timing("Finished settings kernel arguments in ");
	
}

void FoMo::RegularGridRenderer::render(const float l, const float b, const float view_width, const float view_height, unsigned char *bytes, float *floats) {
	
	// Used by both public render methods
	// If the display mode is AllIntensities, floats cannot be NULL, else bytes cannot be NULL
	// Indexing is done first on y, then on x, so index = y*x_pixel + x. Format of single pixel depends on display mode.
	// The returned intensities are already converted from per Mm to per cm
	
	// Calculate and write frame parameters
	float temp[3];
	float inx[] = {1, 0, 0};
	float iny[] = {0, 1, 0};
	float inz[] = {0, 0, 1};
	float rx[3];
	float ry[3];
	float rz[3];
	rotateAroundY(inx, b, temp);
	rotateAroundZ(temp, -l, rx);
	rotateAroundZ(iny, -l, ry);
	rotateAroundY(inz, b, temp);
	rotateAroundZ(temp, -l, rz);
	parameters->rxx = rx[0]; parameters->rxy = rx[1]; parameters->rxz = rx[2];
	parameters->ryx = ry[0]; parameters->ryy = ry[1]; parameters->ryz = ry[2];
	parameters->rzx = rz[0]; parameters->rzy = rz[1]; parameters->rzz = rz[2];
	parameters->pixel_width = view_width/x_pixel; parameters->pixel_height = view_height/y_pixel;
	queues[0].enqueueWriteBuffer(cl_buffer_parameters, CL_FALSE, 0, sizeof(Parameters), parameters);
	
	// Ping-pong in between two kernels until work is finished
	queues[0].finish();
	queues[1].finish();
	int index = 0;
	int pixels = x_pixel*y_pixel;
	enqueueKernel(index, 0, std::min(chunk_size, pixels));
	for(int offset = chunk_size; offset < pixels; offset += chunk_size) {
		enqueueKernel(1 - index, offset, std::min(chunk_size, pixels - offset)); // Queue other kernel execution
		// Wait for this kernel to finish execution and extract its output, other kernel is already queued for execution so GPU is not waiting
		extractData(index, chunk_size, offset - chunk_size, bytes, floats);
		index = 1 - index;
	}
	// Wait for the last kernel to finish execution and extract its output
	extractData(index, pixels%chunk_size, pixels/chunk_size*chunk_size, bytes, floats);
	
	#if (GPU_REGULAR_GRID_DEBUG == 1)
		err = queues[0].enqueueReadBuffer(cl_buffer_debug, CL_TRUE, 0, GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE*sizeof(cl_float), debug_buffer);
		if(err != CL_SUCCESS) {
			std::cerr << "Error: Could not enqueue debug buffer read: " << err << std::endl;
			exit(1);
		}
		for(int i = 0; i < GPU_REGULAR_GRID_DEBUG_BUFFER_SIZE; i++) {
			std::cout << i << "\t" << debug_buffer[i] << std::endl << std::flush;
		}
	#endif
	
}

// Helper methods

inline std::chrono::time_point<std::chrono::high_resolution_clock> FoMo::RegularGridRenderer::time_now() {
	return std::chrono::high_resolution_clock::now();
}

inline void FoMo::RegularGridRenderer::start_timing() {
	start = time_now();
}

inline void FoMo::RegularGridRenderer::finish_timing(std::string message) {
	// Also resets start
	if (commrank==0) std::cout << message << std::chrono::duration<double>(time_now() - start).count() << " seconds." << std::endl << std::flush;
	start_timing();
}

inline void FoMo::RegularGridRenderer::rotateAroundZ(float* in, float angle, float* out) {
	
	// in/out are 3-component vectors
	// angle is in radians
	
	float cosa = cos(angle);
	float sina = sin(angle);
	out[0] = cosa*in[0] - sina*in[1];
	out[1] = sina*in[0] + cosa*in[1];
	out[2] = in[2];
	
}

inline void FoMo::RegularGridRenderer::rotateAroundY(float* in, float angle, float* out) {
	
	// in/out are 3-component vectors
	// angle is in radians
	
	float cosa = cos(angle);
	float sina = sin(angle);
	out[0] = cosa*in[0] + sina*in[2];
	out[1] = in[1];
	out[2] = -sina*in[0] + cosa*in[2];
	
}

inline std::string FoMo::RegularGridRenderer::float_to_string(float val) {
	// Replaces commas with dots
	std::string s = std::to_string(val);
	std::replace(s.begin(), s.end(), ',', '.');
	return s;
}

inline void FoMo::RegularGridRenderer::enqueueKernel(int index, int offset, int size) {
	// Enqueues the kernel with the given index for execution
	err = queues[index].enqueueNDRangeKernel(kernels[index], cl::NDRange(offset), cl::NDRange(size), cl::NullRange);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue the kernel for execution: " << err << std::endl;
        exit(1);
    }
}

inline void FoMo::RegularGridRenderer::extractData(int index, int pixels_in_job, int offset, unsigned char *bytes, float *floats) {
	
	// Extracts data into the correct buffer depending on the display mode
	// pixels is the total amount of pixels being processed in the corresponding job
	// offset is the pixel index of the first pixel being processed in the corresponding job
    
	if (displayMode == DisplayMode::AllIntensities) {
		// Load floats
		enqueueRead(cl_buffer_floats_out, index, pixels_in_job*lambda_pixel*sizeof(cl_float), floats_out[index]);
		int local_offset = offset*lambda_pixel;
		for(int i = 0; i < pixels_in_job*lambda_pixel; i++)
			floats[local_offset + i] = floats_out[index][i];
	} else {
		// Load bytes
		enqueueRead(cl_buffer_bytes_out, index, pixels_in_job*bytes_per_pixel*sizeof(cl_uchar), bytes_out[index]);
		int local_offset = offset*bytes_per_pixel;
		for(int i = 0; i < pixels_in_job*bytes_per_pixel; i++)
			bytes[local_offset + i] = bytes_out[index][i];
	}
	
}

inline void FoMo::RegularGridRenderer::enqueueRead(cl::Buffer *cl_buffers, int index, int size, void *buffer) {
	// Blocks and reads output from the output buffer with the given index
	err = queues[index].enqueueReadBuffer(cl_buffers[index], CL_TRUE, 0, size, buffer);
	if(err != CL_SUCCESS) {
        std::cerr << "Error: Could not enqueue buffer read: " << err << std::endl;
        exit(1);
    }
}

namespace FoMo {
	
	FoMo::RenderCube RenderWithGPURegularGrid(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
		std::vector<double> lvec, std::vector<double> bvec, std::string outfile) {
		
		// Renders frames using the RegularGridRenderer API
		// Only returns the last generated RenderCube
		
		// Pre-processing
		FoMo::RenderCube renderCube(goftcube);
		RegularGridRenderer renderer = RegularGridRenderer(&goftcube);
		float minx, maxx, miny, maxy, minz, maxz;
		renderer.readBounds(minx, maxx, miny, maxy, minz, maxz);
		float maxdistance = 2*std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/.3; // Taken from NearestNeighbour
		renderer.constructRegularGrid(x_pixel, y_pixel, z_pixel, maxdistance, maxdistance, maxdistance);
		renderer.setRenderingSettings(x_pixel, y_pixel, lambda_pixel, lambda_width, DisplayMode::AllIntensities);
		
		for (std::vector<double>::iterator lit = lvec.begin(); lit != lvec.end(); ++lit) {
			for (std::vector<double>::iterator bit = bvec.begin(); bit != bvec.end(); ++bit) {
				
				// Render one frame if an output file is given or this 
				double l = *lit;
				double b = *bit;
				std::stringstream ss;
				ss << outfile;
				ss << "l";
				ss << std::setfill('0') << std::setw(3) << std::round(l/M_PI*180.);
				ss << "b";
				ss << std::setfill('0') << std::setw(3) << std::round(b/M_PI*180.);
				ss << ".txt";
				renderer.renderToCube(l, b, maxx - minx, maxy - miny, outfile.empty() ? "" : ss.str(), std::next(lit) == lvec.end() && std::next(bit) == bvec.end() ? &renderCube : NULL);
				
			}
		}
		
		return renderCube;
		
	}
	
}

#endif
