#include <config.h>

#include <vector>
#include <string>

#ifdef HAVE_CL_CL_HPP
#include <CL/cl.hpp>
#include <chrono>
#include <iostream>
#include <boost/geometry/index/rtree.hpp>
#endif

namespace FoMo
{
	double readgoftfromchianti(const std::string chiantifile);
	DataCube readgoftfromchianti(const std::string chiantifile, std::string & ion, double & lambda0, double & atweight);
	GoftCube emissionfromdatacube(DataCube datacube, std::string chiantifile, std::string abundfile, const FoMoObservationType observationtype);
	
#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H	
	FoMo::RenderCube RenderWithCGAL(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithCGAL2D(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, const std::string outfile);
#endif
	
	FoMo::RenderCube RenderWithNearestNeighbour(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithProjection(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
#ifdef HAVE_CL_CL_HPP
	//FoMo::RenderCube RenderWithGPUNearestNeighbour(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	//std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithGPURegularGrid(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	typedef RegularGridRendererDisplayMode DisplayMode;
	
	class RegularGridRenderer {
	public:
		RegularGridRenderer(FoMo::GoftCube *goftcube);
		~RegularGridRenderer();
		void readBounds(float &minx, float &maxx, float &miny, float &maxy, float &minz, float &maxz);
		void constructRegularGrid(const int gridx, const int gridy, const int gridz, const float max_distance_x, const float max_distance_y,
			const float max_distance_z);
		void setRenderingSettings(const int x_pixel, const int y_pixel, const int lambda_pixel, const float lambda_width,  DisplayMode displayMode,
			const float max_emissivity = 1.0, const float max_doppler_shift = 1.0, const float max_spectral_width = 1.0);
		void renderToBuffer(const float l, const float b, const float view_width, const float view_height, unsigned char *data);
		void renderToCube(const float l, const float b, const float view_width, const float view_height, std::string fileName,
			FoMo::RenderCube *renderCubePointer = NULL);
	private:
		
		// Types
		typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
		typedef boost::geometry::model::box<point> box;
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
			cl_float pixel_width;
			cl_float pixel_height;
		} Parameters;
		
		// Constants
		const int chunk_size = 1024*2; // Amount of jobs submitted to the GPU simultaneously
		const int bytes_per_pixel = 9; // 3 images, 3 channels. The alpha channel is added during extraction
		
		// Temporary variables
		std::chrono::time_point<std::chrono::high_resolution_clock> start;
		
		// OpenCL variables
		cl_int err;
		cl::Context cl_context;
		std::vector<cl::Device> cl_devices;
		// Index: y*gridx*gridz + x*gridz + z
		// The peak and emissivity values stored here are already converted from per Mm to per cm (multiplied by 1e8)
		cl::Buffer cl_buffer_points;
		cl::Buffer cl_buffer_gaussian_parameters;
		cl::Buffer cl_buffer_lambdaval;
		cl::Buffer cl_buffer_parameters;
		cl::Buffer cl_buffer_bytes_out[2];
		cl::Buffer cl_buffer_floats_out[2];
		cl::Buffer cl_buffer_debug;
		cl::CommandQueue queues[2];
		cl_float *lambdaval;
		Parameters *parameters;
		cl_uchar *bytes_out[2];
		cl_float *floats_out[2];
		cl_float *debug_buffer;
		
		// Variables known from construction onward
		int commrank;
		FoMo::GoftCube *goftCube;
		float minx, maxx, miny, maxy, minz, maxz;
		boost::geometry::index::rtree<value, boost::geometry::index::quadratic<16>> rtree_var;
		
		// Variables known from regular grid construction onward
		bool hasRegularGrid = false;
		int gridx, gridy, gridz;
		float grid_size_x, grid_size_y, grid_size_z;
		float grid_mid_x, grid_mid_y, grid_mid_z;
		
		// Variables known from setting rendering settings onwards
		bool hasRenderingSettings = false;
		int x_pixel, y_pixel, lambda_pixel;
		float view_width, view_height, lambda_width;
		float max_emissivity; // Used for GaussianParameters display mode
		float max_doppler_shift; // Used for GaussianParameters display mode
		float max_spectral_width; // Used for GaussianParameters display mode
		float ox, oy;
		DisplayMode displayMode;
		cl::Kernel kernels[2];
		
		// Internal methods
		void setDisplayMode(DisplayMode displayMode);
		void render(const float l, const float b, const float view_width, const float view_height, unsigned char *bytes, float *floats);
		
		// Helper methods
		inline std::chrono::time_point<std::chrono::high_resolution_clock> time_now();
		inline void start_timing();
		inline void finish_timing(std::string message);
		std::string readKernelSource();
		inline void rotateAroundZ(float* in, float angle, float* out);
		inline void rotateAroundY(float* in, float angle, float* out);
		inline std::string float_to_string(float val);
		inline void enqueueKernel(int index, int offset, int size);
		inline void extractData(int index, int pixels_in_job, int offset, unsigned char *bytes, float *floats);
		inline void enqueueRead(cl::Buffer *cl_buffers, int index, int size, void *buffer);
	};
#endif
	
}
