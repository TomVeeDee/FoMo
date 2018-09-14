#include <config.h>
#include "FoMo.h"
#include "FoMo-internal.h"

#include <string>

#ifdef HAVE_CL_CL_HPP

// This file simply stores the OpenCL kernel used by the RegularGridRenderer

std::string FoMo::RegularGridRenderer::readKernelSource() {
	return R"(
typedef struct __attribute__ ((packed)) Parameters {
	float rxx;
	float rxy;
	float rxz;
	float ryx;
	float ryy;
	float ryz;
	float rzx;
	float rzy;
	float rzz;
	float pixel_width;
	float pixel_height;
} Parameters;

// Macro variables that should be defined at compile-time:
// DEBUG: Specifies whether or not a debug buffer will be present
// ALL_INTENSITIES: Specifies whether or not the intensity for every wavelength should be calculated. If so, a float buffer should be used as output instead of a byte buffer.
// GAUSSIAN_PARAMETERS: Specifies whether or not only the integrated intensity needs to be calculated. If so, a float buffer (emissivity) should be given instead of a float8 buffer (points).
// MAX_EMISSIVITY: Used for GAUSSIAN_PARAMETERS display mode.
// MAX_DOPPLER_SHIFT: Used for GAUSSIAN_PARAMETERS display mode.
// MAX_SPECTRAL_WIDTH: Used for GAUSSIAN_PARAMETERS display mode.
// X_PIXEL: Amount of pixels along axes.
// LAMBDA_PIXEL, LAMBDA0
// MINX, MINY, MINZ, MAXX, MAXY, MAXZ: (non-integer) bounds of grid in global coordinates.
// GX, GY, GZ: size of a grid cell along each axis.
// GSX, GSY, GSZ: amount of grid cells along each axis.
// OX, OY: Represents the coordinates of the global origin in pixel coordinates (so usually X_PIXEL/2.0 and Y_PIXEL/2.0). Pixels are sampled at half-coordinates.
#define SPEEDOFLIGHT 299792458.0

inline int calculate_intersection(float4 pos, float* distance, float4 rz) {
	
	// Calculates the intersection in between the ray and box
	// If there is one, returns 1 and stores the distance along the ray in *distance, if not returns 0.
	// Standard algorithms don't work here because of the combination of macro variables and corner cases, ugly but compiles to fairly short code
	// Based on "An efficient and robust ray-box intersection algorithm" (https://dl.acm.org/citation.cfm?id=1198748)
	
	float tmax, tymin, tymax, tzmin, tzmax;
	if (rz.x == 0) {
		
		// Ray is parallel to X-planes
		if (rz.y == 0) {
			
			// Ray is parallel to X/Y-planes, can't be parallel to Z-planes too
			if ((pos.x > MINX) && (pos.x < MAXX) && (pos.y > MINY) && (pos.y < MAXY)) {
				// Ray intersects box
				if (rz.z > 0)
					*distance = MINZ - pos.z;
				else
					*distance = pos.z - MINZ;
				return 1;
			} else {
				return 0;
			}
			
		} else if (rz.z == 0) {
			
			// Ray is parallel to X/Z-planes but not parallel to Y-planes
			if (pos.x > MINX && pos.x < MAXX && pos.z > MINZ && pos.z < MAXZ) {
				// Ray intersects box
				if (rz.y > 0)
					*distance = MINY - pos.y;
				else
					*distance = pos.y - MINY;
				return 1;
			} else {
				return 0;
			}
			
		} else {
			
			// Ray is parallel to X-planes but not parallel to Y/Z-planes
			if (rz.y > 0) {
				*distance = (MINY - pos.y)/rz.y;
				tymax = (MAXY - pos.y)/rz.y;
			} else {
				*distance = (MAXY - pos.y)/rz.y;
				tymax = (MINY - pos.y)/rz.y;
			}
			if (rz.z > 0) {
				tzmin = (MINZ - pos.z)/rz.z;
				tzmax = (MAXZ - pos.z)/rz.z;
			} else {
				tzmin = (MAXZ - pos.z)/rz.z;
				tzmax = (MINZ - pos.z)/rz.z;
			}
			if ((*distance > tzmax) || (tzmin > tymax))
				return 0;
			if (tzmin < *distance)
				*distance = tzmin;
			return 1;
			
		}
		
	} else {
		
		// Ray is not parallel to X-planes
		if (rz.y == 0) {
			
			if (rz.z == 0) {
				
				// Ray is parallel to Y/Z-planes, but not parallel to X-planes
				if (pos.y > MINY && pos.y < MAXY && pos.z > MINZ && pos.z < MAXZ) {
					// Ray intersects box
					if (rz.x > 0)
						*distance = MINX - pos.x;
					else
						*distance = pos.x - MINX;
					return 1;
				} else {
					return 0;
				}
				
			} else {
				
				// Ray is parallel to Y-planes, but not parallel to X/Z-planes
				if (rz.x > 0) {
					*distance = (MINX - pos.x)/rz.x;
					tmax = (MAXX - pos.x)/rz.x;
				} else {
					*distance = (MAXX - pos.x)/rz.x;
					tmax = (MINX - pos.x)/rz.x;
				}
				if (rz.z > 0) {
					tzmin = (MINZ - pos.z)/rz.z;
					tzmax = (MAXZ - pos.z)/rz.z;
				} else {
					tzmin = (MAXZ - pos.z)/rz.z;
					tzmax = (MINZ - pos.z)/rz.z;
				}
				if ((*distance > tzmax) || (tzmin > tmax))
					return 0;
				if (tzmin < *distance)
					*distance = tzmin;
				return 1;
				
			}
			
		} else if (rz.z == 0) {
				
			// Ray is parallel to Z-planes, but not parallel to X/Y-planes
			if (rz.x > 0) {
				*distance = (MINX - pos.x)/rz.x;
				tmax = (MAXX - pos.x)/rz.x;
			} else {
				*distance = (MAXX - pos.x)/rz.x;
				tmax = (MINX - pos.x)/rz.x;
			}
			if (rz.y > 0) {
				tymin = (MINY - pos.y)/rz.y;
				tymax = (MAXY - pos.y)/rz.y;
			} else {
				tymin = (MAXY - pos.y)/rz.y;
				tymax = (MINY - pos.y)/rz.y;
			}
			if ((*distance > tymax) || (tymin > tmax))
				return 0;
			if (tymin < *distance)
				*distance = tymin;
			return 1;
			
		} else {
				
			// Ray is not parallel to X/Y/Z-planes
			if (rz.x > 0) {
				*distance = (MINX - pos.x)/rz.x;
				tmax = (MAXX - pos.x)/rz.x;
			} else {
				*distance = (MAXX - pos.x)/rz.x;
				tmax = (MINX - pos.x)/rz.x;
			}
			if (rz.y > 0) {
				tymin = (MINY - pos.y)/rz.y;
				tymax = (MAXY - pos.y)/rz.y;
			} else {
				tymin = (MAXY - pos.y)/rz.y;
				tymax = (MINY - pos.y)/rz.y;
			}
			if ((*distance > tymax) || (tymin > tmax))
				return 0;
			if (tymin > *distance)
				*distance = tymin;
			if (tymax < tmax)
				tmax = tymax;
			if (rz.z > 0) {
				tzmin = (MINZ - pos.z)/rz.z;
				tzmax = (MAXZ - pos.z)/rz.z;
			} else {
				tzmin = (MAXZ - pos.z)/rz.z;
				tzmax = (MINZ - pos.z)/rz.z;
			}
			if ((*distance > tzmax) || (tzmin > tmax))
				return 0;
			if (tzmin > *distance)
				*distance = tzmin;
			return 1;
			
		}
		
	}

}

kernel void calculate_ray(
#if (GAUSSIAN_PARAMETERS == 1)
	global const float8* gaussian_parameters,
#else
	global const float8* points,
#endif
#if (ALL_INTENSITIES == 1)
	constant const float* lambdaval,
#endif
	constant const Parameters* parameters,
#if (ALL_INTENSITIES == 1)
	global float* floats_out
#else
	global uchar* bytes_out
#endif
#if (DEBUG == 1)
	, global float* debug_buffer
#endif
) {
	
	// Calculate screen coordinates
	size_t id = get_global_id(0); // y*x_pixel + x
	int local_id = id - get_global_offset(0);
	float local_x = (convert_float(id%X_PIXEL) + (0.5 - OX))*parameters->pixel_width;
	float local_y = (convert_float(id/X_PIXEL) + (0.5 - OY))*parameters->pixel_height;
	
	// Initialize quantities
	#if (GAUSSIAN_PARAMETERS == 1)
		float emissivity = 0;
		float mu = 0; // We calculate the average relative to lambda0, so 0 means no Doppler shift
		float variance = 0;
	#else
		float intensity[LAMBDA_PIXEL];
		#pragma unroll
		for(int i = 0; i < LAMBDA_PIXEL; i++) {
			intensity[i] = 0;
		}
	#endif
	
	// Calculate first intersection with box
	float4 rx = (float4) (parameters->rxx, parameters->rxy, parameters->rxz, 0);
	float4 ry = (float4) (parameters->ryx, parameters->ryy, parameters->ryz, 0);
	float4 rz = (float4) (parameters->rzx, parameters->rzy, parameters->rzz, 0);
	float4 pos = rx*local_x + ry*local_y; // Position of ray at z_local = 0
	float t;
	#if (DEBUG == 1)
		int counter = 0;
	#endif
	if (calculate_intersection(pos, &t, rz) == 1) {
		
		// Ray intersects box at distance tmin
		
		// Initialization
		float4 intersection = pos + t*rz;
		int x = max(0, min(GSX - 1, convert_int((intersection.x - MINX)*(1.0/GX))));
		int y = max(0, min(GSY - 1, convert_int((intersection.y - MINY)*(1.0/GY))));
		int z = max(0, min(GSZ - 1, convert_int((intersection.z - MINZ)*(1.0/GZ))));
		float steps[3] = {fabs(GX/rz.x), fabs(GY/rz.y), fabs(GZ/rz.z)};
		float event_distance[3];
		if (rz.x > 0)
			event_distance[0] = (MINX + GX*(x + 1) - pos.x)/rz.x;
		else if (rz.x == 0)
			event_distance[0] = FLT_MAX;
		else
			event_distance[0] = (MINX + GX*x - pos.x)/rz.x;
		if (rz.y > 0)
			event_distance[1] = (MINY + GY*(y + 1) - pos.y)/rz.y;
		else if (rz.y == 0)
			event_distance[1] = FLT_MAX;
		else
			event_distance[1] = (MINY + GY*y - pos.y)/rz.y;
		if (rz.z > 0)
			event_distance[2] = (MINZ + GZ*(z + 1) - pos.z)/rz.z;
		else if (rz.z == 0)
			event_distance[2] = FLT_MAX;
		else
			event_distance[2] = (MINZ + GZ*z - pos.z)/rz.z;
		int in_bounds;
		int axis;
		do {
			
			// Find next event axis
			int index = y*(GSX*GSZ) + x*GSZ + z;
			if (event_distance[0] <= event_distance[1] && event_distance[0] <= event_distance[2]) {
				axis = 0;
				if (rz.x > 0) {
					x++;
					in_bounds = (x < GSX);
				} else {
					x--;
					in_bounds = (x >= 0);
				}
			} else if (event_distance[1] <= event_distance[2]) {
				axis = 1;
				if (rz.y > 0) {
					y++;
					in_bounds = (y < GSY);
				} else {
					y--;
					in_bounds = (y >= 0);
				}
			} else {
				axis = 2;
				if (rz.z > 0) {
					z++;
					in_bounds = (z < GSZ);
				} else {
					z--;
					in_bounds = (z >= 0);
				}
			}
			//if (in_bounds)
			//	prefetch(&(points[y*(GSX*GSZ) + x*GSZ + z]), 1);
			
			float distance = event_distance[axis] - t;
			
			#if (GAUSSIAN_PARAMETERS == 1)
				float8 point = gaussian_parameters[index];
				if (point.s0 != 0) {
					// Calculate Gaussian parameters using incremental formulas
					float input_emissivity = point.s0*distance;
					float emissivity2 = emissivity + input_emissivity;
					float input_mu = (rz.x*point.s2 + rz.y*point.s3 + rz.z*point.s4)*(-LAMBDA0/SPEEDOFLIGHT);
					float mu2 = (mu*emissivity + input_mu*input_emissivity)/emissivity2;
					float mu_diff1 = mu - mu2;
					float mu_diff2 = input_mu - mu2;
					variance = ((variance + mu_diff1*mu_diff1)*emissivity + (point.s1 + mu_diff2*mu_diff2)*input_emissivity)/emissivity2;
					emissivity = emissivity2;
					mu = mu2;
				}
				
			#else
			
				#if (LAMBDA_PIXEL > 1)
				
					// Pre-compute values for intensity calculations
					// Some more pre-computations could be done theoretically but the pre-emptive exponentiations sometimes hit the limit of single-precision
					float8 point = points[index];
					float fwhm = point.s1;
					float b = ((rz.x*point.s2 + rz.y*point.s3 + rz.z*point.s4)/fwhm)*(2*LAMBDA0/SPEEDOFLIGHT);
					float factor = point.s0*distance;
					float exponent0 = -b*b;
					float exponent1 = 4*b/fwhm;
					float exponent2 = -4/(fwhm*fwhm);
					
					#pragma unroll
					for(int i = 0; i < LAMBDA_PIXEL; i++) {
						float a = lambdaval[i];
						intensity[i] += factor*exp2(exponent0 + a*(exponent1 + a*exponent2));
					}
				
				#else
					intensity[0] += points[index].s0*distance;
				#endif
			#endif
			
			// Update traversal information
			t = event_distance[axis];
			event_distance[axis] += steps[axis];
			
			#if (DEBUG == 1)
				counter++;
			#endif
			
		} while (in_bounds);
		
	}
	
	#if (GAUSSIAN_PARAMETERS == 1)
	
		// Write out Gaussian parameters
		int index = 9*local_id;
		
		// Write out emissivity
		emissivity *= (3.0/MAX_EMISSIVITY);
		if (emissivity < 1) {
			bytes_out[index++] = max(0, min(255, convert_int(round(emissivity*255))));
			bytes_out[index++] = 0;
			bytes_out[index++] = 0;
		} else if (emissivity < 2) {
			bytes_out[index++] = 255;
			bytes_out[index++] = max(0, min(255, convert_int(round((emissivity - 1)*255))));
			bytes_out[index++] = 0;
		} else {
			bytes_out[index++] = 255;
			bytes_out[index++] = 255;
			bytes_out[index++] = max(0, min(255, convert_int(round((emissivity - 2)*255))));
		}
		
		// Write out Doppler shift
		mu *= (1.0/MAX_DOPPLER_SHIFT);
		if (mu < 0) {
			// Blue shift
			int val = max(0, min(255, convert_int(round((1 + mu)*255))));
			bytes_out[index++] = val;
			bytes_out[index++] = val;
			bytes_out[index++] = 255;
		} else {
			// Red shift
			int val = max(0, min(255, convert_int(round((1 - mu)*255))));
			bytes_out[index++] = 255;
			bytes_out[index++] = val;
			bytes_out[index++] = val;
		}
		
		// Write out spectral line width
		float spectral_width = sqrt(variance)*(3.0/MAX_SPECTRAL_WIDTH);
		if (emissivity == 0) {
			// Spectral line width is undefined, so set color to white
			bytes_out[index++] = 255;
			bytes_out[index++] = 255;
			bytes_out[index++] = 255;
		} else if (spectral_width < 1) {
			bytes_out[index++] = 0;
			bytes_out[index++] = 0;
			bytes_out[index++] = max(0, min(255, convert_int(round(spectral_width*255))));
		} else if (spectral_width < 2) {
			bytes_out[index++] = 0;
			bytes_out[index++] = max(0, min(255, convert_int(round((spectral_width - 1)*255))));
			bytes_out[index++] = 255;
		} else {
			bytes_out[index++] = max(0, min(255, convert_int(round((spectral_width - 2)*255))));
			bytes_out[index++] = 255;
			bytes_out[index++] = 255;
		}
		
	#else
		// Write out intensities
		int index = LAMBDA_PIXEL*local_id;
		#pragma unroll
		for(int i = 0; i < LAMBDA_PIXEL; i++) {
			floats_out[index++] = intensity[i];
		}
	#endif
	
}
)";
}

#endif
