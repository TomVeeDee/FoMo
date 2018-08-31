typedef struct __attribute__ ((packed)) Parameters {
	int offset;
	float rxx;
	float rxy;
	float rxz;
	float ryx;
	float ryy;
	float ryz;
	float rzx;
	float rzy;
	float rzz;
	float x_offset;
	float y_offset;
} Parameters;

typedef struct __attribute__ ((packed)) RegularPoint {
	float peak;
	float fwhm;
	float vx;
	float vy;
	float vz;
} RegularPoint;

// Macro variables that should be defined at compile-time:
// DEBUG: Specifies whether or not a debug buffer will be present
// X_PIXEL: Amount of pixels along axes.
// PIXEL_WIDTH, PIXEL_HEIGHT: Dimensions of pixels in global coordinates.
// LAMBDA_PIXEL, LAMBDA0
// MINX, MINY, MINZ, MAXX, MAXY, MAXZ: (non-integer) bounds of grid in global coordinates.
// GX, GY, GZ: size of a grid cell along each axis.
// GSX, GSY, GSZ: amount of grid cells along each axis.
// X_OFFSET, Y_OFFSET: Represents the offset that should be added to the local coordinates to compensate for the shift that happened during the construction of the regular grid.
// OX, OY: Represents the coordinates of the global origin in pixel coordinates (so usually X_PIXEL/2.0 and Y_PIXEL/2.0). Pixels are sampled at half-coordinates.
// RIJ for I, J = X, Y, Z: (RIX, RIY, RIZ) represents the I-vector from the local coordinate system in global coordinates. All unit vectors.
// SX, SY, SZ: signs of RZX, RZY, RZZ (-1/0/1).
// STEPX, STEPY, STEPZ: amount of distance that should be travelled along this ray in order to reach the next cell along this axis. STEPI = abs(GI/RZI). Not necessary for SI == 0.
#define SPEEDOFLIGHT 299792458.0
#define DATA_OUT_PER_POINT 4

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
			if (pos.x > MINX && pos.x < MAXX && pos.y > MINY && pos.y < MAXY) {
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

#if (DEBUG == 1)
	kernel void calculate_ray(global const RegularPoint* points, constant const float* lambdaval, constant const Parameters* parameters, global float* data_out, global float* debug_buffer) {
#else
	kernel void calculate_ray(global const RegularPoint* points, constant const float* lambdaval, constant const Parameters* parameters, global float* data_out) {
#endif
	
	// Calculate screen coordinates
	size_t id = get_global_id(0); // y*x_pixel + x
	int local_id = id - parameters->offset;
	float local_x = (convert_float(id%X_PIXEL) + (0.5 - OX))*PIXEL_WIDTH;
	float local_y = (convert_float(id/X_PIXEL) + (0.5 - OY))*PIXEL_HEIGHT;
	
	// Initialize intensities
	float intensity[LAMBDA_PIXEL];
	//#pragma unroll
	for(int i = 0; i < LAMBDA_PIXEL; i++) {
		intensity[i] = 0;
	}
	
	// Calculate first intersection with box
	float4 rx = (float4) (parameters->rxx, parameters->rxy, parameters->rxz, 0);
	float4 ry = (float4) (parameters->ryx, parameters->ryy, parameters->ryz, 0);
	float4 rz = (float4) (parameters->rzx, parameters->rzy, parameters->rzz, 0);
	float4 pos = rx*local_x + ry*local_y; // Position of ray at z_local = 0
	float t;
	int bug_index = -1;
	int bug_x = -1;
	int bug_y = -1;
	int bug_z = -1;
	int bug_axis = -1;
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
		#if (DEBUG == 1)
			int counter = 0;
		#endif
		int in_bounds;
		int axis;
		do {
			
			// Find next event axis
			// Condition ordering can possibly be optimized during pre-processing by using S/STEP-macros
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
			
			#if (LAMBDA_PIXEL > 1)
			
				// Pre-compute values for intensity calculations
				// Some more pre-computations could be done theoretically but the pre-emptive exponentiations sometimes hit the limit of single-precision
				float fwhm = points[index].fwhm;
				float b = ((rz.x*points[index].vx + rz.y*points[index].vy + rz.z*points[index].vz)/fwhm)*(2*LAMBDA0/SPEEDOFLIGHT);
				float factor = points[index].peak*(event_distance[axis] - t);
				float exponent0 = -b*b;
				float exponent1 = 4*b/fwhm;
				float exponent2 = -4/(fwhm*fwhm);
				
				//#pragma unroll
				for(int i = 0; i < LAMBDA_PIXEL; i++) {
					float a = lambdaval[i];
					intensity[i] += factor*exp2(exponent0 + a*(exponent1 + a*exponent2));
					#if (DEBUG == 1)
						if (id == 7454 && i == 0) {
							debug_buffer[counter] = points[index].fwhm;
						}
					#endif
				}
			
			#else
				intensity[0] += points[index].peak;
			#endif
			
			// Update traversal information
			t = event_distance[axis];
			event_distance[axis] += steps[axis];
			
			#if (DEBUG == 1)
				counter++;
			#endif
			
		} while (in_bounds);
		
	}
	
	// Write out intensities
	int index = (DATA_OUT_PER_POINT*LAMBDA_PIXEL)*local_id;
	//#pragma unroll
	for(int i = 0; i < LAMBDA_PIXEL; i++) {
		data_out[index++] = local_x + parameters->x_offset;
		data_out[index++] = local_y + parameters->y_offset;
		data_out[index++] = lambdaval[i] + LAMBDA0;
		data_out[index++] = intensity[i];
	}
	
}
