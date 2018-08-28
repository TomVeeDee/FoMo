typedef struct __attribute__ ((packed)) Parameters {
	int offset;
} Parameters;

typedef struct __attribute__ ((packed)) RegularPoint {
	float peak;
	float fwhm;
	float vx;
	float vy;
	float vz;
} RegularPoint;

// Macro variables that should be defined at compile-time:
// X_PIXEL, Z_PIXEL, XZ_PIXEL, LAMBDA_PIXEL, LAMBDA0
// OX, OY: Represents the coordinates of the global origin in the local coordinate system.
// Rij for i, j = X, Y, Z: (RIX, RIY, RIZ) represents the I-vector from the local coordinate system in global coordinates.
// RX and RY do not have to be unit vectors in global coordinates (1*RX represents a distance of 1 pixel on the screen), RZ has to be a unit vector in global coordinates.
// SX, SY, SZ: signs of RZX, RZY, RZZ (-1/0/1).
// MINX, MINY, MINZ, MAXX, MAXY, MAXZ: (non-integer) bounds of grid in global coordinates.
// GX, GY, GZ: size of a grid cell along each axis.
// STEPX, STEPY, STEPZ: amount of distance that should be travelled along this ray in order to reach the next cell along this axis. STEPI = GI/RZI. Not necessary for SI == 0.
// GSX, GSY, GSZ: amount of grid cells along each axis.
#define SPEEDOFLIGHT 299792458.0

inline int calculate_intersection(float4 pos, float* distance) {
	
	// Calculates the intersection in between the ray and box
	// If there is one, returns 1 and stores the distance along the ray in *distance, if not returns 0.
	// Standard algorithms don't work here because of the combination of macro variables and corner cases, ugly but compiles to fairly short code
	// Based on "An efficient and robust ray-box intersection algorithm" (https://dl.acm.org/citation.cfm?id=1198748)
	
	#if (SX == 0)
		
		// Ray is parallel to X-planes
		#if (SY == 0)
			
			// Ray is parallel to X/Y-planes, can't be parallel to Z-planes too
			if (pos.x > MINX && pos.x < MAXX && pos.y > MINY && pos.y < MAXY) {
				// Ray intersects box
				#if (SZ == 1)
					*distance = MINZ - pos.z;
				#else
					*distance = pos.z - MINZ;
				#endif
				return 1;
			} else {
				return 0;
			}
			
		#elif (SZ == 0)
			
			// Ray is parallel to X/Z-planes but not parallel to Y-planes
			if (pos.x > MINX && pos.x < MAXX && pos.z > MINZ && pos.z < MAXZ) {
				// Ray intersects box
				#if (SY == 1)
					*distance = MINY - pos.y;
				#else
					*distance = pos.y - MINY;
				#endif
				return 1;
			} else {
				return 0;
			}
			
		#else
			
			// Ray is parallel to X-planes but not parallel to Y/Z-planes
			#if (SY == 1)
				*distance = (MINY - pos.y)/RZY;
				float tymax = (MAXY - pos.y)/RZY;
			#else
				*distance = (MAXY - pos.y)/RZY;
				float tymax = (MINY - pos.y)/RZY;
			#endif
			#if (SZ == 1)
				float tzmin = (MINZ - pos.z)/RZZ;
				float tzmax = (MAXZ - pos.z)/RZZ;
			#else
				float tzmin = (MAXZ - pos.z)/RZZ;
				float tzmax = (MINZ - pos.z)/RZZ;
			#endif
			if ((*distance > tzmax) || (tzmin > tymax))
				return 0;
			if (tzmin < *distance)
				*distance = tzmin;
			return 1;
			
		#endif
		
	#else
		
		// Ray is not parallel to X-planes
		#if (SY == 0)
			#if (SZ == 0)
				
				// Ray is parallel to Y/Z-planes, but not parallel to X-planes
				if (pos.y > MINY && pos.y < MAXY && pos.z > MINZ && pos.z < MAXZ) {
					// Ray intersects box
					#if (SX == 1)
						*distance = MINX - pos.x;
					#else
						*distance = pos.x - MINX;
					#endif
					return 1;
				} else {
					return 0;
				}
				
			#else
				
				// Ray is parallel to Y-planes, but not parallel to X/Z-planes
				#if (SX == 1)
					*distance = (MINX - pos.x)/RZX;
					float txmax = (MAXX - pos.x)/RZX;
				#else
					*distance = (MAXX - pos.x)/RZX;
					float txmax = (MINX - pos.x)/RZX;
				#endif
				#if (SZ == 1)
					float tzmin = (MINZ - pos.z)/RZZ;
					float tzmax = (MAXZ - pos.z)/RZZ;
				#else
					float tzmin = (MAXZ - pos.z)/RZZ;
					float tzmax = (MINZ - pos.z)/RZZ;
				#endif
				if ((*distance > tzmax) || (tzmin > txmax))
					return 0;
				if (tzmin < *distance)
					*distance = tzmin;
				return 1;
				
			#endif
		#elif (SZ == 0)
				
			// Ray is parallel to Z-planes, but not parallel to X/Y-planes
			#if (SX == 1)
				*distance = (MINX - pos.x)/RZX;
				float txmax = (MAXX - pos.x)/RZX;
			#else
				*distance = (MAXX - pos.x)/RZX;
				float txmax = (MINX - pos.x)/RZX;
			#endif
			#if (SY == 1)
				float tymin = (MINY - pos.y)/RZY;
				float tymax = (MAXY - pos.y)/RZY;
			#else
				float tymin = (MAXY - pos.y)/RZY;
				float tymax = (MINY - pos.y)/RZY;
			#endif
			if ((*distance > tymax) || (tymin > txmax))
				return 0;
			if (tymin < *distance)
				*distance = tymin;
			return 1;
			
		#else
				
			// Ray is not parallel to X/Y/Z-planes
			#if (SX == 1)
				*distance = (MINX - pos.x)/RZX;
				tmax = (MAXX - pos.x)/RZX;
			#else
				*distance = (MAXX - pos.x)/RZX;
				tmax = (MINX - pos.x)/RZX;
			#endif
			#if (SY == 1)
				tymin = (MINY - pos.y)/RZY;
				tymax = (MAXY - pos.y)/RZY;
			#else
				tymin = (MAXY - pos.y)/RZY;
				tymax = (MINY - pos.y)/RZY;
			#endif
			if ((*distance > tymax) || (tymin > tmax))
				return 0;
			if (tymin > *distance)
				*distance = tymin;
			if (tymax < tmax)
				tmax = tymax;
			#if (SZ == 1)
				tzmin = (MINZ - pos.z)/RZZ;
				tzmax = (MAXZ - pos.z)/RZZ;
			#else
				tzmin = (MAXZ - pos.z)/RZZ;
				tzmax = (MINZ - pos.z)/RZZ;
			#endif
			if ((*distance > tzmax) || (tzmin > tmax))
				return 0;
			if (tzmin > *distance)
				*distance = tzmin;
			return 1;
			
		#endif
		
	#endif

}

kernel void calculate_ray(global const RegularPoint* points, local const float* lambda_val, local const Parameters* parameters, global float* data_out) {
	
	// Calculate screen coordinates
	size_t id = get_global_id(0); // y*x_pixel + x
	int local_id = id - parameters->offset;
	float screen_x = convert_float(id%X_PIXEL) - OX;
	float screen_y = convert_float(id/X_PIXEL) - OY;
	
	// Initialize intensities
	float intensity[LAMBDA_PIXEL];
	//#pragma unroll
	for(int i = 0; i < LAMBDA_PIXEL; i++) {
		intensity[i] = 0;
	}
	
	// Calculate first intersection with box
	// Needed: cell[3], event_distance[3], t
	float4 rx = (float4) (RXX, RXY, RXZ, 0);
	float4 ry = (float4) (RYX, RYY, RYZ, 0);
	float4 ry = (float4) (RZX, RZY, RZZ, 0);
	float4 pos = rx*screen_x + ry*screen_y; // Position of ray at z_local = 0
	float t;
	if (calculate_intersection(pos, &t) == 1) {
		
		// Ray intersects box at distance tmin
		
		// Initialization
		float4 intersection = pos + t*rz;
		int x = convert_int((intersection.x - MINX)/GX);
		int y = convert_int((intersection.y - MINY)/GY);
		int z = convert_int((intersection.z - MINZ)/GZ)};
		float steps[3] = {STEPX, STEPY, STEPZ};
		float event_distance[3];
		#if (SX == 1)
			event_distance[0] = ((GX + 1)*cell[0] - pos.x)*STEPX;
		#elif (SX == 0)
			event_distance[0] = FLT_MAX;
		#else
			event_distance[0] = (pos.x - GX*cell[0])*STEPX;
		#endif
		#if (SY == 1)
			event_distance[1] = ((GY + 1)*cell[1] - pos.y)*STEPY;
		#elif (SY == 0)
			event_distance[1] = FLT_MAX;
		#else
			event_distance[1] = (pos.y - GY*cell[1])*STEPY;
		#endif
		#if (SZ == 1)
			event_distance[2] = ((GZ + 1)*cell[2] - pos.z)*STEPZ;
		#elif (SZ == 0)
			event_distance[2] = FLT_MAX;
		#else
			event_distance[2] = (pos.z - GZ*cell[2])*STEPZ;
		#endif
		uchar8 in_bounds = 1;
		int axis;
		do {
			
			// Find next event axis
			// Condition ordering can possibly be optimized during pre-processing by using S/STEP-macros
			int index = y*XZ_PIXEL + x*Z_PIXEL + z;
			if (event_distance[0] <= event_distance[1] && event_distance[0] <= event_distance[2]) {
				axis = 0;
				#if (SX == 1)
					x++;
					in_bounds = (x < GSX);
				#else
					x--;
					in_bounds = (x >= 0);
				#endif
			} else if (event_distance[1] <= event_distance[2]) {
				axis = 1;
				#if (SY == 1)
					y++;
					in_bounds = (y < GSY);
				#else
					y--;
					in_bounds = (y >= 0);
				#endif
			} else {
				axis = 2;
				#if (SZ == 1)
					z++;
					in_bounds = (z < GSZ);
				#else
					z--;
					in_bounds = (z >= 0);
				#endif
			}
			
			float distance = event_distance[axis] - t;
			t = event_distance[axis];
			event_distance[axis] += steps[axis];
			
			// Get values from cell
			float peak = points[index]->peak;
			float fwhm = points[index]->fwhm;
			float losvel = RZX*points[index]->vx + RZY*points[index]->vy + RZZ*points[index]->vz;
			
			// Pre-compute values for intensity calculations
			float b = -(losvel/fwhm)*2*LAMBDA0/SPEEDOFLIGHT;
			float factor = peak*distance;
			
			//#pragma unroll
			for(int i = 0; i < LAMBDA_PIXEL; i++) {
				float a = lambdaval[i]/fwhm;
				intensity[i] = 0;
			}
			
		} while (in_bounds);
		
	}
	
	// Write out intensities
	int index = DATA_OUT_PER_POINT*LAMBDA_PIXEL*local_id;
	//#pragma unroll
	for (int i = 0; i < LAMBDA_PIXEL; i++) {
		data_out[index++] = pos.x;
		data_out[index++] = pos.y;
		data_out[index++] = lambda_val[i] + LAMBDA0;
		data_out[index++] = intensity[i];
	}
	
}
	
/*	// Reading input
	size_t id = get_global_id(0); // y*x_pixel + x
	int local_id = id - params->offset;
	int x_pixel = params->x_pixel;
	int ng = params->ng;
	float pos[3];
	pos[0] = convert_float(id%x_pixel)*params->dx + params->minx;
	pos[1] = convert_float(id/x_pixel)*params->dy + params->miny;
	float minz = params->minz;
	float dz = params->dz;
	int z_pixel = params->z_pixel;
	int lambda_pixel = params->lambda_pixel;
	float lambda_width_in_A = params->lambda_width_in_A;
	float lambda0 = params->lambda0;
	float lambda1 = lambda_width_in_A/(lambda_pixel - 1);
	float lambda2 = lambda_width_in_A/2.;
	float lambda3 = lambda0/params->speedoflight;
	
	// Initialize intensity vector
	float intens[MAX_LAMBDA_PIXEL];
	float lambdaval[MAX_LAMBDA_PIXEL];
	for(int il = 0; il < lambda_pixel; il++) {
		intens[il] = 0;
		lambdaval[il] = convert_float(il)*lambda1 - lambda2;
	}
	
	// Iterate through samples
	uchar state[MAX_DEPTH]; // Keeps track of path
	float diffs[MAX_DEPTH]; // Stores non-absolute distances to splitting planes
	uchar axises[MAX_DEPTH]; // Stores non-absolute distances to splitting planes
	int single_child_node = 0;
	if (ng%2 == 0) single_child_node = ng/2;
	float best_dist2 = 0;
	int min_index = -1;
	for(int k = 0; k < z_pixel; k++) {
		
		pos[2] = convert_float(k)*dz + minz;
		
		// Find closest data point - KD-tree
		
		// Initialize KD-tree NN
		int stack_pointer = 0;
		int current_node = 1;
		if (k > 0) {
			// We're reusing the last nearest neighbour as best candidate so far
			best_dist2 = 0;
			for(int l = 0; l < 3; l++) {
				float temp = pos[l] - coords[3*min_index + l];
				best_dist2 += temp*temp;
			}
		}
		
		// Continue going down and up the tree until we go up the root, indicating that the search has finished
		while (current_node > 0) {
			
			// Descend until we reach a leaf
			while (1) {
				
				// Stop at leaf
				axises[stack_pointer] = nodes[current_node];
				int axis = axises[stack_pointer];
				if (axis == LEAF)
					break;
				
				// Continue descending and keep track of path on stack
				diffs[stack_pointer] = pos[axis] - coords[3*current_node + axis];
				if (diffs[stack_pointer] < 0) {
					// Left child
					state[stack_pointer++] = (current_node == single_child_node ? WENT_DOWN_BOTH : WENT_DOWN_LEFT);
					current_node = current_node << 1;
				} else {
					// Right child
					state[stack_pointer] = WENT_DOWN_RIGHT;
					if (current_node != single_child_node) {
						current_node = current_node << 1;
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
				// Excessive work if state[stack_pointer] == WENT_DOWN_BOTH, however keeping this in seems to be faster, possibly due to lower divergence?
				float temp_dist = 0;
				for(int l = 0; l < 3; l++) {
					float temp = pos[l] - coords[3*current_node + l];
					temp_dist += temp*temp;
				}
				if (temp_dist < best_dist2 || min_index == -1) {
					min_index = current_node;
					best_dist2 = temp_dist;
				}
				
				// Consider branching out
				if (axises[stack_pointer] != LEAF && state[stack_pointer] != WENT_DOWN_BOTH && diffs[stack_pointer]*diffs[stack_pointer] < best_dist2) {
					// Node is not a leaf, has an unexplored branch and splitting plane is within best distance, so branch out
					if (state[stack_pointer] == WENT_DOWN_RIGHT) {
						// Left child
						current_node = current_node << 1;
					} else {
						// Right child
						current_node = current_node << 1;
						current_node++;
					}
					state[stack_pointer++] = WENT_DOWN_BOTH;
					break;
				} else {
					// Continue going up
					current_node = current_node >> 1;
					stack_pointer--;
				}
				
			}
			
		}
		
		// End of KD-tree NN-search
		
		// Get data from closest point
		int index = DATA_IN_PER_POINT*min_index;
		float intpolpeak = data_in[index++];
		float intpolfwhm = data_in[index++];
		float intpollosvel = data_in[index];
		
		// Calculate intensities
		if (lambda_pixel > 1) {
			// Spectroscopic study
			float lambda4 = intpollosvel*lambda3;
			for (int il = 0; il < lambda_pixel; il++) {
				// if intpolpeak is not zero then the correct expression is used. otherwise, the intensity is just 0
				// it remains to be tested if this is faster than just the direct computation
				float temp = (lambdaval[il] - lambda4)/intpolfwhm;
				intens[il] += intpolpeak*exp(-temp*temp*2.772588722); // Last number is 4*ln(2)
			}
		} else {
			// Imaging study
			intens[0] += intpolpeak;
		}
	}
	
	// Write out intensity for every wavelength
	int index = DATA_OUT_PER_POINT*lambda_pixel*local_id;
	for (int il = 0; il < lambda_pixel; il++) {
		data_out[index++] = pos[0];
		data_out[index++] = pos[1];
		data_out[index++] = lambdaval[il] + lambda0;
		data_out[index++] = intens[il];
	}
	
}*/
