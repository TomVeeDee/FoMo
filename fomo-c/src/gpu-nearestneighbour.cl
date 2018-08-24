struct __attribute__ ((packed)) parameters_struct {
	int offset;
	int ng;
	float minx;
	float dx;
	int x_pixel;
	float miny;
	float dy;
	float minz;
	float dz;
	int z_pixel;
	int lambda_pixel;
	float lambda0;
	float lambda_width_in_A;
	float speedoflight;
};

// Macro variables that should be defined at compile-time: MAX_LAMBDA_PIXEL, MAX_DEPTH
#define DATA_IN_PER_POINT 3
#define DATA_OUT_PER_POINT 4
#define LEAF 3
#define WENT_DOWN_LEFT 1
#define WENT_DOWN_RIGHT 2
#define WENT_DOWN_BOTH 3

kernel void calculate_ray(global const float* coords, global const uchar* nodes, global const float* data_in, constant struct parameters_struct* params, global float* data_out,
	global double* debug_buffer) {
	
	// Reading input
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
	
	float sum = 0;
	/*if (id == 0) {
		for(int i = 1; i < ng + 1; i++) {
			for(int j = 0; j < 3; j++) {
				sum += coords[3*i + j];
				//sum /= 2;
			}
		}
		for(int i = 1; i < ng + 1; i++) {
			sum += nodes[i];
			//sum /= 2;
		}
	}*/
	//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 0;
	
	// Iterate through samples
	uchar stack[MAX_DEPTH]; // Keeps track of path
	float diffs[MAX_DEPTH]; // Stores non-absolute distances to splitting planes
	int single_child_node = 0;
	if (ng%2 == 0) single_child_node = ng/2;
	float best_dist2 = 0;
	for(int k = 0; k < z_pixel; k++) {
		
		pos[2] = convert_float(k)*dz + minz;
		
		// Find closest data point - KD-tree
		
		// Initialize KD-tree NN
		//for(int i = 0; i < MAX_DEPTH; i++) {
		//	stack[i] = 0;
		//	diffs[i] = 0;
		//}
		int stack_pointer = 0;
		int current_node = 1;
		int min_index = -1;
		
		// Continue going down and up the tree until we go up the root, indicating that the search has finished
		if (1) {
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 1] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 2] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 3] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 4] = convert_float(min_index);
			while (current_node > 0) {
				
				// Descend until we reach a leaf
				while (1) {
					
					/*if (id == 11485 && current_node >= ng + 1) {
						// Read too far
						debug_buffer[0] = 1;
						debug_buffer[1] = id;
						debug_buffer[2] = stack[stack_pointer - 1];
						debug_buffer[3] = current_node - 13000000;
						debug_buffer[4] = ng + 1 - 13000000;
						return;
					}*/
					
					// Stop at leaf
					int axis = nodes[current_node];
					if (axis == LEAF)
						break;
					
					/*if (stack_pointer >= MAX_DEPTH) {
						// Read too far
						debug_buffer[0] = 2;
						return;
					}
					if (axis >= 3) {
						// Read too far
						debug_buffer[0] = 3;
						return;
					}
					if (3*current_node + axis >= 3*(ng + 1)) {
						// Read too far
						debug_buffer[0] = 4;
						return;
					}*/
					
					// Continue descending and keep track of path on stack
					diffs[stack_pointer] = pos[axis] - coords[3*current_node + axis];
					if (diffs[stack_pointer] < 0) {
						// Left child
						stack[stack_pointer++] = (current_node == single_child_node ? WENT_DOWN_BOTH : WENT_DOWN_LEFT);
						current_node = current_node << 1;
						/*if (id == 11485 && current_node == 6553600) {
							// Read too far
							debug_buffer[10] = 1;
							debug_buffer[11] = id;
							debug_buffer[12] = stack[stack_pointer - 1];
							debug_buffer[13] = current_node - 6000000;
							debug_buffer[14] = ng + 1 - 13000000;
						}*/
					} else {
						// Right child
						stack[stack_pointer] = WENT_DOWN_RIGHT;
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
					
					/*if (current_node >= ng + 1) {
						// Read too far
						debug_buffer[10] = 5;
						return;
					}*/
					
					// Check if current node is best candidate
					float temp_dist = 0;
					for(int l = 0; l < 3; l++) {
						float temp = pos[l] - coords[3*current_node + l];
						temp_dist += temp*temp;
					}
					if (temp_dist < best_dist2 || min_index == -1) {
						min_index = current_node;
						best_dist2 = temp_dist;
					}
					
					int axis = nodes[current_node];
					/*if (stack_pointer >= MAX_DEPTH && (axis != LEAF || stack_pointer > MAX_DEPTH)) {
						// Read too far
						debug_buffer[0] = stack_pointer;
						debug_buffer[1] = MAX_DEPTH;
						debug_buffer[2] = axis;
						debug_buffer[3] = LEAF;
						debug_buffer[4] = 2 + (axis != LEAF && stack[stack_pointer] != WENT_DOWN_BOTH);
						return;
					}*/
					
					// Consider branching out
					if (axis != LEAF && stack[stack_pointer] != WENT_DOWN_BOTH && diffs[stack_pointer]*diffs[stack_pointer] < best_dist2) {
						// Node is not a leaf, has an unexplored branch and splitting plane is within best distance, so branch out
						if (stack[stack_pointer] == WENT_DOWN_RIGHT) {
							// Left child
							current_node = current_node << 1;
						} else {
							// Right child
							current_node = current_node << 1;
							current_node++;
							/*if (id == 11485 && current_node >= ng + 1) {
								// Read too far
								debug_buffer[15] = 1;
								debug_buffer[16] = id;
								debug_buffer[17] = stack[stack_pointer];
								debug_buffer[18] = current_node - 13000000;
								debug_buffer[19] = ng + 1 - 13000000;
							}*/
						}
						stack[stack_pointer++] = WENT_DOWN_BOTH;
						break;
					} else {
						// Continue going up
						current_node = current_node >> 1;
						stack_pointer--;
					}
					
				}
				
			}
			
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 5] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 6] = convert_float(min_index - 1000000);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 7] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 8] = convert_float(min_index);
			//data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 9] = convert_float(min_index);
			//return;
		}
		
		// End of KD-tree NN-search
		
		/*if (min_index < 1 || min_index >= ng + 1) {
			debug_buffer[0] = 7;
			debug_buffer[1] = debug_buffer[1] + 1;
			debug_buffer[2] = id;
			debug_buffer[3] = z_pixel;
			debug_buffer[4] = min_index;
			return;
		}*/
		
		// Get data from closest point
		int index = DATA_IN_PER_POINT*min_index;//DATA_IN_PER_POINT*min_index;
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
	
}

/*		if (id == 1) {
			while (current_node > 0) {
				
				// Descend until we reach a leaf
				while (1) {
					
					if (current_node >= ng + 1) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 1;
						return;
					}
					
					// Stop at leaf
					int axis = nodes[current_node];
					if (axis == LEAF)
						break;
					
					if (stack_pointer >= MAX_DEPTH) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 2;
						return;
					}
					if (axis >= 3) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 3;
						return;
					}
					if (3*current_node + axis >= 3*(ng + 1)) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 4;
						return;
					}
					
					// Continue descending and keep track of path on stack
					diffs[stack_pointer] = pos[axis] - coords[3*current_node + axis];
					if (diffs[stack_pointer] < 0) {
						// Left child
						current_node = current_node << 1;
						stack[stack_pointer++] = (current_node == single_child_node ? WENT_DOWN_BOTH : WENT_DOWN_LEFT);
					} else {
						// Right child
						stack[stack_pointer] = WENT_DOWN_RIGHT;
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
					
					if (current_node >= ng + 1) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = 5;
						return;
					}
					
					// Check if current node is best candidate
					float temp_dist = 0;
					for(int l = 0; l < 3; l++) {
						float temp = pos[l] - coords[3*current_node + l];
						temp_dist += temp*temp;
					}
					if (temp_dist < best_dist2 || min_index == -1) {
						min_index = current_node;
						best_dist2 = temp_dist;
					}
					
					int axis = nodes[current_node];
					if (stack_pointer >= MAX_DEPTH && (axis != LEAF || stack_pointer > MAX_DEPTH)) {
						// Read too far
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id] = stack_pointer;
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 1] = MAX_DEPTH;
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 2] = axis;
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 4] = LEAF;
						data_out[DATA_OUT_PER_POINT*lambda_pixel*local_id + 5] = 2 + (axis != LEAF && stack[stack_pointer] != WENT_DOWN_BOTH);
						return;
					}
					
					// Consider branching out
					if (axis != LEAF && stack[stack_pointer] != WENT_DOWN_BOTH && diffs[stack_pointer]*diffs[stack_pointer] < best_dist2) {
						// Node is not a leaf, has an unexplored branch and splitting plane is within best distance, so branch out
						if (stack[stack_pointer] == WENT_DOWN_RIGHT) {
							// Left child
							current_node = current_node << 1;
						} else {
							// Right child
							current_node = current_node << 1;
							current_node++;
						}
						stack[stack_pointer++] = WENT_DOWN_BOTH;
						break;
					} else {
						// Continue going up
						current_node = current_node >> 1;
						stack_pointer--;
					}
					
				}
				
			}
		}*/
