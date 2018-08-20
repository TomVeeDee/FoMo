struct __attribute__ ((packed)) parameters_struct {
	int ng;
	float minx;
	float dx;
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

#define MAX_LAMBDA_PIXEL 200 // TODO: Try to make this variable by changing the kernel at compile-time (so run-time for the host)?
#define DATA_IN_PER_POINT 3
#define DATA_OUT_PER_POINT 4

kernel void calculate_ray(global const float* coords, global const float* nodes, global const float* data_in, constant struct parameters_struct* params, global float* data_out) {
	
	// Reading input
	size_t i = get_global_id(0); // X-coordinate of current pixel
	size_t j = get_global_id(1); // Y-coordinate of current pixel
	size_t x_pixel = get_global_size(0);
	int ng = params->ng;
	float x = convert_float(i)*params->dx + params->minx;
	float y = convert_float(j)*params->dy + params->miny;
	int lambda_pixel = params->lambda_pixel;
	float lambda0 = params->lambda0;
	float lambda1 = params->lambda_width_in_A/(lambda_pixel - 1);
	float lambda2 = params->lambda_width_in_A/2.;
	float lambda3 = lambda0/params->speedoflight;
	
	// Initialize intensity vector
	float intens[MAX_LAMBDA_PIXEL];
	float lambdaval[MAX_LAMBDA_PIXEL];
	for(int il = 0; il < lambda_pixel; il++) {
		intens[il] = 0;
		lambdaval[il] = convert_float(il)*lambda1 - lambda2;
	}
	
	// Iterating through samples
	for(int k = 0; k < params->z_pixel; k++) {
		
		float z = convert_float(k)*params->dz + params->minz;
		
		// Find closest data point
		int min_index;
		float min_dist2;
		int index = 0;
		for(int l = 0; l < 100; l++) {
			float diffx = x - coords[index++];
			float diffy = y - coords[index++];
			float diffz = z - coords[index++];
			float dist2 = diffx*diffx + diffy*diffy + diffz*diffz;
			if (dist2 < min_dist2 || l == 0) {
				min_index = l;
				min_dist2 = dist2;
			}
		}
		
		// Get data from closest point
		index = DATA_IN_PER_POINT*min_index;
		float intpolpeak = data_in[index++];
		float intpolfwhm = data_in[index++];
		float intpollosvel = data_in[index];
		
		// Calculate intensities
		if (lambda_pixel > 1) {
			// Spectroscopic study
			for (int il = 0; il < 1; il++) {
				// lambda the relative wavelength around lambda0, with a width of lambda_width
				// if intpolpeak is not zero then the correct expression is used. otherwise, the intensity is just 0
				// it remains to be tested if this is faster than just the direct computation
				float temp = (lambdaval[il] - intpollosvel/lambda3)/intpolfwhm;
				intens[il] += intpolpeak*exp(-temp*temp*2.772588722); // Last number is 4*ln(2)
			}
		} else {
			// Imaging study
			intens[0] += intpolpeak;
		}
	}
	
	// Write out intensity for every wavelength
	int index = DATA_OUT_PER_POINT*(j*x_pixel + i)*lambda_pixel;
	for (int il = 0; il < lambda_pixel; il++) {
		data_out[index++] = x;
		data_out[index++] = y;
		data_out[index++] = lambdaval[il] + lambda0;
		data_out[index++] = intens[il];
	}
	
}
