// writing file with grid points to check inside/outside loop

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
//double xl, xr, yl, yr, zl, zr;
//int dimx, dimy, dimz;

void writegrid(double xl, double xr, int dimx, double yl, double yr, int dimy, double zl, double zr, int dimz) {

	ofstream outfile ("/users/cpa/sgijsen/fomo/version_patrick_nov13/fomo-c/tools/grids/pointgrid.dat");
	
	for (int i=0; i< dimx; i++) {
		for (int j=0; j<dimy; j++) {
			for (int k=0; k<dimz; k++) {
				outfile << xl + (xr - xl) * i/(dimx-1) << "\t" << yl + (yr - yl) *j/(dimy-1) << "\t" << zl + (zr - zl) *k/(dimz-1) << endl;
			}
		}
	}
	outfile.close();
}

