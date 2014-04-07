#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

double co, vao, ct, ro, bo, po, ce, vae, cte, re, be, pe, aa, amp, kafix, kart, wkrt, moaa, meaa, aa1, tu, dt;
int sigi=-1, sigk=1, aa0=1, dimr, diml, dimz, dimt;

void readequilibrium() {

// Read in the equilibrium parameters, the eigenfrequency selected for the kink (general: m=1) mode and the data corresponding to the simulation box (dimensions, resolutions)

string mystr;
ifstream myfile ("/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/alldata.dat");

getline (myfile,mystr);  stringstream(mystr) >> co;
getline (myfile,mystr);  stringstream(mystr) >> vao;
getline (myfile,mystr);  stringstream(mystr) >> ct;
getline (myfile,mystr);  stringstream(mystr) >> ro;
getline (myfile,mystr);  stringstream(mystr) >> bo;
getline (myfile,mystr);  stringstream(mystr) >> po;
getline (myfile,mystr);  stringstream(mystr) >> ce;
getline (myfile,mystr);  stringstream(mystr) >> vae;
getline (myfile,mystr);  stringstream(mystr) >> cte;
getline (myfile,mystr);  stringstream(mystr) >> re;
getline (myfile,mystr);  stringstream(mystr) >> be;
getline (myfile,mystr);  stringstream(mystr) >> pe;
getline (myfile,mystr);  stringstream(mystr) >> aa;
getline (myfile,mystr);  stringstream(mystr) >> amp;
getline (myfile,mystr);  stringstream(mystr) >> kafix;
getline (myfile,mystr);  stringstream(mystr) >> kart;
getline (myfile,mystr);  stringstream(mystr) >> wkrt;
getline (myfile,mystr);  stringstream(mystr) >> moaa;
getline (myfile,mystr);  stringstream(mystr) >> meaa;
getline (myfile,mystr);  stringstream(mystr) >> aa1;
getline (myfile,mystr);  stringstream(mystr) >> tu;
getline (myfile,mystr);  stringstream(mystr) >> dt;
getline (myfile,mystr);  stringstream(mystr) >> dimr;
getline (myfile,mystr);  stringstream(mystr) >> diml;
getline (myfile,mystr);  stringstream(mystr) >> dimz;
getline (myfile,mystr);  stringstream(mystr) >> dimt;

}

int main() {

	readequilibrium();

	double xr = 2*aa;  double xl = -xr;
	int dimx = dimr;
	double yl = -xr;  double yr = xr;
	int dimy = dimr;
	double zl = 0;  double zr = M_PI*aa/(2*kart);  		// Half the loop length

	double * advectedpos = new double [dimx][dimy][dimz][3];
	signed int * loc = new signed int [dimx][dimy][dimz];


// Construct initial state (t=0) of positions and loc indicators

	for (int i=0; i< dimx; i++) {
		for (int j=0; j<dimy; j++) {
			for (int k=0; k<dimz; k++) {
				double xpos = xl + (xr - xl) * i/(dimx-1);
				double ypos = yl + (yr - yl) * j/(dimy-1);
				double zpos = zl + (zr - zl) * k/(dimz-1);
				advectedpos[i][j][k][0] = xpos;
				advectedpos[i][j][k][1] = ypos;
				advectedpos[i][j][k][2] = zpos;
				if (xpos*xpos+ypos*ypos < aa*aa) loc[i][j][k]=1; else loc[i][j][k]=-1;
			}
		}
	}

	cout << " position: " << advectedpos[43][46][48][0] << ", " << advectedpos[43][46][48][1] << ", " << advectedpos[43][46][48][2] << endl;

return 0;
}
