#include "header.h"
#include <cmath>

const int eqx=150;
const int eqy=170;
const int eqz=190;

void builtingrid(const int x, const int y, const int z, tgrid grid)
{
// x, y, z are the number of pixels in the respective direction (usually you'd want to use eqx, eqy, eqz)
	double mmperpixel=.72;
	for (int i=0; i<x; i++)
		for (int j=0; j<y; j++)
			for (int k=0; k<z; k++)
			{
				int index=i*(y*z)+j*z+k;
				grid[0].push_back((i-x/2)*mmperpixel*1000.);
				grid[1].push_back((j-y/2)*mmperpixel*1000.);
				grid[2].push_back((k-z/2)*mmperpixel*1000.);
			}
};

double density(const double r, const double phi, const double z)
{
	if (r/0.02<1-thickness/2) return rhoint*(1-alpha*sin(M_PI*z));
	else if (r/0.02>1+thickness/2) return rhoint/contrast*(1-alpha*sin(M_PI*z));
	else return rhoint/2*(1-alpha*sin(M_PI*z))*((1+1/contrast)-(1-1/contrast)*sin(M_PI/thickness*(r/0.02-1)));
};

double temperature(const double r, const double phi, const double z)
{
	return 1000000.;
};

