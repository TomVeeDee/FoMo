#include "header.h"
#include <cmath>

void builtingrid(const int x_pixel, const int y_pixel, const int z_pixel, tgrid grid)
{
	int ng=x_pixel*y_pixel*z_pixel;
	for (int i=0; i<x_pixel; i++)
		for (int j=0; j<y_pixel; j++)
			for (int k=0; k<z_pixel; k++)
			{
				int index=i*(y_pixel*z_pixel)+j*z_pixel+k;
				grid[0][index]=i;
				grid[1][index]=j;
				grid[2][index]=k;
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

