#include "../config.h"
#include "header.h"
#include <cstdlib>
#include <cstdio>

cube::cube(const int invars, const int ingrid)
{
	dim=3;
	qtype=empty;
	nvars=invars;
	ng=ingrid;
	grid = new tcoord[dim];
	for (int k=0; k<dim; k++)
		grid[k].reserve(ng);
//	cout << "Maximum size vector" << grid[0].max_size() << "\n Current size" << grid[0].size() << "\n Capacity" << grid[0].capacity() << "\n";
	vars = new tphysvar[nvars];
	for (int k=0; k<nvars; k++)
		vars[k].reserve(ng);
};

cube::~cube()
{};

int cube::readngrid()
{
	return ng;
};

int cube::readnvars()
{
	return nvars;
};

void cube::settype(EqType input)
{
	qtype=input;
};

EqType cube::readtype()
{
	return qtype;
};

void cube::setgrid(tgrid ingrid)
{
	grid=ingrid;
};

tgrid cube::readgrid()
{
	return grid;
};

void cube::fillcube()
{
	switch (qtype) {
		case builtin:
			builtingrid(x_pixel, y_pixel, z_pixel, grid);
			for (int i=0; i<ng; i++)
			{
				double x=grid[0][i];
				double y=grid[1][i];
				double z=grid[2][i];
				double R = length*1000./M_PI;
		                double r = sqrt(pow(sqrt(pow(x,2)+pow(z,2))-R,2)+pow(y,2));
//		                r/=width*1000./0.02; // normalize r to 0 -> 1
		                double phi = atan(y/(sqrt(pow(x,2)+pow(z,2))-R));
		                if (pow(x,2)+pow(z,2)<pow(R,2)) phi+=M_PI; // (x,z) lies inside the circle in the XZ-plane
		                double z_or;
		                if ((x>0.00001)||(x<-0.00001)) // x lies far enough from 0, no division problems
		                {
		                        z_or = atan(z/x)/M_PI;
		                }
		                else
		                {
		                        z_or = .5; // does not matter which value is taken here, because it is the singular point of the coordinate system
		                }
		                if (z_or<0) z_or++;  // atan returns a value between -pi/2 and pi/2
				vars[0][i]=density(r,phi,z_or);
				vars[1][i]=temperature(r,phi,z_or);
			}
			break;
		case empty:
			cout << "Error: first set the type of the equilibrium with datacube.settype(char*)\n";
			exit(EXIT_FAILURE);
			break;
		default:
			cout << "Error: unknown equilibrium type";
			exit(EXIT_FAILURE);
			break;
	}
};
