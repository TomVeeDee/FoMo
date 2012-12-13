#include "../config.h"
#include "header.h"
#include <cstdlib>
#include <cstdio>

tphysvar log10(tphysvar const & in)
{
	tphysvar out;
	tphysvar::const_iterator init=in.begin();
		cout << "a" << *init;
	for (; init != in.end(); ++init)
	{
		cout << "a" << *init;
		out.push_back(log10(*init));
	}
	return out;
}

tphysvar sqrt(tphysvar const & in)
{
	tphysvar out;
	tphysvar::const_iterator init=in.begin();
	for (; init != in.end(); ++init)
	{
		out.push_back(sqrt(*init));
	}
	return out;
}

tphysvar operator*(double const& c, tphysvar const& in)
{
	tphysvar out;
	tphysvar::const_iterator init=in.begin();
	for (; init != in.end(); ++init)
	{
		out.push_back(c*(*init));
	}
	return out;
}

cube::cube(const int invars, const int ingrid, const int indim )
{
	dim=indim;
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

int cube::readngrid() const
{
	return ng;
};

int cube::readnvars() const
{
	return nvars;
};

void cube::settype(EqType input)
{
	qtype=input;
};

EqType cube::readtype() const
{
	return qtype;
};

void cube::setgrid(tgrid ingrid)
{
	grid=ingrid;
};

tgrid cube::readgrid() const
{
	return grid;
};

void cube::setvar(const int nvar, const tphysvar var)
{
	vars[nvar]=var;
}

tphysvar cube::readvar(const int nvar) const
{
	tphysvar var=vars[nvar];
	return var;
}

void cube::fillcube()
{
	switch (qtype) {
		case builtineq:
			builtingrid(eqx, eqy, eqz, grid);
			for (int i=0; i<ng; i++)
			{
				double xacc=grid[0][i];
				double yacc=grid[1][i];
				double zacc=grid[2][i];
				const double psi = 0;
		                double x = (cos(psi)*cos(l)-sin(psi)*sin(l)*sin(b))*xacc+(-sin(psi)*cos(b))*yacc+(-cos(psi)*sin(l)-sin(psi)*cos(l)*sin(b))*zacc;
		                double y = (sin(psi)*cos(l)+cos(psi)*sin(l)*sin(b))*xacc+(cos(psi)*cos(b))*yacc+(-sin(psi)*sin(l)+cos(psi)*cos(l)*sin(b))*zacc;
		                double z = (sin(l)*cos(b))*xacc+(-sin(b))*yacc+(cos(l)*cos(b))*zacc;
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
			cout << "Error: first set the type of the equilibrium with datacube.settype(EqType)\n";
			exit(EXIT_FAILURE);
			break;
		default:
			cout << "Error: unknown equilibrium type";
			exit(EXIT_FAILURE);
			break;
	}
};
