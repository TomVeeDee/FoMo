#include "../config.h"
#include "header.h"
#include <cstdlib>
#include <cstdio>

tphysvar pow(const double base, tphysvar const & in)
{
	tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=pow(base,in[i]);
	}

	return out;
}

tphysvar operator/(tphysvar const & a, tphysvar const & b)
{
	tphysvar out;
	int s=a.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=a[i]/b[i];
	}

	return out;
}

tphysvar operator*(tphysvar const & a, tphysvar const & b)
{
	tphysvar out;
	int s=a.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=a[i]*b[i];
	}

	return out;
}

tphysvar log10(tphysvar const & in)
{
	tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=log10(in[i]);
	}
	return out;
}

tphysvar sqrt(tphysvar const & in)
{
	tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=sqrt(in[i]);
	}
	return out;
}

tphysvar operator*(double const& c, tphysvar const& in)
{
	tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=c*in[i];
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
	{
		grid[k].resize(ng);
	}
	vars = new tphysvar[nvars];
	for (int k=0; k<nvars; k++)
	{
		vars[k].resize(ng);
	}	
};

cube::~cube()
{};

int cube::readdim() const
{
	return dim;
};

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
			cout << "Setting up built-in equilibrium... " << flush;
			builtingrid(eqx, eqy, eqz, grid);
			for (int i=0; i<nvars; i++)
			{
				vars[i].resize(ng);
			}

#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (int i=0; i<ng; i++)
			{
				double x=grid[0][i];
				double y=grid[1][i];
				double z=grid[2][i];
				// curvature radius in Mm
				double R = length*1000./M_PI;
				// map the torus to a cylinder, using simple toroidal coordinates 
		                double r = sqrt(pow(sqrt(pow(x,2)+pow(z,2))-R,2)+pow(y,2));
		                double phi = atan2(y,(sqrt(pow(x,2)+pow(z,2))-R));
		                double z_or;
		                if ((x>0.00001)||(x<-0.00001)) // x lies far enough from 0, no division problems
		                {
		                        z_or = atan2(z,x)/M_PI;
		                }
		                else
		                {
		                        z_or = .5; // does not matter which value is taken here, because it is the singular point of the coordinate system
		                }
				vars[0][i]=density(r,phi,z_or);
				vars[1][i]=temperature(r,phi,z_or);
				vars[2][i]=0.;
				vars[3][i]=30.*sin(M_PI*z_or);
				vars[4][i]=0.;
			}
			cout << "Done!" << endl << flush;
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
}
