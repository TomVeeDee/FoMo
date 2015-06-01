#include "../config.h"
#include "FoMo.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

FoMo::tphysvar pow(const double base, FoMo::tphysvar const & in)
{
	FoMo::tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=std::pow(base,in[i]);
	}

	return out;
}

FoMo::tphysvar operator/(FoMo::tphysvar const & a, FoMo::tphysvar const & b)
{
	FoMo::tphysvar out;
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

FoMo::tphysvar operator*(FoMo::tphysvar const & a, FoMo::tphysvar const & b)
{
	FoMo::tphysvar out;
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

FoMo::tphysvar log10(FoMo::tphysvar const & in)
{
	FoMo::tphysvar out;
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

FoMo::tphysvar sqrt(FoMo::tphysvar const & in)
{
	FoMo::tphysvar out;
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

FoMo::tphysvar operator*(double const& c, FoMo::tphysvar const& in)
{
	FoMo::tphysvar out;
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

FoMo::cube::cube(const int invars, const int ingrid, const int indim )
{
	dim=indim;
	nvars=invars;
	ng=ingrid;
	grid = new tcoord[dim];
	for (int k=0; k<dim; k++)
	{
		grid[k].reserve(ng);
	}
	vars = new tphysvar[nvars];
	for (int k=0; k<nvars; k++)
	{
		vars[k].reserve(ng);
	}	
};

FoMo::cube::~cube()
{};

int FoMo::cube::readdim() const
{
	return dim;
};

int FoMo::cube::readngrid() const
{
	return ng;
};

int FoMo::cube::readnvars() const
{
	return nvars;
};

void FoMo::cube::setgrid(FoMo::tgrid ingrid)
{
	grid=ingrid;
};

FoMo::tgrid FoMo::cube::readgrid() const
{
	return grid;
};

void FoMo::cube::setvar(const int nvar, const FoMo::tphysvar var)
{
	vars[nvar]=var;
}

FoMo::tphysvar FoMo::cube::readvar(const int nvar) const
{
	FoMo::tphysvar var=vars[nvar];
	return var;
}
