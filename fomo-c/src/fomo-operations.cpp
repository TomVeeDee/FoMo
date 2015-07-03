#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <cmath>
#include <cassert>
#include <iostream>

FoMo::tphysvar FoMo::pow(const double base, FoMo::tphysvar const & in)
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

FoMo::tphysvar FoMo::operator/(FoMo::tphysvar const & a, FoMo::tphysvar const & b)
{
	FoMo::tphysvar out;
	unsigned int s=a.size();
	assert(b.size() == s);
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int i=0; i<s; i++)
	{
		out[i]=a[i]/b[i];
	}

	return out;
}

template<typename T> FoMo::tphysvar FoMo::operator*=(FoMo::tphysvar const & a, T b)
{
	return a*b;
}

FoMo::tphysvar FoMo::operator*(FoMo::tphysvar const & a, FoMo::tphysvar const & b)
{
	FoMo::tphysvar out;
	unsigned int s=a.size();
	assert(b.size() == s);
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (unsigned int i=0; i<s; i++)
	{
		out[i]=a[i]*b[i];
	}

	return out;
}

FoMo::tphysvar FoMo::log10(FoMo::tphysvar const & in)
{
	FoMo::tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		if (in[i] <= 0)
		{
			std::cerr << "Warning: some densities were <= 0 at position " << i << std::endl;
			out[i]=0;
		}
		else
		{
			out[i]=std::log10(in[i]);
		}
	}
	return out;
}

FoMo::tphysvar FoMo::sqrt(FoMo::tphysvar const & in)
{
	FoMo::tphysvar out;
	int s=in.size();
	out.resize(s);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<s; i++)
	{
		out[i]=std::sqrt(in[i]);
	}
	return out;
}

FoMo::tphysvar FoMo::operator*(double const& c, FoMo::tphysvar const& in)
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
