#include "../config.h"
#include "FoMo.h"
#include <cmath>

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

FoMo::tphysvar FoMo::operator*(FoMo::tphysvar const & a, FoMo::tphysvar const & b)
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
		out[i]=std::log10(in[i]);
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
