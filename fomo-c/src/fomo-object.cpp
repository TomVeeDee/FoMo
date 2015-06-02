#include "FoMo.h"
#include <map>
#include <iostream>

FoMo::FoMoObject::FoMoObject()
{
}

FoMo::FoMoObject::~FoMoObject()
{
}

void FoMo::FoMoObject::setrendermethod(const std::string inrendermethod)
{
	rendermethod=inrendermethod;
}

std::string FoMo::FoMoObject::readrendermethod()
{
	return rendermethod;
}

enum FoMoRenderValue
{
	RenderMethodNotDefined,
	CGAL,
	// add more methods between these
	LastVirtualRenderMethod
};

typedef std::map<std::string, FoMoRenderValue> StringMap;

typedef StringMap::value_type     StringMapValue;

static const std::map<std::string, FoMoRenderValue>::value_type RenderMapEntries[]=
{
	std::map<std::string, FoMoRenderValue>::value_type("CGAL",CGAL),
};

static std::map<std::string, FoMoRenderValue> RenderMap( &RenderMapEntries[CGAL], &RenderMapEntries[LastVirtualRenderMethod] );

void FoMo::FoMoObject::render(const std::vector<double> lvec, const std::vector<double> bvec)
{
	FoMo::DataCube tmprendering;
	switch ( RenderMap[rendermethod] )
	{
		// add other rendermethods here
		case CGAL:
			tmprendering=FoMo::RenderWithCGAL(this->datacube,lvec,bvec);
			break;
		default:
			std::cerr << "Error: unknown rendering method." << std::endl << std::flush;
			exit(EXIT_FAILURE);
			break;
	}
}

void FoMo::FoMoObject::render(const double l, const double b)
{
	std::vector<double> lvec{l};
	std::vector<double> bvec{b};
	this->render(lvec,bvec);

}