#include "../config.h"
#include "FoMo.h"
#include <map>
#include <iostream>

FoMo::FoMoObject::FoMoObject(const int indim):
	datacube(indim), goftcube(datacube), rendering(goftcube)
{
}

FoMo::FoMoObject::~FoMoObject()
{
}

void FoMo::FoMoObject::setrendermethod(const std::string inrendermethod)
{
	rendering.setrendermethod(inrendermethod);
}

std::string FoMo::FoMoObject::readrendermethod()
{
	return rendering.readrendermethod();
}

void FoMo::FoMoObject::setabundfile(const std::string inabund)
{
	goftcube.setabundfile(inabund);
	rendering.setabundfile(inabund);
}

std::string FoMo::FoMoObject::readabundfile()
{
	return goftcube.readabundfile();
}

void FoMo::FoMoObject::setchiantifile(const std::string inchianti)
{
	goftcube.setlambda0(readgoftfromchianti(inchianti));
	goftcube.setchiantifile(inchianti);
	rendering.setlambda0(readgoftfromchianti(inchianti));
	rendering.setchiantifile(inchianti);
}

std::string FoMo::FoMoObject::readchiantifile()
{
	return goftcube.readchiantifile();
}

void FoMo::FoMoObject::setrenderingdata(tgrid ingrid, tvars invars)
{
	rendering.setdata(ingrid,invars);
}

enum FoMoRenderValue
{
	RenderMethodNotDefined,
	CGAL,
	// add more methods between these
	LastVirtualRenderMethod
};

static const std::map<std::string, FoMoRenderValue>::value_type RenderMapEntries[]=
{
	std::map<std::string, FoMoRenderValue>::value_type("CGAL",CGAL),
	std::map<std::string, FoMoRenderValue>::value_type("ThisIsNotARealRenderMethod",LastVirtualRenderMethod)
};

// static int nmethods=2;
// if the modular adding is too clumsy, just introduce the above variable.

static std::map<std::string, FoMoRenderValue> RenderMap{ &RenderMapEntries[0], &RenderMapEntries[LastVirtualRenderMethod-1] };

void FoMo::FoMoObject::render(const std::vector<double> lvec, const std::vector<double> bvec)
{
	FoMo::RenderCube tmprender(this->goftcube);
	int x_pixel, y_pixel, z_pixel, lambda_pixel;
	double lambda_width;
	this->rendering.readresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	
	switch (RenderMap[rendering.readrendermethod()])
	{
		// add other rendermethods here
		case CGAL:
			std::cout << "Using CGAL for rendering." << this;
			//tmprender=FoMo::RenderWithCGAL(datacube,lvec,bvec);
			tmprender=FoMo::RenderWithCGAL(this->datacube,this->rendering.readchiantifile(),this->rendering.readabundfile(),this->rendering.readobservationtype(), 
			x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width,
			lvec,bvec);
			break;
		case LastVirtualRenderMethod: // this should not be reached, since it is excluded from the map
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