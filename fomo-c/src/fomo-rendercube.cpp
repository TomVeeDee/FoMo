#include "../config.h"
#include "FoMo.h"

FoMo::RenderCube::RenderCube(FoMo::GoftCube goftcube)
{
	chiantifile=goftcube.readchiantifile();
	abundfile=goftcube.readabundfile();
//	ion=goftcube.readion();
	lambda0=goftcube.readlambda0();
	grid=goftcube.readgrid();
	x_pixel=101;
	y_pixel=102;
	z_pixel=300;
	lambda_pixel=30;
	lambda_width=.13;
	rendermethod="CGAL";
	observationtype=Spectroscopic;
}

void FoMo::RenderCube::readresolution(int nx, int ny, int nz, int nlambda, double lambdawidth)
{
	nx=x_pixel;
	ny=y_pixel;
	nz=z_pixel;
	nlambda=lambda_pixel;
	lambdawidth=lambda_width;
}

void FoMo::RenderCube::setresolution(const int nx, const int ny, const int nz, const int nlambda, const double lambdawidth)
{
	x_pixel=nx;
	y_pixel=ny;
	z_pixel=nz;
	lambda_pixel=nlambda;
	lambda_width=lambdawidth;
}

void FoMo::RenderCube::setrendermethod(const std::string inrendermethod)
{
	rendermethod=inrendermethod;
}

std::string FoMo::RenderCube::readrendermethod()
{
	return rendermethod;
}

void FoMo::RenderCube::setobservationtype(FoMo::FoMoObservationType obtype)
{
	observationtype=obtype;
}

FoMo::FoMoObservationType FoMo::RenderCube::readobservationtype()
{
	return observationtype;

}