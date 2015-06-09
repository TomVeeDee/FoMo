#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"

FoMo::GoftCube::GoftCube(const int indim)
{
	dim=indim;
	ng=0;
	nvars=0;
	chiantifile="../chiantitables/goft_table_fe_12_0194small_abco.dat";
	abundfile="/empty";
//	ion="fe_12";
	lambda0=193.509;
}

FoMo::GoftCube::GoftCube(FoMo::DataCube incube)
{
	dim=incube.readdim();
	ng=incube.readngrid();
	grid=incube.readgrid();
	nvars=0;
	chiantifile="../chiantitables/goft_table_fe_12_0194small_abco.dat";
	abundfile="/empty";
//	ion="fe_12";
	lambda0=193.509;
}

std::string FoMo::GoftCube::readabundfile()
{
	return abundfile;
}

std::string FoMo::GoftCube::readchiantifile()
{
	return chiantifile;
}

void FoMo::GoftCube::setabundfile(const std::string inabund)
{
	abundfile=inabund;
}

void FoMo::GoftCube::setchiantifile(const std::string inchianti)
{
	chiantifile=inchianti;
}

void FoMo::GoftCube::setlambda0(const double lambda)
{
	lambda0=lambda;
}

double FoMo::GoftCube::readlambda0()
{
	return lambda0;
}

/*std::string FoMo::GoftCube::readion()
{
	return ion;
}

void FoMo::GoftCube::setion(const std::string inion)
{
	ion=inion;
}*/