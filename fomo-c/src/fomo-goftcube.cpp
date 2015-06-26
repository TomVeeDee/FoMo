#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"

/**
 * @brief The default constructor for the GoftCube object.
 * 
 * A tgrid is created with length indim. The number of grid points and variables is set to 0. The default
 * chiantifile is "../chiantitables/goft_table_fe_12_0194small_abco.dat", and the default abundance
 * file is "sun_coronal.abund" (which is available as internal variable in the code). The wavelength
 * for this spectral line is set accordingly.
 * @param indim The dimension of the GoftCube is set to indim. It defaults to 3.
 */
FoMo::GoftCube::GoftCube(const int indim)
{
	dim=indim;
	grid.resize(indim);
	ng=0;
	nvars=0;
	chiantifile="../chiantitables/goft_table_fe_12_0194small_abco.dat";
	abundfile="/empty";
//	ion="fe_12";
	lambda0=193.509;
}

/**
 * @brief A specialised constructor from an already existing and populated DataCube.
 * 
 * The specialised constructor loads the dimension and corresponding grid from the DataCube. The
 * dimension, number of grid points and grid is set from the information in the argument.
 * The number of variables is set to 0. The default
 * chiantifile is "../chiantitables/goft_table_fe_12_0194small_abco.dat", and the default abundance
 * file is "sun_coronal.abund" (which is available as internal variable in the code). The wavelength
 * for this spectral line is set accordingly.
 * @param incube The datacube from which the GoftCube must be constructed.
 */
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