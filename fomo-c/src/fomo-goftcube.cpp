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

/**
 * @brief This displays the currently set abundance file.
 * 
 * In case it has not been set before by setabundfile(), then the default value of "/empty" is 
 * returned. In this case, the built-in sun_coronal.abund is being used.
 * @return It returns the path and filename to the abundance file.
 */
std::string FoMo::GoftCube::readabundfile()
{
	return abundfile;
}

/**
 * @brief This returns the chiantifile in use.
 * 
 * In case the chiantifile has not been set with setchiantifile(), the default value
 * "../chiantitables/goft_table_fe_12_0194small_abco.dat" will be returned.
 * @return The path to the chiantifile being used.
 */
std::string FoMo::GoftCube::readchiantifile()
{
	return chiantifile;
}

/**
 * @brief Use this routine to set the abundance file.
 * @param inabund A string containing the path and filename of a valid abundance file 
 * in CHIANTI format.
 */
void FoMo::GoftCube::setabundfile(const std::string inabund)
{
	abundfile=inabund;
}

/**
 * @brief This routine sets the chiantifile.
 * @param inchianti A string containing the path and filename of a valid chiantifile. These 
 * are the tabulated G(T) values to be found in the FoMo repository under chiantitables.
 */
void FoMo::GoftCube::setchiantifile(const std::string inchianti)
{
	chiantifile=inchianti;
}

/**
 * @brief With this member, the wavelength \f$\lambda_0\f$ is set.
 * @param lambda A double value of the wavelength (in \f$\AA{}\f$).
 */
void FoMo::GoftCube::setlambda0(const double lambda)
{
	lambda0=lambda;
}

/**
 * @brief This reads the currently set wavelength.
 * 
 * If the wavelength has not been set, the default value of 193.509 is returned.
 * @return The return value is the wavelength in \f$\AA{}\f$.
 */
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