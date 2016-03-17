#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"

/**
 * @brief The default constructor for a RenderCube.
 * 
 * As part of the initialisation, the information in the GoftCube is copied (such as chiantifile, 
 * abundfile, lambda0, grid). The rendermethod defaults to "NearestNeighbour", and the observationtype to 
 * Spectroscopic.\n
 * The resolution is also set to initial values of 101 (x resolution), 102 (y resolution), 300
 * (resolution along LOS). The number of pixels in the wavelength direction is set to 30, and the width
 * of the spectral window to \f$.13\AA{}\f$.
 * @param goftcube The GoftCube from which the RenderCube must be constructed.
 */
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
	rendermethod="NearestNeighbour";
	observationtype=Spectroscopic;
}

/**
 * @brief This reads the current resolution in the RenderCube.
 * 
 * If the resolution has not been set by setresolution() before, the default values of 
 * 101, 102, 300, 30, .13 are returned.
 * @param nx This is the resolution in the x direction (POS).
 * @param ny This is the resolution in the y direction (POS).
 * @param nz This is the resolution along the LOS.
 * @param nlambda This is the number of wavelength pixels.
 * @param lambdawidth This is the width of the spectral window. It is given in \f$\AA{}\f$.
 */
void FoMo::RenderCube::readresolution(int & nx, int & ny, int & nz, int & nlambda, double & lambdawidth)
{
	nx=x_pixel;
	ny=y_pixel;
	nz=z_pixel;
	nlambda=lambda_pixel;
	lambdawidth=lambda_width;
}

/**
 * @brief This sets the resolution of the rendering to be performed.
 * 
 * It should be obvious that these values should be positive.
 * @param nx The resolution in the x direction (POS).
 * @param ny The resolution in the y direction (POS).
 * @param nz The resolution along the LOS.
 * @param nlambda The number of wavelength points.
 * @param lambdawidth The width of the spectral window in \f$\AA{}\f$.
 */
void FoMo::RenderCube::setresolution(const int & nx, const int & ny, const int & nz, const int & nlambda, const double & lambdawidth)
{
	x_pixel=nx;
	y_pixel=ny;
	z_pixel=nz;
	lambda_pixel=nlambda;
	lambda_width=lambdawidth;
}

/**
 * @brief This sets the rendermethod of the RenderCube.
 * @param inrendermethod A string that describes the rendermethod.
 */
void FoMo::RenderCube::setrendermethod(const std::string inrendermethod)
{
	rendermethod=inrendermethod;
}

/**
 * @brief This returns the currently stored rendermethod.
 * 
 * If it has not been set with setrendermethod() before, it returns the default value of "CGAL".
 * @return This returns the currently stored rendermethod.
 */
std::string FoMo::RenderCube::readrendermethod()
{
	return rendermethod;
}

/**
 * @brief This sets the observation type of the RenderCube.
 * @param obtype The value of observation type that should be used.
 */
void FoMo::RenderCube::setobservationtype(FoMo::FoMoObservationType obtype)
{
	observationtype=obtype;
}

/**
 * @brief This reads the observation type of the RenderCube.
 * 
 * If it has not been set, the default value of Spectroscopic is returned. If render() has been used,
 * it returns Spectroscopic if the number of wavelength pixels was greater than 1, otherwise it will
 * be Imaging.
 * @return It returns the observation type.
 */
FoMo::FoMoObservationType FoMo::RenderCube::readobservationtype()
{
	return observationtype;
}

/**
 * @brief This sets the angles of the RenderCube.
 * @param lin The l-angle.
 * @param bin The b-angle.
 */
void FoMo::RenderCube::setangles(const double lin, const double bin)
{
	l=lin;
	b=bin;
}

/**
 * @brief This reads the viewing angles of the RenderCube.
 * @param lout The l-angle.
 * @param bout The b-angle.
 */
void FoMo::RenderCube::readangles(double & lout, double & bout)
{
	lout=l;
	bout=b;
}