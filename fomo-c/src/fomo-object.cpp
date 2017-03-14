#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <map>
#include <iostream>

/**
 * @brief This member is the constructor of the FoMoObject.
 * 
 * The default constructor for the FoMo::FoMoObject constructs in order 
 * FoMo::FoMoObject.datacube(indim), FoMo::FoMoObject.goftcube(datacube), 
 * FoMoObject.rendering(goftcube).
 * 
 * @param indim The integer indim sets the dimension of the datacube. It defaults to 3.
 */
FoMo::FoMoObject::FoMoObject(const int indim):
	datacube(indim), goftcube(datacube), rendering(goftcube)
{
}

/**
 * @brief The default destructor for the FoMoObject. 
 */
FoMo::FoMoObject::~FoMoObject()
{
}

/**
 * @brief This member sets the rendermethod to be used.
 * 
 * Use this method to set the rendermethod. At the moment (version 3.3), there 
 * are three rendermethods: "CGAL", "CGAL2D" and "NearestNeighbour".
 * It should be read before the render() is called, because that used the information here.
 * @param inrendermethod The function takes a string as an argument, which is 
 * then internally connected to a rendermethod.
 */
void FoMo::FoMoObject::setrendermethod(const std::string inrendermethod)
{
	rendering.setrendermethod(inrendermethod);
}

/**
 * @brief This reads the rendermethod that is currently set for the FoMoObject.

The function returns the string for the rendermethod that is currently set in the FoMoObject.  
 * @return The return string contains the rendermethod.
 */
std::string FoMo::FoMoObject::readrendermethod()
{
	return rendering.readrendermethod();
}

/**
 * @brief This routine returns the GoftCube of the FoMoObject.
 * 
 * FoMoObject.goftcube is a private member, thus its access can only be done through this function.
 * @return The return value is the goftcube of the FoMoObject and is of type GoftCube. It contains the calculated
 * emission at each datapoint of the original mesh.
 */
FoMo::GoftCube FoMo::FoMoObject::readgoftcube()
{
	 return goftcube;
}

/**
 * @brief This sets the abundance file to be used for the CHIANTI conversion to emissivity.
 * 
 * The abundance file is set to the string argument. The string argument should 
 * point to a valid abundance file: these are the CHIANTI abundance files (such as sun_coronal_2012_schmelz.abund).
 * The latter is available internally in the code, and is the standard used by the code. 
 * @param inabund A string with the path to a CHIANTI abundance file.
 */
void FoMo::FoMoObject::setabundfile(const std::string inabund)
{
	goftcube.setabundfile(inabund);
	rendering.setabundfile(inabund);
}

/**
 * @brief This function returns the abundance file currently used.
 * 
 * If the abundance file has not been set previously, the default value "/empty" is returned. In this case,
 * it means that the built-in sun_coronal_2012_schmelz.abund is being used. 
 * @return This returns the path to the abundance file currently in use. 
 */
std::string FoMo::FoMoObject::readabundfile()
{
	return goftcube.readabundfile();
}

/**
 * @brief This routine sets the chiantifile.
 * 
 * The chiantifile contains the tabulated G(T) from CHIANTI, to be used 
 * for converting \f$\rho\f$ and T to emissivity \f$\rho^2G(\rho,T)\f$.
 * @param inchianti The input string should contain a path to a valid chiantifile. 
 * These are the tabulated intensities included in the FoMo archive, under the directory chiantifiles.
 */
void FoMo::FoMoObject::setchiantifile(const std::string inchianti)
{
	goftcube.setlambda0(readgoftfromchianti(inchianti));
	goftcube.setchiantifile(inchianti);
	rendering.setlambda0(readgoftfromchianti(inchianti));
	rendering.setchiantifile(inchianti);
}

/**
 * @brief This returns the currently used chiantifile.
 * 
 * If chiantifile has not been set before with setchiantifile, then the default 
 * value "../chiantitables/goft_table_fe_12_0194_abco.dat" is returned.
 * @return The return string is the path to the currently used chiantifile.
 */
std::string FoMo::FoMoObject::readchiantifile()
{
	return goftcube.readchiantifile();
}

/**
 * @brief This routines gets additional data into the datacube.
 * 
 * The FoMoObject adds the argument data to its protected member datacube.
 * @param coordinate This is a vector of length dimension of the FoMoObject.
 * @param variables This is a vector of length nvars (i.e. number of variables required by rendermethod, e.g. \f$\rho\f$, T, vx, vy, vz for CHIANTI renderings)
 */
void FoMo::FoMoObject::push_back_datapoint(std::vector<double> coordinate, std::vector<double> variables)
{
	this->datacube.push_back(coordinate,variables);
}

/**
 * @brief This stores the grid and data into the FoMoObject.
 * 
 * The protected member datacube is filled with the argument grid and data. tgrid and tvars should have columns of
 * equal length. It sets the dimension, number of gridpoints and number of variables in the FoMoObject.
 * @param ingrid A tgrid of the desired dimension, containing a number of data points.
 * @param indata The physical variables at the corresponding grid points.
 */
void FoMo::FoMoObject::setdata(tgrid& ingrid, tvars& indata)
{
	this->datacube.setdata(ingrid,indata);
}

/**
 * @brief This sets the observation type.
 * @param observationtype The observation type of the FoMoObject is set to the argument.
 */
void FoMo::FoMoObject::setobservationtype(FoMoObservationType observationtype)
{
	this->rendering.setobservationtype(observationtype);
}

/**
 * @brief This reads the currently set observation type.
 * 
 * If the observation type has not been set, it returns the default Spectroscopic.
 * @return It returns the FoMoObservationType. 
 */
FoMo::FoMoObservationType FoMo::FoMoObject::readobservationtype()
{
	return this->rendering.readobservationtype();
}

/**
 * @brief This sets the resolution of the rendering.
 * 
 * It should be done before render() is called, because the resolution is read in that routine.
 * @param x_pixel This is the resolution of the x-direction of the rendercube.
 * @param y_pixel This is the resolution in the y-direction of the rendercube.
 * @param z_pixel This specifies how many points there are along the LOS. It should 
 * probably be around the finest resolution of the datacube grid, to achieve the optimum rendering.
 * @param lambda_pixel This specifies how many wavelength points are required.
 * @param lambda_width This sets the width of the wavelength window. It is in \f$m/s\f$.
 */
void FoMo::FoMoObject::setresolution(const int & x_pixel, const int & y_pixel, const int & z_pixel, const int & lambda_pixel, const double & lambda_width)
{
	this->rendering.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
}

/**
 * @brief This reads the current resolution of the rendering.
 * @param x_pixel The number of pixels in the x-direction.
 * @param y_pixel The number of pixels in the y-direction.
 * @param z_pixel The number of pixels along the LOS.
 * @param lambda_pixel The number of wavelength pixels.
 * @param lambda_width The width of the wavelength window.
 */
void FoMo::FoMoObject::readresolution(int & x_pixel, int & y_pixel, int & z_pixel, int & lambda_pixel, double & lambda_width)
{
	this->rendering.readresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
}


/**
 * @brief This sets the grid and data of the rendering.
 * 
 * @param ingrid This is a parameter of type tgrid, and contains the grid of the rendering. Usually it is x,y,lambda for a spectroscopic rendering.
 * @param invars This sets the variables for the rendering and is of type tvars. Usually it would be intensity.
 */
void FoMo::FoMoObject::setrenderingdata(tgrid ingrid, tvars invars)
{
	rendering.setdata(ingrid,invars);
}

/**
 * @brief This routine returns the rendering of the FoMoObject.
 * 
 * At the moment, it is equivalent to accessing FoMoObject.rendering, because that member is public.
 * @return The return value is the rendering of the FoMoObject and is of type RenderCube. It contains the forward model.
 */
FoMo::RenderCube FoMo::FoMoObject::readrendering()
{
	return rendering;
}

/**
 * @brief This routine sets the base of the file where the output is written.
 * 
 * The outfile is the base of the file name which will be used for writing the rendering to. It could be 
 * e.g. "./fomo-output" and the default value is "". When writing, the render() routine will write the 
 * results to outfile, appended with "l"+l+"b"+b, where l and b are the respective viewing angles in degrees.
 * @param instring The base path and filename of the output files into which the rendering will be written.
 */
void FoMo::FoMoObject::setoutfile(const std::string instring)
{
	outfile=instring;
}

/** 
 * I think it would be great if we could document these in the manual. I'm not sure how to do  that.
 */
/// [RenderMethods]
enum FoMoRenderValue
{
	RenderMethodNotDefined,
#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H
	CGAL,
	CGAL2D,
#endif
	NearestNeighbour,
	Projection,
	// add more methods here
	LastVirtualRenderMethod
};
/// [RenderMethods]

static const std::map<std::string, FoMoRenderValue>::value_type RenderMapEntries[]=
{
	/// [Rendermethods]
#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H
	std::map<std::string, FoMoRenderValue>::value_type("CGAL",CGAL),
	std::map<std::string, FoMoRenderValue>::value_type("CGAL2D",CGAL2D),
#endif
	std::map<std::string, FoMoRenderValue>::value_type("NearestNeighbour",NearestNeighbour),
	std::map<std::string, FoMoRenderValue>::value_type("Projection",Projection),
	/// [Rendermethods]
	std::map<std::string, FoMoRenderValue>::value_type("ThisIsNotARealRenderMethod",LastVirtualRenderMethod)
};

// static int nmethods=2;
// if the modular adding is too clumsy, just introduce the above variable.

static std::map<std::string, FoMoRenderValue> RenderMap{ &RenderMapEntries[0], &RenderMapEntries[LastVirtualRenderMethod-1] };

/**
 * @brief This is the main render routine for the FoMo::FoMoObject.
 * 
 * The routine starts from the FoMo::FoMoObject.datacube, and renders it using the 
 * rendermethod (set with FoMo::FoMoObject.setrendermethod), using the resolution set in 
 * FoMoObject.rendering.setresolution. It does this for each possible combination of angles
 * in lvec and bvec (and so in total n*m renderings will be performed for lvec of length n and 
 bvec of length m).\n 
  * The routine sets the FoMoObservationType to Spectroscopic if the number of
  * wavelength pixels (lambda_pixel) is larger than 1, otherwise it is Imaging. \n
   The routine writes out all the renderings done (i.e. each combination of lvec and bvec) to a file 
 * outfile (set with setoutfile()) appended with "l"+l+"b"+b (where l and b are the viewing angles 
 * converted into 
 * degrees). When the routine returns, the last rendering is stored (with viewing angles lvec.back() 
 * and bvec.back()).\n
 * If the render is done with CHIANTI, the called rendermethods convert the \f$\rho\f$ 
 * and T in the FoMoObject.datacube to emissivity 
 * (proportional to \f$\rho^2G(\rho,T)\f$) using the chiantifile (set with setchiantifile()) and abundfile
 * (set with setabundfile()). 
 * @param lvec This is vector of l-angles which need to be considered by the rendering. The values should be in radians.
 * @param bvec This is vector of b-angles which need to be considered by the rendering. The values should be in radians.
 */
void FoMo::FoMoObject::render(const std::vector<double> lvec, const std::vector<double> bvec)
{
	FoMo::GoftCube tmpgoft;
	FoMo::RenderCube tmprender(this->goftcube);
	int x_pixel, y_pixel, z_pixel, lambda_pixel;
	double lambda_width;
	this->rendering.readresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	
	tmpgoft=FoMo::emissionfromdatacube(this->datacube,this->rendering.readchiantifile(),this->rendering.readabundfile(),this->rendering.readobservationtype());
	this->goftcube=tmpgoft;
	
	switch (RenderMap[rendering.readrendermethod()])
	{
		// add other rendermethods here
#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H
		case CGAL2D:
			std::cout << "Using CGAL-2D for rendering." << std::endl << std::flush;
			if (bvec.size()>0) std::cout << "Warning: the bvec-values are not used in this 2D routine." << std::endl << std::flush;
			tmprender=FoMo::RenderWithCGAL2D(this->datacube,this->goftcube,this->rendering.readobservationtype(),
			x_pixel, y_pixel, lambda_pixel, lambda_width, lvec, this->outfile);
			break;
		case CGAL:
			std::cout << "Using CGAL for rendering." << std::endl << std::flush;
			tmprender=FoMo::RenderWithCGAL(this->datacube,this->goftcube,this->rendering.readobservationtype(),
			x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width,lvec,bvec,this->outfile);
			break;
#endif
		case NearestNeighbour:
			std::cout << "Using nearest-neighbour rendering." << std::endl << std::flush;
			tmprender=FoMo::RenderWithNearestNeighbour(this->goftcube,x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width, lvec, bvec, this->outfile);
			break;
		case Projection:
			std::cout << "Using projection for rendering." << std::endl << std::flush;
			tmprender=FoMo::RenderWithProjection(this->goftcube,x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width,lvec,bvec, this->outfile);
			break;
		case LastVirtualRenderMethod: // this should not be reached, since it is excluded from the map
		default:
			std::cerr << "Error: unknown rendering method." << std::endl << std::flush;
			exit(EXIT_FAILURE);
			break;
	}

	tmprender.setrendermethod(rendering.readrendermethod());
	tmprender.setobservationtype(rendering.readobservationtype());
	this->rendering=tmprender;
}

/**
 * @brief A single rendering of a datacube.
 * 
 * For practical reasons, it results in a call to render(std::vector<double>, std::vector<double>),
 * with as first argument {l} and second argument {b}. First, the simulation box is rotated with an angle -b around the y-axis, then it is rotated around an angle -l around the z-axis.
 * @param b First, the FoMoObject.goftcube is rotated an angle -b around the y-axis. It defaults to 0.
 * @param l Then, the FoMoObject.goftcube is rotated an angle -l around the z-axis. It defaults to 0.
 */
void FoMo::FoMoObject::render(const double l, const double b)
{
	std::vector<double> lvec{l};
	std::vector<double> bvec{b};
	this->render(lvec,bvec);
}
