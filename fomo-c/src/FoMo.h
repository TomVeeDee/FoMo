#include <vector>
#include <string>
#ifndef FOMO_H
#define FOMO_H 
/**
 * @file 
 * This file contains the specifications for the FoMo code, which should 
 * be included in codes linking against FoMo.
 */

/**
 * This is the namespace that the FoMo code and objects live in.
 */
namespace FoMo
{
	/**
	 * The type tcoord is defined as a vector of doubles. It has length DataCube.ng.
	 */
	typedef std::vector<double> tcoord;
	/**
	 * The type tgrid is the type of a grid. It is a vector of ::tcoord. It has length DataCube.dim.
	 */
	typedef std::vector<tcoord> tgrid;
	/**
	 * The type tphysvar stores the physical variables in the datapoints in a ::tcoord. It has length DataCube.ng.
	 */
	typedef std::vector<double> tphysvar;
	/**
	 * The type tvars is a vector of ::tphysvar. It has length DataCube.nvars.
	 */
	typedef std::vector<tphysvar> tvars;
	
	tphysvar pow(const double, tphysvar const&);
	tphysvar operator/(tphysvar const&, tphysvar const&);
	tphysvar operator*(tphysvar const&, tphysvar const&);
	tphysvar operator*=(tphysvar const&, tphysvar const&);
	tphysvar log10(tphysvar const&);
	tphysvar operator*(double const &, tphysvar const &);
	tphysvar sqrt(tphysvar const&);
	
	/**
	 * @brief The DataCube is the structure in which the model data needs to be loaded.
	 * 
	 * This DataCube object is used to represent a simulation in FoMo. The first part of a code using
	 * FoMo should consist of loading data into the DataCube.
	 */
	class DataCube 
	{
	protected:
		/** This is the dimension of the DataCube. */
		unsigned int dim;
		/** This is the number of variables stored in the DataCube. */
		unsigned int nvars;
		/** This is the number of grid points in the DataCube. */
		unsigned int ng;
		FoMo::tgrid grid;
		FoMo::tvars vars;
		void setgrid(tgrid ingrid);
		void setvar(const unsigned int, const tphysvar);
	public:
		DataCube(const int = 3);
		~DataCube();
		int readdim() const;
		int readngrid() const;
		int readnvars() const;
		tgrid readgrid() const;
		tphysvar readvar(const unsigned int) const;
		void setdim(const int indim);
		void setnvars(const int innvars);
		void setngrid(const int inngrid);
		void setdata(tgrid& ingrid, tvars& indata);
		void push_back(std::vector<double> coordinate, std::vector<double> variables);
	};
	
	/**
	 * @brief GoftCube contains the processed data, ready for rendering.
	 * 
	 * GoftCube is a class derived from the DataCube class, with additional information on wavelength
	 * (and corresponding data files) that were used for constructing the G(T) from a FoMo::DataCube
	 */
	class GoftCube: public DataCube
	{
	protected:
		std::string chiantifile;
		std::string abundfile;
//		std::string ion;
		double lambda0;
//		void setion(const std::string ion);
	public:
		GoftCube(const int = 3);
		GoftCube(DataCube datacube);
		void setchiantifile(const std::string inchianti);
		void setabundfile(const std::string inabund);
		std::string readchiantifile();
		std::string readabundfile();
		double readlambda0();
		void setlambda0(const double lambda0);
//		std::string readion();
		// keep the legacy ability to write Delaunay_triangulation
//		void writegoftcube(const std::string, const Delaunay_triangulation_3 *);
//		void readgoftcube(const std::string, Delaunay_triangulation_3*);
		void writegoftcube(const std::string);
		void readgoftcube(const std::string);
	};
	
	/**
	 * This enum allows the selection of the type of forward modelling that needs to be done. 
	 * Careful though, the value of this in the RenderCube is set by FoMoObject.render(), based on the resolution
	 * in the \f$\lambda\f$ direction (set with FoMoObject.setresolution()).
	*/
	enum FoMoObservationType
	{
		ObservationTypeNotDefined, /*!< This value should not be used: for coding purposes only!*/
		Spectroscopic, /*!< In this case, spectroscopic information is generated.*/
		Imaging /*!< In this case, only imaging information is obtained (e.g. using AIA filters).*/
	};
	
	/**
	 * @brief This class will contain resulting cubes from the rendering. 
	 * 
	 * The RenderCube class is a specialisation of the GoftCube class. It additionally contains information
	 * on the viewing angle, the resolution used and to-be-used for the rendering. Also, the rendermethod
	 * and observationtype is part of the RenderCube, and can be set and read with members of the class.
	 */
	class RenderCube: public GoftCube
	{
	protected:
		double l;
		double b;
		int x_pixel;
		int y_pixel;
		int z_pixel;
		int lambda_pixel;
		double lambda_width;
		std::string rendermethod;
		FoMoObservationType observationtype;
	public:
		RenderCube(GoftCube goftcube);
		void setresolution(const int & x_pixel, const int & y_pixel, const int & z_pixel, const int & lambda_pixel, const double & lambda_width);
		void readresolution(int & x_pixel, int & y_pixel, int & z_pixel, int & lambda_pixel, double & lambda_width);
		void setangles(const double l, const double b);
		void readangles(double & l, double & b);
		void setrendermethod(const std::string inrendermethod);
		std::string readrendermethod();
		void setobservationtype(FoMoObservationType);
		FoMoObservationType readobservationtype();
	};
	
	/**
	 * @brief FoMoObject is the main class of the FoMo library.
	 * 
	 * The class stores the output from the simulation results (internally as a FoMo::DataCube, 
	 * accessible via member FoMoObject.datacube), and converts it with the 
	 * FoMoObject::render()
	 * function to a FoMo::RenderCube, which can be read by the user as member FoMoObject.rendering.
	 */
	class FoMoObject
	{
	protected:
		std::string outfile;
		/**
		 @brief The datacube member if of type DataCube, storing the simulation data.
		  * 
		  * The FoMoObject.datacube member contains the data that is about to be rendered. Usually it 
		  * comes from a numerical simulation or a (semi-) analytical model. The FoMoObject.goftcube is 
		  * constructed from this data.
		*/
		FoMo::DataCube datacube;
		/**
		 @brief The goftcube member is of type GoftCube, with the G(T) values computed from datacube.
		  * 
		  * It contains the processed data from the DataCube, using the emission file and abundance 
		  * set with FoMoObject.setabundfile and FoMoObject.setchiantifile.
		*/
		FoMo::GoftCube goftcube;
		/**
		 @brief The rendering member is the final result of the FoMo simulation.
		  * 
		  * In the FoMoObject.rendering, the results of FoMoObject.render are stored. 
		*/
		FoMo::RenderCube rendering;
	public:
		FoMoObject(const int =3);
		~FoMoObject();
		void render(const double = 0, const double = 0); // l and b are arguments
		void render(const std::vector<double> lvec, const std::vector<double> bvec);
		void setrenderingdata(tgrid ingrid, tvars invars);
		FoMo::RenderCube readrendering();
		FoMo::GoftCube readgoftcube();
		void setrendermethod(const std::string inrendermethod);
		std::string readrendermethod();
		void setchiantifile(const std::string inchianti);
		void setabundfile(const std::string inabund);
		std::string readchiantifile();
		std::string readabundfile();
		void setoutfile(const std::string outfile);
		void push_back_datapoint(std::vector<double> coordinate, std::vector<double> variables);
		void setdata(tgrid& ingrid, tvars& indata);
		void setobservationtype(FoMoObservationType);
		FoMoObservationType readobservationtype();
		void setresolution(const int & x_pixel, const int & y_pixel, const int & z_pixel, const int & lambda_pixel, const double & lambda_width);
		void readresolution(int & x_pixel, int & y_pixel, int & z_pixel, int & lambda_pixel, double & lambda_width);
	};
}

#endif
