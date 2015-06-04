#include <vector>
#include <string>

namespace FoMo
{
	typedef std::vector<double> tcoord;
	typedef std::vector<tcoord> tgrid;
	typedef std::vector<double> tphysvar;
	typedef std::vector<tphysvar> tvars;
	
	tphysvar pow(const double, tphysvar const&);
	tphysvar operator/(tphysvar const&, tphysvar const&);
	tphysvar operator*(tphysvar const&, tphysvar const&);
	tphysvar log10(tphysvar const&);
	tphysvar operator*(double const &, tphysvar const &);
	tphysvar sqrt(tphysvar const&);
	
	class DataCube 
	{
	protected:
		unsigned int dim;
		unsigned int nvars;
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
		void setdata(tgrid& ingrid, tvars& indata);
		void push_back(std::vector<double> coordinate, std::vector<double> variables);
	};
	// included for backwards compatibility
	typedef DataCube cube;
	
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
	};
	
	enum FoMoObservationType
	{
		ObservationTypeNotDefined,
		Spectroscopic,
		Imaging
	};
	
	class RenderCube: public GoftCube
	{
	protected:
		double l;
		double b;
		int x_pixel;
		int y_pixel;
		int z_pixel;
		int lambda_pixel;
		int lambda_width;
		std::string rendermethod;
		FoMoObservationType observationtype;
	public:
		RenderCube(GoftCube goftcube);
		void setresolution(const int & nx, const int & ny, const int & nz, const int & nlambda, const double & lambdawidth);
		void readresolution(int & nx, int & ny, int & nz, int & nlambda, double & lambdawidth);
		void setangles(const double l, const double b);
		void readangles(const double l, const double b);
		void setrendermethod(const std::string inrendermethod);
		std::string readrendermethod();
		void setobservationtype(FoMoObservationType);
		FoMoObservationType readobservationtype();
	};
	
	class FoMoObject
	{
	public:
		FoMo::DataCube datacube;
	protected:
		FoMo::GoftCube goftcube;
	public:
		FoMo::RenderCube rendering;
		FoMoObject(const int =3);
		~FoMoObject();
		void render(const double = 0, const double = 0); // l and b are arguments
		void render(const std::vector<double> lvec, const std::vector<double> bvec);
		void setrenderingdata(tgrid ingrid, tvars invars);
		FoMo::RenderCube readrendering();
		void setrendermethod(const std::string inrendermethod);
		std::string readrendermethod();
		void setchiantifile(const std::string inchianti);
		void setabundfile(const std::string inabund);
		std::string readchiantifile();
		std::string readabundfile();
	};
	
	double readgoftfromchianti(const std::string chiantifile);
	DataCube readgoftfromchianti(const std::string chiantifile, std::string & ion, double & lambda0, double & atweight);
	GoftCube emissionfromdatacube(DataCube, std::string, std::string, const FoMoObservationType);
	
	FoMo::RenderCube RenderWithCGAL(FoMo::DataCube datacube, std::string chiantifile, std::string abundfile, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec);
}