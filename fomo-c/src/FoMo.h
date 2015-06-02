#include <vector>
#include <string>
#include <cassert>

namespace FoMo
{
	typedef std::vector<double> tcoord;
	typedef std::vector<tcoord> tgrid;
	typedef std::vector<double> tphysvar;
	typedef std::vector<tphysvar> tvars;
	
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
		void setdata(tgrid ingrid, tvars indata);
		void push_back(std::vector<double> coordinate, std::vector<double> variables);
	};
	
	// included for backwards compatibility
	typedef DataCube cube;

	class FoMoObject
	{
	protected:
		std::string rendermethod;
		FoMo::DataCube rendering;
	public:
		FoMoObject();
		~FoMoObject();
		FoMo::DataCube datacube;
		void setrendermethod(const std::string inrendermethod);
		std::string readrendermethod();
		void render(const double = 0, const double = 0); // l and b are arguments
		void render(const std::vector<double> lvec, const std::vector<double> bvec);
		FoMo::DataCube readrendering();
	};
	
	FoMo::DataCube RenderWithCGAL(FoMo::DataCube datacube, const std::vector<double> lvec, const std::vector<double> bvec);
}