#include <vector>
#include <string>

namespace FoMo
{
	typedef std::vector<double> tcoord;
	typedef tcoord * tgrid;
	typedef std::vector<double> tphysvar;
	typedef tphysvar * tvars;
	
	class SimulationPoint
	{
	protected:
		int dim;
		std::vector<double> position;
		// ordered vector of physical properties (rho, T, vx, vy, vz)
		std::vector<double> physicalvariables;
	public:
		SimulationPoint(const int = 3);
		~SimulationPoint();
		void setdimension(const int);
		int dimension();
		// maybe better to use std::array<T,dim> coords as argument, instead of overloading. This could also be used as internal representation.
		void setposition(const double x);
		void setposition(const double x, const double y);
		void setposition(const double x, const double y, const double z);
		void readposition(double x, double y, double z) const;
		void setvars(const double rho, const double T, const double vx, const double vy, const double vz);
		void readvars(double rho, double T, double vx, double vy, double vz) const;
	};
	
	class cube 
	{
	protected:
		int dim;
		int ng;
		int nvars;
		FoMo::tgrid grid;
		FoMo::tvars vars;
	public:
		cube(const int invars, const int ingrid, const int = 3);
		~cube();
		int readdim() const;
		int readngrid() const;
		int readnvars() const;
		void setgrid(tgrid ingrid);
		tgrid readgrid() const;
		void setvar(const int, const tphysvar);
		tphysvar readvar(const int) const;
	};

	class FoMoObject
	{
	protected:
		FoMo::cube simulation;
		FoMo::cube rendering;
		std::string rendermethod;
	public:
		FoMoObject();
		~FoMoObject();
		void setvar();
		void setgrid();
		void push_back(FoMo::SimulationPoint simulationpoint);
		void setrendermethod(const std::string rendermethod);
		std::string readrendermethod();
		void render(const double l, const double b);
		FoMo::cube readrendering();
	};
	
}