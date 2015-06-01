#include "FoMo.h"
#include <cassert>

/*
This files defines the member functions of the FoMo::SimulationPoint class. Theoretically, it defines a number of variables at a location (of which the dimension is stored in int dim).
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
		void setposition(const double x, const double y, const double z);
		void readposition(double x, double y, double z) const;
		void setvars(const double rho, const double T, const double vx, const double vy, const double vz);
		void readvars(double rho, double T, double vx, double vy, double vz) const;
	};*/

FoMo::SimulationPoint::SimulationPoint(const int dimension)
{
	dim=dimension;
}

FoMo::SimulationPoint::~SimulationPoint()
{
}

void FoMo::SimulationPoint::setdimension(const int dimension)
{
	dim=dimension;
}

int FoMo::SimulationPoint::dimension()
{
	return dim;
}

void FoMo::SimulationPoint::setposition(const double x)
{
	assert(dim==1);
	position.push_back(x);
}

void FoMo::SimulationPoint::setposition(const double x, const double y)
{
	assert(dim==2);
	position.push_back(x);
	position.push_back(y);
}

void FoMo::SimulationPoint::setposition(const double x, const double y, const double z)
{
	assert(dim==3);
	position.push_back(x);
	position.push_back(y);
	position.push_back(z);
}

void FoMo::SimulationPoint::readposition(double x, double y, double z) const
{
	

}