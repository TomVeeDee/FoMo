#include "../config.h"
#include "FoMo.h"
#include <cassert>

FoMo::DataCube::DataCube(const int indim)
{
	dim=indim;
	ng=0;
	nvars=0;
}

FoMo::DataCube::~DataCube()
{};

int FoMo::DataCube::readdim() const
{
	return dim;
};

int FoMo::DataCube::readngrid() const
{
	return ng;
};

int FoMo::DataCube::readnvars() const
{
	return nvars;
};

void FoMo::DataCube::setdim(const int indim)
{
	assert(dim > 0);
	dim=indim;
}

void FoMo::DataCube::setgrid(FoMo::tgrid ingrid)
{
	grid=ingrid;
};

FoMo::tgrid FoMo::DataCube::readgrid() const
{
	return grid;
};

void FoMo::DataCube::setvar(const unsigned int nvar, const FoMo::tphysvar var)
{
	assert(nvar < nvars);
	vars.reserve(nvars);
	vars[nvar]=var;
}

FoMo::tphysvar FoMo::DataCube::readvar(const unsigned int nvar) const
{
	assert(nvar < nvars);
	FoMo::tphysvar var=vars[nvar];
	return var;
}

void FoMo::DataCube::push_back(std::vector<double> coordinate, std::vector<double> variables)
{
	if (this->readngrid() == 0)
	{
		tgrid startgrid;
		for (unsigned int i=0; i<coordinate.size(); i++)
		{
			tcoord tempcoord{coordinate[i]};
			startgrid.push_back(tempcoord);
		}
		tvars startvars;
		for (unsigned int i=0; i<variables.size(); i++)
		{
			tphysvar tempvar{variables[i]};
			startvars.push_back(tempvar);
		}
		this->setdata(startgrid,startvars);
	}
	else
	{
		assert(coordinate.size() == dim);
		assert(variables.size() == nvars);

		// tgrid grid=this->readgrid();
		for (unsigned int i=0; i<dim; i++)
		{
			grid[i].push_back(coordinate[i]);
		}
	// this->setgrid(grid);
	
		for (unsigned int i=0; i<nvars; i++)
		{
			vars[i].push_back(variables[i]);
		}
	
		ng++;
	}
}

void FoMo::DataCube::setdata(tgrid ingrid, tvars indata)
{
	/* This routine sets the data of the DataCube from a simulation 
	 * (which could have been read in, for example)
	 * It loads the grid (tgrid ingrid) and variables (tvars indata).
	 * It checks if the tgrid and tvars all have columns of equal length.
	 * It sets the number of grid points (this->readngrid()) to the length of the column. 
	 * It sets the dimension of the grid (this->readdim()) to the number of columns in ingrid.
	 * It sets the number of physical variables (this->readnvars())
	 * to the number of columns in indata.
	*/
	
	//perform checks
	assert(ingrid.size() != 0); // The dimension is 0. This shouldn't work.
	dim=ingrid.size();
	
	unsigned int tempng=ingrid[0].size();
	for (unsigned int i=1; i<dim; i++)
	{
		// all the coordinate vectors should be equally long
		assert(ingrid[i].size() == tempng);
	}
	
	ng=tempng;
	
	this->setgrid(ingrid);
	

	// first check if there are variables (i.e. if indata is not empty)
	if (indata.size() !=0)
	{
		// set nvars
		nvars=indata.size();
		for (unsigned int i=0; i<nvars; i++)
		{
			assert(indata[i].size() == ng);
			this->setvar(i,indata[i]);
		}
	}
}