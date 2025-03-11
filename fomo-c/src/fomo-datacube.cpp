#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <cassert>

/**
 * @brief The default constructor for a DataCube.
 * 
 * As part of the constructor, the dimension is set to indim, the grid is resized to contain indim 
 * coordinate vectors, the number of grid points and variables is set to 0. 
 * @param indim The dimension of the DataCube. It defaults to 3.
 */
FoMo::DataCube::DataCube(const int indim)
{
	dim=indim;
	grid.resize(indim);
	ng=0;
	nvars=0;
}

/**
 * @brief The destructor for DataCube.
 */
FoMo::DataCube::~DataCube()
{};

/**
 * @brief This returns the current dimension of the DataCube. 
 * @return This returns the dimension of the DataCube as an integer.
 */
int FoMo::DataCube::readdim() const
{
	return dim;
};

/**
 * @brief This returns the number of data points in the DataCube.
 * @return This returns the number of data points in the DataCube as an integer.
 */
int FoMo::DataCube::readngrid() const
{
	return ng;
};

/**
 * @brief This returns the number of variables in the DataCube.
 * @return The number of variables present in the DataCube is returned as an integer.
 */
int FoMo::DataCube::readnvars() const
{
	return nvars;
};

/**
 * @brief This function sets the dimension of the DataCube.
 * 
 * The dimension is set to the positive integer (the argument). The grid is resized to have indim 
 * columns. This means that if the DataCube contained already a smaller grid, extra coordinates are 
 * created. If 
 * the grid was larger than indim, coordinates will be removed.
 * @param indim The parameter must be a positive integer and is the dimension of the DataCube.
 */
void FoMo::DataCube::setdim(const int indim)
{
	assert(dim > 0);
	dim=indim;
	grid.resize(indim);
}

/**
 * @brief This sets the number of variables in the DataCube.
 * 
 * If the DataCube already contained variables, then the data is lost (if nvars is decreased from its
 * previous state), or extra variables are created.
 * @param innvars A positive integer that the number of variables will be set to.
 */
void FoMo::DataCube::setnvars(const int innvars)
{
	assert(nvars >= 0);
	nvars=innvars;
	vars.resize(innvars);
}

/**
 * @brief This sets the number of grid points.
 * 
 * If the grid present in the DataCube was longer than inngrid, data will be lost. If the grid was
 * shorter, extra spaces are created. 
 * @param inngrid The number of grid points. It should be a positive integer.
 */
void FoMo::DataCube::setngrid(const int inngrid)
{
	assert(inngrid >= 0);
	ng=inngrid;
	for (unsigned int i=0; i<dim; i++)
	{
		grid.at(i).resize(inngrid);
	}
}

/**
 * @brief This routine sets the grid of the DataCube.
 * 
 * Additionally, the number of grid points is updated to the length of the vectors in ingrid. The 
 * dimension is set to the length of ingrid. \n
 * This routine is protected. It should not be used by end
 * users. The data should be imported by setdata() instead. 
 * @param ingrid The parameter is a vector of vectors of doubles (tgrid). It should have at least one 
 * element.
 * The elements of ingrid should have equal length.
 */
void FoMo::DataCube::setgrid(FoMo::tgrid ingrid)
{
	//perform checks
	assert(ingrid.size() != 0); // The dimension is 0. This shouldn't work.
	dim=ingrid.size();
	
	unsigned int tempng=ingrid[0].size();
	for (unsigned int i=1; i<dim; i++)
	{
		// all the coordinate vectors should be equally long
		assert(ingrid.at(i).size() == tempng);
	}
	
	grid=ingrid;
	ng=tempng;
};

/**
 * @brief This function allows reading of the grid.
 * @return The grid is returned as a tgrid.
 */
FoMo::tgrid FoMo::DataCube::readgrid() const
{
	return grid;
};

/**
 * @brief This sets the specific variable to the argument vector.
 * 
 * Since not a lot of checks are implemented, this routine is protected. It should not be used by end
 * users. The data should be imported by setdata() instead. 
 * @param nvar This is the variable that should be set. It's value can be 0, and goes to nvars-1. 
 * @param var This is the content of the variable to be set. Even though it is not checked, the length should
 * be ng (number of grid points). The responsibility is left to the programmer to allow greater flexibility
 * in importing the numerical data. If the check would be enforced, first the grid needs to be set with
 * setgrid().
 */
void FoMo::DataCube::setvar(const unsigned int nvar, const FoMo::tphysvar var)
{
	assert(nvar < vars.size());
	vars.at(nvar)=var;
}

/**
 * @brief This reads a specific variable in the DataCube.
 * @param nvar The variable that should be read. It should be between 0 and nvars-1.
 * @return The physical variable is returned as a vector.
 */
FoMo::tphysvar FoMo::DataCube::readvar(const unsigned int nvar) const
{
	assert(nvar < nvars);
	FoMo::tphysvar var=vars.at(nvar);
	return var;
}

/**
 * @brief This adds a data point to the DataCube.
 * 
 * If it is the first data point added to the DataCube, the ngrid is set to 1, the dimension is set to 
 * the length of coordinate, and the nvars is set to the length of variables. \n
 * If there is already data in the DataCube, then it is checked if the length of coordinates matches the 
 * dimension, and if the length of the variables matches nvars. The point is added, and ngrid is increased by 1.
 * @param coordinate This contains the coordinates of the data point. It should have length dimension,
 * unless there is no data yet in DataCube.
 * @param variables This contains the physical variables at the data point. It should have length
 * nvars, unless there is no data yet in the DataCube.
 */
void FoMo::DataCube::push_back(std::vector<double> coordinate, std::vector<double> variables)
{
	if (this->readngrid() == 0)
	{
		tgrid startgrid;
		for (unsigned int i=0; i<coordinate.size(); i++)
		{
			tcoord tempcoord{coordinate.at(i)};
			startgrid.push_back(tempcoord);
		}
		tvars startvars;
		for (unsigned int i=0; i<variables.size(); i++)
		{
			tphysvar tempvar{variables.at(i)};
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
			grid.at(i).push_back(coordinate.at(i));
		}
	// this->setgrid(grid);
	
		for (unsigned int i=0; i<nvars; i++)
		{
			vars.at(i).push_back(variables.at(i));
		}
	
		ng++;
	}
}

/**
 * @brief This routine sets the data of the DataCube from a simulation.
 * 
 * It loads the grid (tgrid ingrid) and variables (tvars indata).
 * It checks if the tgrid and tvars all have columns of equal length.
 * It sets the number of grid points (this->readngrid()) to the length of the column. 
 * It sets the dimension of the grid (this->readdim()) to the number of columns in ingrid.
 * It sets the number of physical variables (this->readnvars())
 * to the number of columns in indata.
 * @param ingrid The grid to be imported in DataCube.
 * @param indata The data to be imported in DataCube.
 */
void FoMo::DataCube::setdata(tgrid& ingrid, tvars& indata)
{
	this->setgrid(ingrid);
	
	// first check if there are variables (i.e. if indata is not empty)
	if (indata.size() !=0)
	{
		// set nvars
		nvars=indata.size();
		vars.resize(nvars,tphysvar(ng,0));
		for (unsigned int i=0; i<nvars; i++)
		{
			assert(indata[i].size() == ng);
			this->setvar(i,indata[i]);
		}
	}
}