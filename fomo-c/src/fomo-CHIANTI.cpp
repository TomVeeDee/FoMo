#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>

// CGAL stuff
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Physical constants
const double alphaconst=1.0645; // alpha=\sqrt{2\pi}/2/\sqrt{2ln2}
const double boltzmannconstant=GSL_CONST_MKSA_BOLTZMANN; // Boltzmann constant
const double massproton=GSL_CONST_MKSA_MASS_PROTON; // hydrogen mass
const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light

const std::string ______chiantitables_sun_coronal_abund_string;

int MPE_Decomp1d(int n, int size, int rank, int *s, int *e )
/*
~  This file contains a routine for producing a decomposition of a 1-d array
~  when given a number of processors.  It may be used in "direct" product
~  decomposition.  The values returned assume a "global" domain in [1:n]
~ */
/*@
~  MPE_Decomp1d - Compute a balanced decomposition of a 1-D array

~  Input Parameters:
+ n  - Length of the array
. size - Number of processors in decomposition
- rank - Rank of this processor in the decomposition (0 <= rank < size)

~  Output Parameters:
. s,e - Array indices are s:e, with the original array considered as 1:n.
@*/
{
	int nlocal, deficit;

	nlocal      = n / size;
	*s  = rank * nlocal + 1;
	deficit     = n % size;
	*s  = *s + ((rank < deficit) ? rank : deficit);
	if (rank < deficit) nlocal++;
	*e      = *s + nlocal - 1;
	if (*e > n || rank == size-1) *e = n;
#ifdef HAVEMPI
	return MPI_SUCCESS;
#else
	return EXIT_SUCCESS;
#endif
}

double abundfromchianti(std::istream & in, const std::string & ion)
{
	// read the abundfile and return the value that goes with the ion string
	if (!in) {
		std::cerr << "Error: no CHIANTI abundance file exists\n";
		exit(EXIT_FAILURE);
	}

	int nelement;
	double value,hvalue,ionvalue;
	hvalue=0;
	ionvalue=0;
	std::string element;
	size_t endion=ion.find_first_of("_");
	std::string ionname=ion.substr(0,endion);
	in >> nelement;
	while (nelement >= 0)
	{
		in >> value;
		in >> element;
		if ((element.size() == 1) && (strncasecmp(element.c_str(),"h",1) == 0))
		{
			hvalue=value;
		}
		if ((element.size() == ionname.size()) && (strncasecmp(element.c_str(),ionname.c_str(),ionname.size()) == 0))
		{
			ionvalue=value;
		}
		in >> nelement;
	}
	if ( (hvalue == 0.) || (ionvalue == 0.) ) 
	{
		std::cerr << "Error: invalid CHIANTI abundance file\n";
		exit(EXIT_FAILURE);
	}

	value = pow(10.,(ionvalue - hvalue));

	return value;
}

FoMo::tphysvar goft(const FoMo::tphysvar logT, const FoMo::tphysvar logrho, const FoMo::DataCube gofttab)
{
	// uses the log(T), because the G(T) is also stored using those values.

	FoMo::tgrid grid=gofttab.readgrid();
	// this is sorted by a slowly varying density and quickly varying temperature
	FoMo::tphysvar tempgrid=grid[0];
	FoMo::tphysvar rhogrid=grid[1];
	// we can find the elementary vector length by using std::count
	unsigned int nrho=std::count(tempgrid.begin(),tempgrid.end(),tempgrid.at(0));
	unsigned int nt=std::count(rhogrid.begin(),rhogrid.end(),rhogrid.at(0));
	assert(tempgrid.size()==nt*nrho);
	FoMo::tphysvar goftvec=gofttab.readvar(0);

	FoMo::tphysvar g;
	int ng=logT.size();
	// the reservation of the size is necessary, otherwise the vector reallocates in the openmp threads with segfaults as consequence
	g.resize(ng);
	std::cout << "Doing G(T) interpolation: " << std::flush;

	unsigned int floortemp, ceiltemp, rhoindex,goftindex;
	double res, x1, x2, y1, y2;
	boost::progress_display show_progress((ng+1)/10);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i=0; i<ng; i++)
	{
		// find index of logT[i] in tempgrid[0:nt-1]
		// find index of logrho[i] in rhogrid[0:nrho-1]
		// tempindex=round((logT.at(i)-tempgrid.at(0))/(tempgrid.at(1)-tempgrid.at(0)));
		floortemp=floor((logT.at(i)-tempgrid.at(0))/(tempgrid.at(1)-tempgrid.at(0)));
		ceiltemp=floortemp+1;
		rhoindex=round((logrho.at(i)-rhogrid.at(0))/(rhogrid.at(nt)-rhogrid.at(0)));
		
		// if these indices are outside the domain, then the G(T) is 0
		if (floortemp<0 || rhoindex<0 || ceiltemp>nt-1 || rhoindex>nrho-1)
		{
			res=0;
		}
		else
		{
			goftindex=rhoindex*nt+ceiltemp;
		// we do a linear interpolation
		// Patrick does a linear interpolation in temperature, and a nearest neighbour in density
			x1=tempgrid.at(ceiltemp);
			y1=goftvec.at(goftindex);
			x2=tempgrid.at(floortemp);
			y2=goftvec.at(goftindex-1);
			res=(logT.at(i)-x1)*(y2-y1)/(x2-x1)+y1;
		}
		
		g.at(i)=res;
		
		if (i % 10 == 0) ++show_progress;
	}
	std::cout << " Done!" << std::endl << std::flush;

	return g;
}

double FoMo::readgoftfromchianti(const std::string chiantifile)
{
	std::ifstream in(chiantifile);
	if (!in) {
		std::cerr << "Error: no CHIANTI G(T) file exists: " << chiantifile << "\n";
		exit(EXIT_FAILURE);
	}

	double lambda0;
	std::string ion;

	// read in first two lines (ion name, rest wavelength)
	in >> ion;
	in >> lambda0;

	in.close();

	return lambda0;
}

FoMo::DataCube FoMo::readgoftfromchianti(const std::string chiantifile, std::string & ion, double & lambda0, double & atweight)
{
	std::ifstream in(chiantifile);
	if (!in) {
		std::cerr << "Error: no CHIANTI G(T) file exists: " << chiantifile << "\n";
		exit(EXIT_FAILURE);
	}

	double templogrho;
	double nrhotemp;
	int nrho, nt;
	double field;
	FoMo::tphysvar temprho, tempt, tempgoft;

	std::cout << "Reading G(T) from " << chiantifile << "... " << std::flush;
	// read in header (ion name, rest wavelength, atomic weight, number of grid points in density direction, number of grid points in temperature direction)
	in >> ion;
	in >> lambda0;
	in >> atweight;
	in >> nrhotemp;
	nrho=round(nrhotemp);
	in >> nt;

	// read in the temperature grid points
	FoMo::tphysvar tvec;

	for (int i=0; i<nt; i++)
	{
		in >> field;
		tvec.push_back(field);
	}

	while (in >> templogrho)
	{
		for (int i=0; i<nt; i++)
		{
			in >> field;
			// we also need to push back the density and temperature
			tempt.push_back(tvec[i]);
			temprho.push_back(templogrho);
			tempgoft.push_back(field);
		}
	}

	unsigned int ntnrho=nt*nrho;
	if (tempt.size() != ntnrho) 
	{
		std::cerr << "Error: CHIANTI table not read in properly!"<< std::endl;
		exit(EXIT_FAILURE);
	}
	FoMo::DataCube gofttab(2); //only 2 dimensions
	
	FoMo::tgrid tempgrid(2);
	// Warning this assumes that there is a straightforward mapping between tcoord and tphysvar!
	tempgrid[0]=tempt;
	tempgrid[1]=temprho;
	
	FoMo::tvars tempvars{tempgoft};
	
	gofttab.setdata(tempgrid,tempvars);

	in.close();

	std::cout << "Done!" << std::endl << std::flush;

	return gofttab;
}

FoMo::tphysvar linefwhm(const FoMo::tphysvar T, const double lambda0, const double atweight)
{
// uses the absolute temperature! i.e. not log(T), but T
// for spectroscopic study return, line width
	FoMo::tphysvar w=FoMo::operator*(2*sqrt(2*log(2))*sqrt(boltzmannconstant/massproton/atweight)*lambda0/speedoflight,FoMo::sqrt(T));
 	return w;
}

FoMo::GoftCube FoMo::emissionfromdatacube(FoMo::DataCube datacube, std::string chiantifile, std::string abundfile, const FoMoObservationType observationtype)
{
// construct the G(T) using CHIANTI tables
	int nvars=datacube.readnvars(); // G(T), width, vx, vy, vz
	// take the last 3 from datacube
	//std::cout << nvars << datacube.vars[0][0] << std::flush;
	
	FoMo::GoftCube emission;
	
	FoMo::tphysvar logrho = FoMo::log10(datacube.readvar(0));
	FoMo::tphysvar T = datacube.readvar(1);
	FoMo::tphysvar logT = FoMo::log10(T);
	
	std::string ion="";// better to be aia for AIA tables this will be checked for correct table [DY]
	double lambda0=.0;
	double atweight=1;
	FoMo::DataCube gofttab=readgoftfromchianti(chiantifile,ion,lambda0,atweight);


// fit the G(T) function to the data points
	FoMo::tphysvar fittedgoft=goft(logT,logrho,gofttab);
// calculate the line width at each data point
	FoMo::tphysvar fittedwidth=linefwhm(T,lambda0,atweight);

// calculate the maximum of the emission profile at each grid point
// The emission in the CHIANTItables is calculated with the sun_coronal.abund file
// There are 2 cases:
// 	- we are doing spectroscopic calculations: the fittedemission should normalised with the alphaconst, and if a different abundance is used, we should renormalise to the new abundance
// 	- we are doing intensity calculations (for e.g. AIA): the tables are already in the correct units, and the fittedemission needs to only be multiplied with 1. 
	double normaliseconst=1.;
	double abundratio=1.;
	if (observationtype == Spectroscopic) // for spectroscopic study
	{
		std::cout << "This code is configured to do spectroscopic modelling." << std::endl;        
		// FoMo::tphysvar fittedwidth=linefwhm(T,lambda0,atweight);
		normaliseconst=1./alphaconst;
		if (abundfile != std::string("/empty"))
		{
			std::ifstream abundfilestream (abundfile);
			std::cout << "Reading abundances from " << abundfile << "... " << std::flush;
			double abundnew=abundfromchianti(abundfilestream, ion);
			std::cout << "Done!" << std::endl << std::flush;
			std::istringstream standardabundfile (______chiantitables_sun_coronal_abund_string);
			double abundold=abundfromchianti(standardabundfile,ion);
			abundratio=abundnew/abundold;
		}
	}
        if (observationtype == Imaging) // for imaging study 
        {
        
        if (atoi(ion.c_str())!=int(lambda0+0.5)) 
          {// check if it is AIA GOFT table
           std::cerr << "GOFT table is not correct!" << std::endl; 
           exit(EXIT_FAILURE);
         }
 
          fittedwidth=FoMo::operator /(T,T);// =1.0
        }

// Spectral modelling: ne^2 [cm^-3] * G(ne,Te) [erg cm^3 s^1] = emis [erg cm^-3 s^-1] 
// for integration *[cm] [D.Y 12 Nov 2014]
// AIA modelling: Temperature response function K(ne,Te)[cm^5 DN S^-1]*ne[cm^-3]=emis [DN cm^-1 s^-1]; 
// for integration *[cm] [D.Y 12 Nov 20
//        std::cout << "normaliseconst" <<normaliseconst << std::endl;
//        std::cout << "abundratio"<<abundratio<<std::endl;
	FoMo::tphysvar fittedemission=FoMo::operator *(normaliseconst*abundratio,FoMo::operator *(FoMo::pow(10.,FoMo::operator *(2.,logrho)),fittedgoft));
	fittedemission=FoMo::operator /(fittedemission,fittedwidth);

// load the emission and width into the emission-cube variable
	FoMo::tgrid grid=datacube.readgrid();
	FoMo::tvars exporteddata;
	exporteddata.push_back(fittedemission);
	exporteddata.push_back(fittedwidth);
	
	// copy velocity vectors
	for (int i=2; i<nvars; i++)
	{
			FoMo::tphysvar velcomp=datacube.readvar(i);
			exporteddata.push_back(velcomp);
	};
	emission.setdata(grid,exporteddata);
	emission.setchiantifile(chiantifile);
	emission.setabundfile(abundfile);
	emission.setlambda0(lambda0);

	return emission;
};

