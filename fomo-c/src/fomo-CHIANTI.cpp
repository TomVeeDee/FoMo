#include "../config.h"
#include "FoMo.h"
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
	typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
	typedef CGAL::Interpolation_traits_2<K>               Traits;
	typedef K::FT                                         Coord_type;
	typedef K::Point_2                                    Point;
        std::map<Point, Coord_type, K::Less_xy_2> function_values;
        typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2 > >  Value_access;

	FoMo::tgrid grid=gofttab.readgrid();
	FoMo::tphysvar tempgrid=grid[0];
	FoMo::tphysvar rhogrid=grid[1];
	FoMo::tphysvar goftvec=gofttab.readvar(0);

	// put the temperature and density grid from the gofttab into the CGAL data structures
	std::vector<Point> delaunaygrid;
	Point temporarygridpoint;

	FoMo::tphysvar::const_iterator tempit=tempgrid.begin();
	FoMo::tphysvar::const_iterator rhoit=rhogrid.begin();
	FoMo::tphysvar::const_iterator goftit=goftvec.begin();

	for (; tempit != tempgrid.end(); ++tempit)
	{
		temporarygridpoint=Point(*tempit,*rhoit);
		delaunaygrid.push_back(temporarygridpoint);
		function_values.insert(std::make_pair(temporarygridpoint,Coord_type(*goftit)));
		++rhoit;
		++goftit;
	}
	Value_access tempmap=Value_access(function_values);

	// compute the Delaunay triangulation
        Delaunay_triangulation DT;
	DT.insert(delaunaygrid.begin(),delaunaygrid.end());

	FoMo::tphysvar g;
	int ng=logT.size();
	// the reservation of the size is necessary, otherwise the vector reallocates in the openmp threads with segfaults as consequence
	g.resize(ng);
	std::cout << "Doing G(T) interpolation: " << std::flush;

	// If MPI is available, we should divide the work between nodes, we can just use a geometric division by the commsize
	int commrank,commsize;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
        MPI_Comm_size(MPI_COMM_WORLD,&commsize);
#else
	commrank = 0;
	commsize = 1;
#endif
	int mpimin, mpimax;
	// MPE_Decomp1d decomposes a vector into more or less equal parts
	// The only restriction is that it takes the vector from 1:ng, rather than 0:ng-1
	// Therefore, the mpimin and mpimax are decreases afterwards.
	MPE_Decomp1d(ng,commsize, commrank, &mpimin, &mpimax);
	mpimin--;
	mpimax--;

	boost::progress_display show_progress((mpimax-mpimin+1)/10);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i=mpimin; i<mpimax; i++)
	{
		Point p(logT[i],logrho[i]);
// a possible linear interpolation
// see http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Interpolation/Chapter_main.html#Subsection_70.3.3
// make sure to use small interpolation file. If not, this takes hours!
  
  		std::vector< std::pair< Point, Coord_type > > coords;
		Coord_type norm = CGAL::natural_neighbor_coordinates_2(DT, p,std::back_inserter(coords)).second;

		Coord_type res =  CGAL::linear_interpolation(coords.begin(), coords.end(),
                                               norm,
                                               Value_access(function_values));	
		g[i]=res;

// let's not do nearest neighbour
// it is much faster than the linear interpolation
// but the spacing in temperature is too broad too handle this

		/*Delaunay_triangulation::Vertex_handle v=DT.nearest_vertex(p);
		Point nearest=v->point();
		std::pair<Coord_type,bool> funcval=tempmap(nearest);
		g[i]=funcval.first;*/
	
		// This introduces a race condition between threads, but since it's only a counter, it doesn't really matter.
		/*#ifdef _OPENMP
		#pragma omp atomic
		#endif*/
		if ((i-mpimin)%10 == 0) ++show_progress;
	}
	std::cout << " Done!" << std::endl << std::flush;

	return g;
}

double FoMo::readgoftfromchianti(const std::string chiantifile)
{
	std::ifstream in(chiantifile);
	if (!in) {
		std::cerr << "Error: no CHIANTI G(T) file exists\n";
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
		std::cerr << "Error: no CHIANTI G(T) file exists\n";
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
	FoMo::tphysvar w=FoMo::operator*(2*sqrt(2*log(2))*sqrt(2*boltzmannconstant/massproton/atweight)*lambda0/speedoflight,FoMo::sqrt(T));
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

