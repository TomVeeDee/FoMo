#include "header.h"
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_const_mksa.h>

// CGAL stuff
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>



// This should become part of the input options, or based on the choice of line position.
const char* chiantifile="chiantitables/goft_table_f2rt_171.dat";

// Physical constants
const double alphaconst=1.0645; // alpha=\sqrt{2\pi}/2/\sqrt{2ln2}
const double boltzmannconstant=GSL_CONST_MKSA_BOLTZMANN; // Boltzmann constant
const double massproton=GSL_CONST_MKSA_MASS_PROTON; // hydrogen mass
const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light


tphysvar goft(const tphysvar logT, const tphysvar logrho, const cube gofttab)
{
	// uses the log(T), because the G(T) is also stored using those values.
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
	typedef CGAL::Interpolation_traits_2<K>               Traits;
	typedef K::FT                                         Coord_type;
	typedef K::Point_2                                    Point;
        std::map<Point, Coord_type, K::Less_xy_2> function_values;
        typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2 > >  Value_access;

	tgrid grid=gofttab.readgrid();
	tphysvar tempgrid=grid[0];
	tphysvar rhogrid=grid[1];
	tphysvar goftvec=gofttab.readvar(0);

	// put the temperature and density grid from the gofttab into the CGAL data structures
	vector<Point> delaunaygrid;
	Point temporarygridpoint;

	tphysvar::const_iterator tempit=tempgrid.begin();
	tphysvar::const_iterator rhoit=rhogrid.begin();
	tphysvar::const_iterator goftit=goftvec.begin();

	for (; tempit != tempgrid.end(); ++tempit)
	{
		temporarygridpoint=Point(*tempit,*rhoit);
		delaunaygrid.push_back(temporarygridpoint);
		function_values.insert(std::make_pair(temporarygridpoint,Coord_type(*goftit)));
		++rhoit;
		++goftit;
	}

	// compute the Delaunay triangulation
        Delaunay_triangulation DT;
	DT.insert(delaunaygrid.begin(),delaunaygrid.end());

	tphysvar::const_iterator logtit=logT.begin();
	tphysvar::const_iterator logrhoit=logrho.begin();
	tphysvar g;
	int i=0;
	int ng=logT.size();
	cout << "Doing interpolation:";
	for (; logtit != logT.end(); ++logtit)
	{
		Point p(*logtit,*logrhoit);
// a possible linear interpolation, perhaps we should check if this eats a lot of time
// see http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Interpolation/Chapter_main.html#Subsection_70.3.3
/*  
  		std::vector< std::pair< Point, Coord_type > > coords;
		Coord_type norm = CGAL::natural_neighbor_coordinates_2(DT, p,std::back_inserter(coords)).second;

		Coord_type res =  CGAL::linear_interpolation(coords.begin(), coords.end(),
                                               norm,
                                               Value_access(function_values));	
		g.push_back(res); */

// let's do nearest neighbour

		Delaunay_triangulation::Vertex_handle v=DT.nearest_vertex(p);
		Value_access tempmap=Value_access(function_values);
		Point nearest=v->point();
		std::pair<Coord_type,bool> funcval=tempmap(nearest);
		g.push_back(funcval.first);
		
		progressbar(i,0,ng-1);
		++logrhoit;
		i++;
	}
	cout << "Done!" << endl;

	return g;
}

cube readgoftfromchianti(const char* chiantifile, string & ion, double & lambda0, double & atweight)
{
	ifstream in(chiantifile);
	if (!in) {
		cout << "Error: no CHIANTI G(T) file exists\n";
		exit(EXIT_FAILURE);
	}

	double templogrho;
	int nt;
	int nrho=0;
	double field;
	tphysvar temprho, tempt, tempgoft;
// this should be inserted after Patrick has updated the G(T) tables.
//	in >> ion;
//	in >> lambda0;
//	in >> atweight;
	ion="feIX";
	lambda0=171.073;
	atweight=56.;
	while (in >> templogrho)
	{
		nrho++;
		in >> nt;
		for (int i=0; i<nt; i++)
		{
			in >> field;
			tempt.push_back(field);
			temprho.push_back(templogrho);
		}
		for (int i=0; i<nt; i++)
		{
			in >> field;
			tempgoft.push_back(field);
		}
	}

	int ntnrho=nt*nrho;
	int nvars=1; // G(T) and width
	cube gofttab(nvars,ntnrho,2); //only 2 dimensions
	gofttab.settype(gofttable);
	
	tgrid tempgrid=new tcoord[2];
	// Warning this assumes that there is a straightforward mapping between tcoord and tphysvar!
	tempgrid[0]=tempt;
	tempgrid[1]=temprho;
	
	gofttab.setgrid(tempgrid);
	gofttab.setvar(0,tempgoft);

	return gofttab;
}

tphysvar linefwhm(const tphysvar T, const double lambda0, const double atweight)
{
	// uses the absolute temperature! i.e. not log(T), but T
	tphysvar w=2*sqrt(2*log(2))*sqrt(boltzmannconstant/massproton/atweight)*lambda0/speedoflight*sqrt(T);
	return w;
}

cube emissionfromdatacube(cube datacube)
{
// construct the G(T) using CHIANTI tables
	int ng=datacube.readngrid();
	int nvars=5; // G(T), width, vx, vy, vz
	// take the last 3 from datacube
	
	cube emission(nvars, ng);
	tgrid grid=datacube.readgrid();
	emission.setgrid(grid);
	emission.settype(emisscube);

	// copy velocity vectors
	for (int i=2; i<nvars; i++)
	{
		tphysvar velcomp=datacube.readvar(i);
		emission.setvar(i,velcomp);
	};
	
	tphysvar logrho = log10(datacube.readvar(0));
	tphysvar T = datacube.readvar(1);
	tphysvar logT = log10(T);
	writearray(logT,"empty");
	
	string ion="";
	double lambda0=.0;
	double atweight=1;
	cube gofttab=readgoftfromchianti(chiantifile,ion,lambda0,atweight);

	tphysvar fittedgoft=goft(logT,logrho,gofttab);
	string varname="goft";
	writearray(fittedgoft,varname);
	varname="logt";
	writearray(logT,varname);
	varname="logrho";
	writearray(logrho,varname);

	emission.setvar(0,fittedgoft);

	tphysvar fittedwidth=linefwhm(T,lambda0,atweight);

	emission.setvar(1,fittedwidth);

	return emission;
};

