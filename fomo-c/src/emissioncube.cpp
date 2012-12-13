#include "header.h"
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_const_mksa.h>

// This should become part of the input options, or based on the choice of line position.
const char* chiantifile="chiantitables/goft_table_f2rt_171.dat";

// Physical constants
const double alphaconst=1.0645; // alpha=\sqrt{2\pi}/2/\sqrt{2ln2}
const double kb=GSL_CONST_MKSA_BOLTZMANN; // Boltzmann constant
const double mh=GSL_CONST_MKSA_MASS_PROTON; // hydrogen mass
const double c=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light


tphysvar goft(const tphysvar T, const tphysvar logrho, const cube gofttab)
{
	// how to efficiently do this? I guess we do a Delaunay as well?
	tphysvar g=T;
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
	int nvars=2; // G(T) and width
	cube gofttab(nvars,ntnrho,2); //only 2 dimensions
	gofttab.settype(gofttable);
	
	tgrid tempgrid=new tcoord[2];
	// Warning this assumes that there is a straightforward mapping between tcoord and tphysvar!
	tempgrid[0]=tempt;
	tempgrid[1]=temprho;
	
	gofttab.setgrid(tempgrid);
	gofttab.setvar(0,tempgoft);

	// This would be the place to read in the width of the lines

	return gofttab;
}

tphysvar linefwhm(const tphysvar T, const double lambda0, const double atweight)
{
	tphysvar w=2*sqrt(2*log(2))*sqrt(kb/mh/atweight)*lambda0/c*sqrt(T);
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
	
	string ion="";
	double lambda0=.0;
	double atweight=1;
	cube gofttab=readgoftfromchianti(chiantifile,ion,lambda0,atweight);

	tphysvar fittedgoft=goft(T,logrho,gofttab);

	emission.setvar(0,fittedgoft);

	tphysvar fittedwidth=linefwhm(T,lambda0,atweight);

	emission.setvar(1,fittedwidth);

	return emission;
};

