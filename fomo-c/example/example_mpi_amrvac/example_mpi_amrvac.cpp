#include "FoMo.h"
#include "FoMo-amrvac.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>

using namespace std;

string amrvac_version("old"); // provide the version of the AMRVAC code. It can be "old" or "gitlab".

int gamma_eqparposition = 0; // provide the index of the eqpar vector in which the gamma value is contained

double L_unit = 1e8; // length normalisation value
double L_unitFoMo = 1e6; // (express length unit in Mm for FoMo)

double rho_unit = 1e-12;
double Teunit = 1e+6;
double n_unit = rho_unit * 1.204 * 1.e21;

int main(int argc, char* argv[])
{
	// First check to see if two filenames are given as argument.
	if (argc != 3)
	{
		cout << "This program should be run with two file names as argument." << endl;
		cout << "The first file is the AMRVAC .dat file to be rendered, the second one is the amrvac.par file." << endl;
		exit(EXIT_FAILURE);
	}
	
	// Initialize the FoMo object and read data from .dat file and .par file
		FoMo::FoMoObject Object = read_amrvac_dat_file(argv[1],argv[2],amrvac_version,gamma_eqparposition, n_unit, Teunit, L_unit/L_unitFoMo);
	
	// data is in structure, now start the rendering
	
	/// [Set rendering options]
	Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat"); // the default value is "../chiantitables/goft_table_fe_12_0194_abco.dat"
	Object.setabundfile("/empty"); //use "/empty" or do not set it at all for the default sun_coronal_2012_schmelz.abund file
	Object.setrendermethod("NearestNeighbour"); // NearestNeighbour is the default rendermethod
	Object.setobservationtype(FoMo::Spectroscopic);
	// adjust the resolution with these parameters
	int x_pixel=149;
	int y_pixel=64;
	int z_pixel=500;
	int lambda_pixel=100;
	double lambda_width=200000; // in m/s
	Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	// determine where the output will be written
	Object.setoutfile("fomo-example-out.");
	
	Object.setwriteoutbinary();
	Object.setwriteoutzip();
	Object.setwriteoutdeletefiles();
	/// [Set rendering options]
	
	/// [Render]
	Object.render(2*atan(1.),4*atan(1.)/6.);
	/// [Render]
}
