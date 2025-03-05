#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

/**
 * @file 
 * This file contains a basic example for using FoMo. The example program takes as argument a file name.
 * The file name in the argument should contain a list of 3D coordinates with associated variables (rho, T, vx, vy, vz).
 */

using namespace std;

int main(int argc, char* argv[])
{
	// First check to see if a filename is given as argument.
	if (argc == 1)
	{
		cout << "This program should be run with one file name as argument." << endl;
		cout << "The text file should contain (x,y,z,rho,T,vx,vy,vz) as elements of 8 columns." << endl ;
		cout << "The file included in example/testfile.txt is provided as an example." << endl;
		exit(EXIT_FAILURE);
	}
	
	/// [Initialize FoMo]
	// Initialize the FoMo object
	FoMo::FoMoObject Object;
	/// [Initialize FoMo]
	
	/// [Read in data]
	// Open the file in the argument as fstream
	// File contains lines with x, y, z (Mm), rho (cm^-3), T (K), vx, vy, vz (m/s)
	ifstream filetoread(argv[1]);
	
	if (filetoread.is_open())
	{
		double tmpvar;
		while ( !filetoread.eof() )
		{
			vector<double> coordinates;
			// The points should be three dimensional, and thus we iterate until 3.
			for (unsigned int i=0; i<3; i++)
			{
				// Read one value of a coordinate at a time.
				filetoread >> tmpvar;
				// push it back into the coordinate vector.
				coordinates.push_back(tmpvar);
			}
			vector<double> variables;
			// There should be 5 variables (n, T, vx, vy, vz), and iterate until 5.
			for (unsigned int i=0; i<5; i++)
			{
				// Read one value of a variable at a time.
				filetoread >> tmpvar;
				// push it back into the variable vector.
				variables.push_back(tmpvar);
			}
			
			// If the file has not reached the end, this is probably a valid data point. Push it back 
			// into the FoMoObject.
			if ( !filetoread.eof() ) Object.push_back_datapoint(coordinates, variables);
		}
	}
	/// [Read in data]
	
	// data is in structure, now start the rendering

	/// [Set rendering options]
	Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat"); // the default value is "../chiantitables/goft_table_fe_12_0194small_abco.dat"
	Object.setabundfile("/empty"); //use "/empty" or do not set it at all for the default sun_coronal_2012_schmelz.abund file
	Object.setrendermethod("Thomson"); // NearestNeighbour is the default rendermethod
	Object.setobservationtype(FoMo::Imaging);
	// adjust the resolution with these parameters
	int x_pixel=1024;
	int y_pixel=1024;
	int z_pixel=1024;
	int lambda_pixel=1000;
	double lambda_width=1000; // in m/s
	Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	// determine where the output will be written
	Object.setoutfile("fomo-example-out-Th.");
	/// [Set rendering options]
	
	/// [Render]
	Object.render(0.0,1.75);
//(1.553,3.246);
//(6.25,0.101);
//0.0149,1.469);
	// alternatively, you could add different angles in radians as argument, e.g. Object.render(1.57/2.,1.57/2.) to obtain some nice Doppler shifts.
	/// [Render]
	
	/// [Details]
	// read the rendering
	FoMo::RenderCube rendercube=Object.readrendering();
	// read the resolution
	int nx,ny,nz,nlambda;
	double lambdawidth;
	rendercube.readresolution(nx,ny,nz,nlambda,lambdawidth);
	cout << "The rendering has resolution of x and y " << nx << " " << ny << ".\n";
	cout << "and " << nlambda << " pixels in the wavelength in a window of width " << lambdawidth << "m/s" << endl;
	cout << "It was done using the rendermethod " << rendercube.readrendermethod() << endl;
	cout << "with a resolution of " << nz << " along the line-of-sight." << endl;
	// This writes out the rendering results to the file fomo-output.dat.gz.
	// These options could be set before Object.render(), too. Then the output would have been written to fomo-example-out.l90b90.dat.gz
	// This option turns on binary writeout.
	rendercube.setwriteoutbinary();
	// This says that it needs to be zipped.
	rendercube.setwriteoutzip();
	// This says that only the zipped file needs to be kept. 
	rendercube.setwriteoutdeletefiles();
	rendercube.writegoftcube("fomo-output.txt");
	/// [Details]
}
