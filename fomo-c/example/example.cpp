#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

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
	// File contains lines with x, y, z (km), rho (cm^-3), T (K), vx, vy, vz (m/s)
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
	Object.setchiantifile("/home/tom/data/idl/FoMo/chiantitables/goft_table_fe_12_0194small.dat");
	Object.setabundfile("/empty"); //use "/empty" for the default sun_coronal.abund file
	Object.setrendermethod("CGAL"); // CGAL is the default rendermethod
	Object.setobservationtype(FoMo::Spectroscopic);
	int x_pixel=101;
	int y_pixel=102;
	int z_pixel=300;
	int lambda_pixel=30;
	double lambda_width=.13;
	Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	Object.setoutfile("fomo-example-out.");
	/// [Set rendering options]
	
	/// [Render]
	Object.render();
	// alternatively, you could add angles as in radians as argument, e.g. Object.render(1.57,0.52).
	/// [Render]
	
	/// [Details]
	FoMo::RenderCube rendercube=Object.readrendering();
	int nx,ny,nz,nlambda;
	double lambdawidth;
	rendercube.readresolution(nx,ny,nz,nlambda,lambdawidth);
	cout << "The rendering has resolution of x and y " << nx << " " << ny << ".\n";
	cout << "and " << nlambda << " pixels in the wavelength." << endl;
	cout << "It was done using the rendermethod " << rendercube.readrendermethod() << endl;
	cout << "with a resolution of " << nz << " along the line-of-sight." << endl;
	// This writes out the rendering results to the file fomo-output.txt.
	rendercube.writegoftcube("fomo-output.txt");
	/// [Details]
}