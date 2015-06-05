#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

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
	
	// Initialize the FoMo object
	FoMo::FoMoObject Object;
	
	// Open the file in the argument as fstream
	// File contains lines with x, y, z, rho (cm^-3), T (K), vx, vy, vz (m/s)
	ifstream filetoread(argv[1]);
	
	if (filetoread.is_open())
	{
		double tmpvar;
		while ( !filetoread.eof() )
		{
			vector<double> coordinates;
			for (unsigned int i=0; i<3; i++)
			{
				filetoread >> tmpvar;
				coordinates.push_back(tmpvar);
			}
			vector<double> variables;
			for (unsigned int i=0; i<5; i++)
			{
				filetoread >> tmpvar;
				variables.push_back(tmpvar);
			}
			
			if ( !filetoread.eof() ) Object.datacube.push_back(coordinates, variables);
		}
	}
	
	// data is in structure, now start the rendering
	Object.setchiantifile("/home/tom/data/idl/FoMo/chiantitables/goft_table_fe_12_0194small.dat");
	Object.setabundfile("/empty"); //use "/empty" for the default sun_coronal.abund file
	Object.rendering.setrendermethod("CGAL"); // CGAL is the default rendermethod
	Object.rendering.setobservationtype(FoMo::Spectroscopic);
	int x_pixel=101;
	int y_pixel=102;
	int z_pixel=300;
	int lambda_pixel=30;
	double lambda_width=.13;
	Object.rendering.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	Object.render();
	
	FoMo::RenderCube rendercube=Object.readrendering();
	int nx,ny,nz,nlambda;
	double lambdawidth;
	rendercube.readresolution(nx,ny,nz,nlambda,lambdawidth);
	cout << "The rendering has resolution of x and y " << nx << " " << ny << ".\n";
	cout << "and " << nlambda << " pixels in the wavelength." << endl;
	cout << "It was done using the rendermethod " << rendercube.readrendermethod() << endl;
	cout << "with a resolution of " << nz << " along the line-of-sight." << endl;
	rendercube.writegoftcube("fomo-output.txt");
}