#include "FoMo.h"
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
}