#define WORKTAG 1001
#define RESULTTAG 1201

#include "header.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

/* Example of main routine to generate synthetic observations from irregularly gridded input data.
This program is adapted to process data written into a text file as follows:
--------------------------------------------------------------------------------
dimension
"intqtype": large integer (although unused), we override 'Eqtype' directly in main function
number of grid points
number of variables
{x, y, z, ne, T, vx, vy, vz} for every grid point
--------------------------------------------------------------------------------
For more information about the underlying geometric operations, please visit the FoMo-C Wiki on https://wiki.esat.kuleuven.be/FoMo/FoMo-C
*/

int main(int argc, char* argv[])
{

// Initialize global paramters, get arguments and set physical parameters
	int commrank,commsize; // For parallellisation

	commrank = 0;
	commsize = 1;
	getarg(argc,argv);  // Read in arguments from call to main function. For an overview of different input parameters, run program with --help argument.
	writefile();

	// Variables used for parallellisation using mpi (as per 26/09/14 not yet implemented)
	int workheight = y_pixel;
	if (commsize>1)	workheight = 1;
	int maxsize = x_pixel*workheight*lambda_pixel;
	double *results;
	// Results is a one dimensional array with all data from rectangle [x1,x2]*[y1,y2]
	// borders included!!!
	results = (double *)malloc(maxsize*sizeof(double));
	double globalmax, globalmin;
	globalmax = 0.; globalmin = 0.;

	// Create G(T) interpolated cube or artificial images (depending on --reuse option) 
        //const int nframes=tend-tstart+1; //Number of time steps = number of simulation snapshots
	double pi=4*atan(1.);
	vector<double> angles={pi/2.,pi*2./3.,pi/3.,pi/6.}; // Rotation angles around y-axis
	int nangles=angles.size();
	int ng = eqx*eqy*eqz;
	int nvars = 5; // \rho, T, vx, vy, vz

	double lambda0=readgoftfromchianti(chiantifile);
	if (lambda0 > 500.) lambda_width=.6;
	stringstream ss;
	Delaunay_triangulation_3 DT;
	for (int t=tstart; t<=tend; t=t+tstep)
	{
		cout << endl << "Doing timestep " << t << endl << flush;
	//	ss << "/users/cpa/dyuan/fomodata/stslow";
	        ss << infileini; // argument input
         	ss << setfill('0') << setw(3) << t;
		ss << ".dat";
		string filename=ss.str(); // input file
		ss.str("");// set to null
		ss << emissionsave;
		ss << setfill('0') << setw(3) << t;
		string snapsave=ss.str();// save file
		ss.str("");
		cube goftcube(1,1,1);
		if (reuse!=1)
		{

// Initialise datacube (member of class 'cube', see header.h) and set
			cube datacube(nvars,ng);
			EqType qtype=stiefkink;
			datacube.settype(qtype);
				if (commrank == 0) cout << "Reading in snapshot " << filename << "... " << flush;
				reademissioncube(datacube, filename, &DT);
				if (commrank == 0) cout << "Done!" << endl << flush;

			datacube.fillcube();
			goftcube=emissionfromdatacube(datacube);

//		In this part the emission data cubes are written to a file, they can be restored to view the emission from different angles with respect to the loop axis. For this purpose the triangulation of the original data are stored together with the emission.
// The uncommented lines below assume an irregular grid that changes with time. Uncomment the following lines for grids that are constant for each time step (this saves time; same triangulation at each time step).
//			{

				DT=triangulationfromdatacube(datacube);
				if (emissionsave.compare("none")!=0)
				{
					if (commrank==0) cout << "Writing data for reuse in file " << snapsave << "... " << flush;
					writeemissioncube(goftcube,snapsave,&DT);
					if (commrank==0) cout << "Done!" << endl << flush;
				}


//			}
//			else	// 
//			{
//				if (emissionsave.compare("none")!=0)
//				{
//					if (commrank==0) cout << "Writing data for reuse in file " << snapsave << "... " << flush;
//					writeemissioncube(goftcube,snapsave);
//					if (commrank==0) cout << "Done!" << endl << flush;
//				}
//			}
		}
		else
		{
			if (commrank==0) cout << "Reading in data from previous run in file " << snapsave << "... " << flush;
//			if (t == 0) 
//			{
				reademissioncube(goftcube, snapsave, &DT);
//			}
//			else
//			{
//				reademissioncube(goftcube, snapsave);
//			}


// Here ends distinction between constant and time-varying grid

			if (commrank==0) cout << "Done!" << endl << flush;
		
		for (int i=0;i<nangles;i++)
		{
			l=0;
			b=angles[i];
			for (int j=0; j<maxsize; j++) results[j]=0;
			mpi_calculatemypart(results,0,x_pixel-1,0,y_pixel-1,t,goftcube,&DT);
			cube observ(1,x_pixel*y_pixel*lambda_pixel);
			fillccd(observ,results,0,x_pixel-1,0,y_pixel-1);
			// write out observ cube
			ss.str("");
			// ss << "/users/cpa/dyuan/fomodata/testemislos";
                        ss << outfileini; // output file initials
			ss << setfill('0') << setw(3) << int((angles[i]*180./pi)+.5);
			ss << "t";
			ss << setfill('0') << setw(3) << t;
			string outfile = ss.str();
			ss.str("");
			cout << "Writing out computed emission to " << outfile << "... " << flush;
			writeemissioncube(observ,outfile);
			cout << "Done!" << endl << flush;

		// Experiment: try to add png images as in the example of the torus.

			double **image;
			image = (double **)malloc(y_pixel*sizeof(double *));
        		for (int row = 0; row < y_pixel; row++)
		        {
                		image[row] = (double *)malloc(x_pixel*sizeof(double));
		// initialize image to black
				for (int j=0; j<x_pixel; j++)
						{
					image[row][j]=0;
		        }
						}
			tgrid grid=observ.readgrid();
			tphysvar intens=observ.readvar(0);
			int row, column;
			for (int k=0; k<x_pixel*y_pixel*lambda_pixel; k++)
			{
				row=int(grid[1][k]);
				column=int(grid[0][k]);
				image[row][column]+=intens[k];
			}
			if ((x_pixel>=42)&&(y_pixel>=42)) writetime(image,t);
			int i,j;
			double max = findmax(image,&i,&j);
			double min = findmin(image,&i,&j);
			if (max>globalmax) globalmax=max;
			if (min<globalmin) globalmin=min;
#ifdef HAVEPNG
			if (png) {
			  cout << "Writing png file... " << flush;
			  writepng(image,t);
			  cout << "Done!" << endl << flush;
			}
#endif
		}		// end angles
		}		// end reuse=1
	}			// end t
	if (commrank==0) cout << "Gelukt\n" ;
	return EXIT_SUCCESS;
}
