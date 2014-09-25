#define WORKTAG 1001
#define RESULTTAG 1201

#include "header.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[])
{

// Initialize global paramters, get arguments and set physical parameters
	int commrank,commsize;

	commrank = 0;
	commsize = 1;
	getarg(argc,argv);  
// Read in arguments from call to main function. For an overview of different input parameters, run program with --help argument.

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
	const int nframes=16; //Number of time steps = number of simulation snapshots
	double pi=4*atan(1.);
	vector<double> angles={pi/2.,pi/3.,pi/6.};
	int nangles=angles.size();
	int ng = eqx*eqy*eqz;
	int nvars = 5; // \rho, T, vx, vy, vz

	double lambda0=readgoftfromchianti(chiantifile);
	if (lambda0 > 500.) lambda_width=.6;

	stringstream ss;
	Delaunay_triangulation_3 DT;
	for (int t=0; t<nframes; t++)
	{
		cout << endl << "Doing timestep " << t << endl << flush;
		//ss << "patrickfiles/datcubes_ka2.24_";
		ss << "/users/cpa/sgijsen/FoMo/eigft/";
		ss << setfill('0') << setw(3) << t;
		ss << ".dat";
		string filename=ss.str();
		ss.str("");
		ss << emissionsave;
		ss << setfill('0') << setw(3) << t;
		string snapsave=ss.str();
		ss.str("");
		cube goftcube(1,1,1);
		if (reuse!=1)
		{
			cube datacube(nvars,ng);
			EqType qtype=stiefkink;
			datacube.settype(qtype);
			if (datacube.readtype() == stiefkink) 
			{
				if (commrank == 0) cout << "Reading in snapshot " << filename << "... " << flush;
				reademissioncube(datacube, filename, &DT);
				if (commrank == 0) cout << "Done!" << endl << flush;
			}
			datacube.fillcube();
			goftcube=emissionfromdatacube(datacube);
//			The uncommented lines below assume an irregular grid that changes with time. Uncomment the following lines for grids that are constant for each time step (this saves time; same triangulation at each time step).
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


// 

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
			ss << "/users/cpa/sgijsen/tmp/fomo-c/tmpadv/data_observ";
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
