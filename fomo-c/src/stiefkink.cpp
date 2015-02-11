#define WORKTAG 1001
#define RESULTTAG 1201

#include "header.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h> /*clock_t, clock, CLOCKS_PER_SEC*/
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
        clock_t t0;
	commrank = 0;
	commsize = 1;
	getarg(argc,argv);  // Read in arguments from call to main function. For an overview of different input parameters, run program with --help argument.
//	writefile();

	// Variables used for parallellisation using mpi (as per 26/09/14 not yet implemented)
//        int workheight = y_pixel;
//        if (commsize>1) workheight = 1;
//        int maxsize = x_pixel*workheight*lambda_pixel;
//        double *results;
//        results = (double *)malloc(maxsize*sizeof(double));//only used in reuse option
        

	double globalmax, globalmin;
	globalmax = 0.; globalmin = 0.;
	// Create G(T) interpolated cube or artificial images (depending on --reuse option) 
        //const int nframes=tend-tstart+1; //Number of time steps = number of simulation snapshots
	// double pi=4*atan(1.);
	vector<double> angles={M_PI/6.}; // M_PI/2.,M_PI/3.,M_PI/4. Rotation angles around y-axis
	int nangles=angles.size();
	int ng = eqx*eqy*eqz;
	int nvars = 5; // \rho, T, vx, vy, vz

        if (lambda_pixel==1) cout << "This code is configured to do imaging modelling!" << endl;
        if (lambda_pixel>1) cout << "This code is configured to do spectral modelling!" << endl;
	stringstream ss;
	Delaunay_triangulation_3 DT;
       // if not dynamicDT the triangulation can be done only once, however, it is saved for every snapshot for restart!
        bool dynamicDT=false;
        bool doDT=true;
        bool readDT=true;
	for (int t=tstart; t<=tend; t=t+tstep)
	{
		cout << endl << "Doing timestep " << t << endl << flush;
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
                if (t==tstart)
                  {
                  doDT=true;
                  readDT=true;
                   }
		if (reuse!=1)
		{

// Initialise datacube (member of class 'cube', see header.h) and set
			cube datacube(nvars,ng);
			EqType qtype=stiefkink;
			datacube.settype(qtype);
                        t0=clock();
            		if (commrank == 0) cout << "Reading in snapshot " << filename << "... " << flush;
			reademissioncube(datacube, filename, &DT);
			if (commrank == 0) cout << "Done!" << endl << flush;
                        cout<<"Reading time used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
                        t0=clock();
			datacube.fillcube();
			goftcube=emissionfromdatacube(datacube);
                        cout<<"Emission calculation used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
                      //steady grid and first time run
                      if(dynamicDT || doDT) // only do DT for dynamic grid or first time
                         { doDT=false;// never repeat
                            t0=clock();
                            DT=triangulationfromdatacube(datacube);// only use the grid info
                            cout<<"Triangulation used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
                         }
		      if (emissionsave.compare("none")!=0)
		 	{
		     	if (commrank==0) cout << "Writing data for reuse in file " << snapsave << "... " << flush;
                        t0=clock();
			writeemissioncube(goftcube,snapsave,&DT);                   
			if (commrank==0) cout << "Done!" << endl << flush;
                        cout<<"Writing data used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
			}


		}
		else // reuse part 
		{       t0=clock();
			if (commrank==0) cout << "Reuse emission data " << snapsave << "... " << flush;
			if (dynamicDT || readDT) // dynamic grid or first time read [DY 19 Nov 2014]
			{       readDT=false; // never repeat for steady grid 
				reademissioncube(goftcube, snapsave, &DT);
                                cout<<"Reading data and DT used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;

			}
			else // not dynamic grid and not first time read
			{
				reademissioncube(goftcube, snapsave);
                               cout<<"Reading data without DT used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
			}



// Here ends distinction between constant and time-varying grid

		if (commrank==0) cout << "Done!" << endl << flush;
                int workheight = y_pixel;
               if (commsize>1) workheight = 1;
                int maxsize = x_pixel*workheight*lambda_pixel;
                double *results;
                results = (double *)malloc(maxsize*sizeof(double));//only used in reuse option
	
		for (int i=0;i<nangles;i++)
		{       t0=clock();
			l=0;
			b=angles[i];
			for (int j=0; j<maxsize; j++) results[j]=0;
			mpi_calculatemypart(results,0,x_pixel-1,0,y_pixel-1,t,goftcube,&DT);
			cube observ(1,x_pixel*y_pixel*lambda_pixel);
			fillccd(observ,results,0,x_pixel-1,0,y_pixel-1);
                        cout<<"Integration used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
			// write out observ cube
			ss.str("");
			// ss << "/users/cpa/dyuan/fomodata/testemislos";
                        ss << outfileini; // output file initials
			ss << setfill('0') << setw(3) << int((angles[i]*180./M_PI)+.5);
			ss << "t";
			ss << setfill('0') << setw(3) << t;
			string outfile = ss.str();
			ss.str("");
                        t0=clock();
			cout << "Writing out integratted emission into " << outfile << "... " << flush;
			writeemissioncube(observ,outfile);
			cout << "Done!" << endl << flush;
                        cout<<"Writing emission data used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;

// Experiment: try to add png images as in the example of the torus.

//			double **image;
//			image = (double **)malloc(y_pixel*sizeof(double *));
//       		for (int row = 0; row < y_pixel; row++)
//		        {
//                		image[row] = (double *)malloc(x_pixel*sizeof(double));
//		// initialize image to black
//				for (int j=0; j<x_pixel; j++)
//						{
//					image[row][j]=0;
//		        }
//						}
//			tgrid grid=observ.readgrid();
//			tphysvar intens=observ.readvar(0);
//			int row, column;
//			for (int k=0; k<x_pixel*y_pixel*lambda_pixel; k++)
//			{
//				row=int(grid[1][k]);
//				column=int(grid[0][k]);
//				image[row][column]+=intens[k];
//			}
//			if ((x_pixel>=42)&&(y_pixel>=42)) writetime(image,t);
//			int i,j;
//			double max = findmax(image,&i,&j);
//			double min = findmin(image,&i,&j);
//			if (max>globalmax) globalmax=max;
//			if (min<globalmin) globalmin=min;
//#ifdef HAVEPNG
//			if (png) {
//			  cout << "Writing png file... " << flush;
//			  writepng(image,t);
//			  cout << "Done!" << endl << flush;
//			}
//#endif
		}		// end angles
		}		// end reuse=1
	}			// end t
	if (commrank==0) cout << "Gelukt\n" ;
	return EXIT_SUCCESS;
}
