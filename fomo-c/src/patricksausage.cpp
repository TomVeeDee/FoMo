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
	writeparameters(cout,'v');
	writefile();

	// here starts mpi
	int workheight = y_pixel;
	if (commsize>1)	workheight = 1;
	int maxsize = x_pixel*workheight*lambda_pixel;
	double *results;
	// results is a one dimensional array with all data from rectangle [x1,x2]*[y1,y2]
	// borders included!!!
	results = (double *)malloc(maxsize*sizeof(double));

	// Let's first get the code working for a single frame. Later on, we can still change to more frames, although that 
	const int nframes=30;
	double pi=4*atan(1.);
	vector<double> angles={pi/2.,pi/3.,pi/4.,pi/6.};
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
		ss << "/volume1/scratch/tomvd/patrickdata/datcubes_ka2.24_";
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
			EqType qtype=patricksausage;
			datacube.settype(qtype);
			if (datacube.readtype() == patricksausage) 
			{
				if (commrank == 0) cout << "Reading in snapshot " << filename << "... " << flush;
				reademissioncube(datacube, filename, &DT);
				if (commrank == 0) cout << "Done!" << endl << flush;
			}
			datacube.fillcube();
			goftcube=emissionfromdatacube(datacube);
			if (t==0)
			{
				DT=triangulationfromdatacube(datacube);
				if (emissionsave.compare("none")!=0)
				{
					if (commrank==0) cout << "Writing data for reuse in file " << snapsave << "... " << flush;
					writeemissioncube(goftcube,snapsave,&DT);
					if (commrank==0) cout << "Done!" << endl << flush;
				}
			}
			else
			{
				if (emissionsave.compare("none")!=0)
				{
					if (commrank==0) cout << "Writing data for reuse in file " << snapsave << "... " << flush;
					writeemissioncube(goftcube,snapsave);
					if (commrank==0) cout << "Done!" << endl << flush;
				}
			}
		}
		else
		{
			if (commrank==0) cout << "Reading in data from previous run in file " << snapsave << "... " << flush;
			if (t == 0) 
			{
				reademissioncube(goftcube, snapsave, &DT);
			}
			else
			{
				reademissioncube(goftcube, snapsave);
			}
			// Perform some checks on the Delaunay triangulation and data cube that have been read in
			// cout << DT.number_of_vertices() << goftcube.readdim() << endl << flush;
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
			ss << "/volume1/scratch/tomvd/patrickdata/fomo-c.observ.b";
			ss << setfill('0') << setw(3) << int((angles[i]*180./pi)+.5);
			ss << "t";
			ss << setfill('0') << setw(3) << t;
			string outfile = ss.str();
			ss.str("");
			cout << "Writing out computed emission to " << outfile << "... " << flush;
			writeemissioncube(observ,outfile);
			cout << "Done!" << endl << flush;
		}
		}
	}
	if (commrank==0) cout << "Gelukt\n" ;
	return EXIT_SUCCESS;
}
