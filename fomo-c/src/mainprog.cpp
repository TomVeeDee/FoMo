#define WORKTAG 1001
#define RESULTTAG 1201

#include "header.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <boost/progress.hpp>

int main(int argc, char* argv[])
{
// Initialize global paramters, get arguments and set physical parameters
#ifdef HAVEMPI
	MPI_Init(&argc,&argv);
#endif
	int commrank,commsize;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
        MPI_Comm_size(MPI_COMM_WORLD,&commsize);
#else
	commrank = 0;
	commsize = 1;
#endif
	getarg(argc,argv);
	writeparameters(cout,'v');
//	writefile(); commented by DY. 11 Feb 2015

#ifdef HAVEMPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	// here starts mpi
	int workheight = y_pixel;
	if (commsize>1)	workheight = 1;
	int maxsize = x_pixel*workheight*lambda_pixel;
	double *results;
	// results is a one dimensional array with all data from rectangle [x1,x2]*[y1,y2]
	// borders included!!!
	results = (double *)malloc(maxsize*sizeof(double));
	cube observ(1,x_pixel*y_pixel*lambda_pixel);

	// if we want a movie, iteration over time t starts here
	char imagefile[] = "./imagefile.XXXXXX";
	mkstemp(imagefile);
	double globalmax, globalmin;
	globalmax = 0.; globalmin = 0.;
	{ // only ofstream s in this scope
	ofstream s(imagefile);
	// Let's first get the code working for a single frame. Later on, we can still change to more frames, although that 
	const int nframes=1;
	int ng = eqx*eqy*eqz;
	int nvars = 5; // \rho, T, vx, vy, vz
	double lambda0=readgoftfromchianti(chiantifile);
	if (lambda0 > 500.) lambda_width=.6;
	for (int t=0.; t<nframes; t++)
	{
		cube goftcube(1,1,1);
		Delaunay_triangulation_3 DT;
		if (reuse!=1) 
		{
			cube datacube(nvars,ng);
			EqType qtype=builtineq;
			datacube.settype(qtype);
			datacube.fillcube();
			goftcube=emissionfromdatacube(datacube);
			DT=triangulationfromdatacube(datacube);
			if (emissionsave.compare("none")!=0)
			{
				if (commrank==0) cout << "Writing data for reuse in file " << emissionsave << "... " << flush;
				writeemissioncube(goftcube,emissionsave,&DT);
				if (commrank==0) cout << "Done!" << endl << flush;
			}
		}
		else
		{
			if (commrank==0) cout << "Reading in data from previous run in file " << emissionsave << "... " << flush;
			reademissioncube(goftcube, emissionsave, &DT);
			// Perform some checks on the Delaunay triangulation and data cube that have been read in
			// cout << DT.number_of_vertices() << goftcube.readdim() << endl << flush;
			if (commrank==0) cout << "Done!" << endl << flush;
		}

		boost::progress_display show_progress(y_pixel);
		if (commsize>1)
		{
#ifdef HAVEMPI
			int done = 0;
			int coords[commsize][4],x1,x2,y1,y2;
			for (int i=0; i<commsize; i++)
				for (int j=0; j<4; j++)
					coords[i][j]=0;
			MPI_Status status;
			if (commrank==0) //we're the master, let's give the slaves some work
			{
				for (int r=1; r<commsize; r++)
				{
					coords[r-1][0]=0;
					coords[r-1][1]=x_pixel-1;
					coords[r-1][2]=done;
					if (done+workheight-1<y_pixel) coords[r-1][3]=done+workheight-1;
					else coords[r-1][3]=y_pixel-1;
					if (coords[r-1][2]<y_pixel)
					{
						MPI_Send(coords[r-1],4,MPI_INT,r,WORKTAG,MPI_COMM_WORLD);
						done+=workheight;
						//cout << "Sent job to node " << r << endl;
					}
				}
				while (done<y_pixel)
				{
					int r;
					// receive the results
					MPI_Recv(results,maxsize,MPI_DOUBLE,MPI_ANY_SOURCE,RESULTTAG,MPI_COMM_WORLD,&status);
					r = status.MPI_SOURCE;
					x1 = coords[r-1][0];
					x2 = coords[r-1][1];
					y1 = coords[r-1][2];
					y2 = coords[r-1][3];
					// send new work to the node
					coords[r-1][0]=0;
					coords[r-1][1]=x_pixel-1;
					coords[r-1][2]=done;
					if (done+workheight-1<y_pixel) coords[r-1][3]=done+workheight-1;
					else coords[r-1][3]=y_pixel-1;
					if (coords[r-1][2]<y_pixel)
					{
						MPI_Send(coords[r-1],4,MPI_INT,r,WORKTAG,MPI_COMM_WORLD);
						done+=workheight;
						//cout << "Sent job to node " << r << endl;
					}
					// process the received results
					fillccd(observ,results,x1,x2,y1,y2);
					++show_progress;
				}
				for (int r=1; r<commsize; r++)
				{
					// receive the results
					MPI_Recv(results,maxsize,MPI_DOUBLE,MPI_ANY_SOURCE,RESULTTAG,MPI_COMM_WORLD,&status);
					r = status.MPI_SOURCE;
					x1 = coords[r-1][0];
					x2 = coords[r-1][1];
					y1 = coords[r-1][2];
					y2 = coords[r-1][3];
					coords[0][0]=-10;
					MPI_Send(coords[r-1],4,MPI_INT,r,WORKTAG,MPI_COMM_WORLD);
					// process the received results
					fillccd(observ,results,x1,x2,y1,y2);
				}
			}
			else //oh no, we're a slave, work to do
			{
				while (coords[0][0]>=0)
				{
					MPI_Recv(coords[0],4,MPI_INT,0,WORKTAG,MPI_COMM_WORLD,&status);
					//cout << "received job" << endl;
					if (coords[0][0]>=0) 
					{
					for (int i=0; i<maxsize; i++) results[i]=0.;
					mpi_calculatemypart(results,coords[0][0],coords[0][1],coords[0][2],coords[0][3],t,goftcube,&DT);
					MPI_Send(results,maxsize,MPI_DOUBLE,0,RESULTTAG,MPI_COMM_WORLD);
					}
				}
			}
#endif
		}
		else
		{
			mpi_calculatemypart(results,0,x_pixel-1,0,y_pixel-1,t,goftcube,&DT);
			fillccd(observ,results,0,x_pixel-1,0,y_pixel-1);
		}
		if (commrank==0)
		{
			// write out observ cube
			stringstream ss;
			ss << "fomo-c" << ".l" << l << ".b" << b << ".observ";
			writeemissioncube(observ,ss.str());
		// allocate image
			double **image;
			image = (double **)malloc(y_pixel*sizeof(double *));
        		for (int row = 0; row < y_pixel; row++)
		        {
                		image[row] = (double *)malloc(x_pixel*sizeof(double));
		// initialize image to black
				for (int j=0; j<x_pixel; j++)
					image[row][j]=0;
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
			// writeout image in outputfile.t*.png
			if (png) writepng(image,t);
#endif
			if (warray) writearray(image,t);
			// write image to disk in temporary file
			for (int j=0; j<y_pixel; j++)
				for (int k=0; k<x_pixel; k++)
					s << image[j][k] << endl;
		}
#ifdef HAVEMPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	} // ofstream s destroyed
#ifdef HAVEMPEG
	stringstream ss;
	ss << "outputfile.t" << seqnum << ".l" << l << ".b" << b << ".png";
	string graphic_out=ss.str();
	if (commrank==0)
	{
		// write movie
		if (mpeg) writemovie(imagefile, globalmin, globalmax);
	}
#endif
	unlink(imagefile);
#ifdef HAVEMPI
	// Finalize mpi
	MPI_Finalize();
#endif
	if (commrank==0) cout << "Gelukt\n" ;
	return EXIT_SUCCESS;
}
