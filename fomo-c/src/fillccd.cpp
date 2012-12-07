#include "header.h"
#include <cmath>
#include <cstdlib>


const int y_pixel = 200;
const int x_pixel = 201;
//pixelsize in kilometers, corresponds to 0.5 arcsecond
// in contradiction with Schrijver et al. (99)? They claim a spatial resolution of 1,25 arcseconds, which corresponds to 900km.
const double size_pixel = 727.22;
// z_pixel is need to determine the resolution in the integration over the line-of-sight, ideally it should be as big as possible
// for 1024, it takes one hour to calculate one frame on one P4 3GHz CPU
const int z_aspect=2; // determines the aspect ratio between the z direction and the other directions
const int z_pixel=x_pixel*z_aspect; // take a factor aspect to make 3D pixel a cube
const double size_z_pixel=size_pixel*x_pixel/z_pixel*z_aspect; 
// this is the size of a pixel in the z direction, the factor aspect expresses the ratio
// size of the considered plasma cube in the z direction versus the size in x and y

double findmax(double * const * const ccd, int * i, int * j)
// returns the maximumvalue of the ccd
// i contains the vertical coordinate
// j contains the horizontal coordinate
{
	*i=0;
	*j=0;
	for (int k=0; k<y_pixel; k++)
		for (int l=0; l<x_pixel; l++)
	{
		if (ccd[k][l]>ccd[*i][*j])
		{
			*i = k;
			*j = l;
		}
	}
	return ccd[*i][*j];
}

double findmin(double * const * const ccd, int * i, int * j)
// returns the minimumvalue of the ccd
// i contains the vertical coordinate
// j contains the horizontal coordinate
{
	*i=0;
	*j=0;
	for (int k=0; k<y_pixel; k++)
		for (int l=0; l<x_pixel; l++)
	{
		if (ccd[k][l]<ccd[*i][*j])
		{
			*i = k;
			*j = l;
		}
	}
	return ccd[*i][*j];
}

int MPE_Decomp1d(int n, int size, int rank, int *s, int *e )
/*
~  This file contains a routine for producing a decomposition of a 1-d array
~  when given a number of processors.  It may be used in "direct" product
~  decomposition.  The values returned assume a "global" domain in [1:n]
~ */
/*@
~  MPE_Decomp1d - Compute a balanced decomposition of a 1-D array

~  Input Parameters:
+ n  - Length of the array
. size - Number of processors in decomposition
- rank - Rank of this processor in the decomposition (0 <= rank < size)

~  Output Parameters:
. s,e - Array indices are s:e, with the original array considered as 1:n.
@*/
{
	int nlocal, deficit;

	nlocal      = n / size;
	*s  = rank * nlocal + 1;
	deficit     = n % size;
	*s  = *s + ((rank < deficit) ? rank : deficit);
	if (rank < deficit) nlocal++;
	*e      = *s + nlocal - 1;
	if (*e > n || rank == size-1) *e = n;
#ifdef HAVEMPI
	return MPI_SUCCESS;
#else
	return EXIT_SUCCESS;
#endif
}

void mpi_getcoords(int & x1, int & x2, int & y1, int & y2)
{
// copied from Dries Kimpe
#ifdef HAVEMPI
	int Dims[2] = {0, 0};               // Choose the dimension yourself
	int Periods[2] = {0, 0};		// Non periodic domain
	int Coords[2];
	int CPUCount;

	MPI_Comm_size(MPI_COMM_WORLD, &CPUCount);
//	MPI_Comm_size(ParentComm, &CPUCount);

	// Do the split
	MPI_Dims_create (CPUCount, 2, Dims);

	// Now Dims[0] contains the number of CPUs in the X direction, Dims[1] in Y

	MPI_Comm CartComm;
	// Create the cartesian communicator
	MPI_Cart_create (MPI_COMM_WORLD, 2, Dims, Periods, 1, &CartComm);

	// Use CartComm from now on

	if (CartComm == MPI_COMM_NULL)
	{
		cerr << "CPU left out!\n";
		exit(EXIT_FAILURE);
	}

	// Returns the coordinates of this CPU in Coords
	MPI_Cart_get (CartComm, 2, Dims, Periods, Coords);
	// We don't need the Comm anymore
	MPI_Comm_free (&CartComm);

	int _LX1, _LX2, _LY1, _LY2;
	// Do the decomposition in the X direction
	MPE_Decomp1d (x_pixel, Dims[0], Coords[0], &_LX1, &_LX2);

	// Do the decomposition in the Y direction
	MPE_Decomp1d (y_pixel, Dims[1], Coords[1], &_LY1, &_LY2);

	// Fix coordinates (MPE_Decomp1d uses 1-based instead of 0-based)
	--_LX1; --_LX2; --_LY1; --_LY2;
	x1 = _LX1;
	x2 = _LX2;
	y1 = _LY1;
	y2 = _LY2;
#else
	x1 = 0;
	x2 = x_pixel-1;
	y1 = 0;
	y2 = y_pixel-1;
#endif
}

void mpi_calculatemypart(double* results, const int x1, const int x2, const int y1, const int y2, const double t, cube datacube)
{
//
// results is an array of at least dimension (x2-x1+1)*(y2-y1+1) and must be initialized to zero
// 
// determine contributions per pixel
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0) cout << "building frame: ";
	double R = length*1000./M_PI;
	cube goftcube=emissionfromdatacube(datacube);
	tgrid grid = datacube.readgrid();
	for (int i=y1; i<y2+1; i++)
		for (int j=x1; j<x2+1; j++)
			for (int k=0; k<z_pixel; k++) // scanning through ccd
	{
		double xacc = (j-x_pixel/2)*size_pixel;
		double yacc = (i-y_pixel/2)*size_pixel;
		double zacc = (k-z_pixel/2)*size_z_pixel; 
		const double psi=0;
		double x = (cos(psi)*cos(l)-sin(psi)*sin(l)*sin(b))*xacc+(-sin(psi)*cos(b))*yacc+(-cos(psi)*sin(l)-sin(psi)*cos(l)*sin(b))*zacc;
		double y = (sin(psi)*cos(l)+cos(psi)*sin(l)*sin(b))*xacc+(cos(psi)*cos(b))*yacc+(-sin(psi)*sin(l)+cos(psi)*cos(l)*sin(b))*zacc;
		double z = (sin(l)*cos(b))*xacc+(-sin(b))*yacc+(cos(l)*cos(b))*zacc;
		double r = sqrt(pow(sqrt(pow(x,2)+pow(z,2))-R,2)+pow(y,2));
		r/=width*1000./0.02; // normalize r to 0 -> 1
		double phi = atan(y/(sqrt(pow(x,2)+pow(z,2))-R));
		if (pow(x,2)+pow(z,2)<pow(R,2)) phi+=M_PI; // (x,z) lies inside the circle in the XZ-plane
		double z_or;
		if ((x>0.00001)||(x<-0.00001)) // x lies far enough from 0, no division problems
		{
			z_or = atan(z/x)/M_PI;
		}
		else
		{
			z_or = .5; // does not matter which value is taken here, because it is the singular point of the coordinate system
		}
		if (z_or<0) z_or++;  // atan returns a value between -pi/2 and pi/2
		if (z>=0.)
		{
			results[(i-y1)*(x2-x1+1)+j-x1]+=pow(density(r,phi,z_or),2);
		}
		// print progress
		if ((commrank==0)&&(j==x1)&&(k==0)) 
		{
			progressbar(i,y1,y2);
			//cout << "i" << i << "y1" << y1 << "y2" << y2;
		}
	}
	if (commrank==0) cout << " finished!" << endl;
}

void fillccd(double * const * const frame, const double *results, const int x1, const int x2, const int y1, const int y2)
{
	for (int i=y1; i<y2+1; i++)
		for (int j=x1; j<x2+1; j++)
			frame[i][j]=results[(i-y1)*(x2-x1+1)+j-x1];
}

