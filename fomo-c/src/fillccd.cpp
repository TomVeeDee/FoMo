#include "header.h"
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>

// CGAL
#include <CGAL/interpolation_functions.h>


const int y_pixel = 45;
const int x_pixel = 300;
const int z_aspect=2; // determines the aspect ratio between the z direction and the other directions
const int z_pixel=x_pixel*z_aspect; // take a factor aspect to make 3D pixel a cube
const int lambda_pixel = 60;
double lambda_width =0.3; // 0.14;
// if lambda0 is larger than 500, then lambda_width=0.3

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light

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

Delaunay_triangulation_3 triangulationfromdatacube(cube goftcube)
{
	typedef K::Point_3                                    Point;
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();
	vector<Point> delaunaygrid;
	delaunaygrid.resize(ng);

	for (int i=0; i<ng; i++)
	{
		delaunaygrid[i]=Point(grid[0][i],grid[1][i],grid[2][i]);
	}

	// compute the Delaunay triangulation
	if (commrank==0) cout << "Doing Delaunay triangulation for interpolation onto rays... " << flush;
	Delaunay_triangulation_3 DT;
	// The triangulation should go quicker if it is sorted
	// CGAL::spatial_sort(delaunaygrid.begin(),delaunaygrid.end());
	// but I don't know how this affects the values in the maps.
	// Apparently, this is already done internally.
	DT.insert(delaunaygrid.begin(),delaunaygrid.end());
	if (commrank==0) cout << "Done!" << endl << flush;
	return DT;
}

void mpi_calculatemypart(double* results, const int x1, const int x2, const int y1, const int y2, const double t, cube goftcube, Delaunay_triangulation_3* DTpointer)
{
//
// results is an array of at least dimension (x2-x1+1)*(y2-y1+1)*lambda_pixel and must be initialized to zero
// 
// determine contributions per pixel

	typedef K::FT                                         Coord_type;
	typedef K::Point_3                                    Point;
        std::map<Point, Coord_type, K::Less_xyz_3> peakmap, fwhmmap, losvelmap;
//        typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xyz_3 > >  Value_access;

	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) cout << "Rotating coordinates to POS reference... " << flush;
	vector<double> xacc, yacc, zacc;
	xacc.resize(ng);
	yacc.resize(ng);
	zacc.resize(ng);
	vector<double> gridpoint;
	gridpoint.resize(3);
	Point temporarygridpoint;
	// Define the unit vector along the line-of-sight
	vector<double> unit = {sin(b)*cos(l), -sin(b)*sin(l), cos(b)};
	// Read the physical variables
	tphysvar peakvec=goftcube.readvar(0);
	tphysvar fwhmvec=goftcube.readvar(1);
	tphysvar vx=goftcube.readvar(2);
	tphysvar vy=goftcube.readvar(3);
	tphysvar vz=goftcube.readvar(4);
	double losvelval;

// No openmp possible here
// Because of the insertions at the end of the loop, we get segfaults :(
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
	for (int i=0; i<ng; i++)
	{
		for (int j=0; j<3; j++)	gridpoint[j]=grid[j][i];
		xacc[i]=gridpoint[0]*cos(b)*cos(l)-gridpoint[1]*cos(b)*sin(l)-gridpoint[2]*sin(b);
		yacc[i]=gridpoint[0]*sin(l)+gridpoint[1]*cos(l);
		zacc[i]=gridpoint[0]*sin(b)*cos(l)-gridpoint[1]*sin(b)*sin(l)+gridpoint[2]*cos(b);
		temporarygridpoint=Point(grid[0][i],grid[1][i],grid[2][i]);
		// also create the map function_values here
		vector<double> velvec = {vx[i], vy[i], vz[i]};
		losvelval = inner_product(unit.begin(),unit.end(),velvec.begin(),0.0);
		losvelmap[temporarygridpoint]=Coord_type(losvelval);
		peakmap[temporarygridpoint]=Coord_type(peakvec[i]);
		fwhmmap[temporarygridpoint]=Coord_type(fwhmvec[i]);
/*		peakmap.insert(make_pair(temporarygridpoint,Coord_type(peakvec[i])));
		fwhmmap.insert(make_pair(temporarygridpoint,Coord_type(fwhmvec[i])));
		losvelmap.insert(make_pair(temporarygridpoint,Coord_type(losvelval)));*/
	}
	double minz=*(min_element(zacc.begin(),zacc.end()));
	double maxz=*(max_element(zacc.begin(),zacc.end()));
	double minx=*(min_element(xacc.begin(),xacc.end()));
	double maxx=*(max_element(xacc.begin(),xacc.end()));
	double miny=*(min_element(yacc.begin(),yacc.end()));
	double maxy=*(max_element(yacc.begin(),yacc.end()));
/*	Value_access peak=Value_access(peakmap);
	Value_access fwhm=Value_access(fwhmmap);
	Value_access losvel=Value_access(losvelmap);*/
	xacc.clear();
	yacc.clear();
	zacc.clear();
	if (commrank==0) cout << "Done!" << endl;

	double lambda0=readgoftfromchianti(chiantifile);
	
	if (commrank==0) cout << "Building frame: " << flush;
	double x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens;
	int li,lj,ind;
	Point p,nearest;
	Delaunay_triangulation_3::Vertex_handle v;
	Delaunay_triangulation_3::Locate_type lt;
	boost::progress_display show_progress((y2-y1+1)*(x2-x1+1));
#ifdef _OPENMP
// it seems as if the first two options below are equivalent and the fastest for this loop
#pragma omp parallel for schedule(dynamic) collapse(2) private (x,y,z,p,lt,li,lj,v,nearest,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens,ind)
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
//#pragma omp parallel for schedule(static,1) collapse(2)
//#pragma omp parallel for schedule(dynamic,(y2-y1+1)*(x2-x1+1)/100) collapse(2)
//#pragma omp parallel for schedule(guided,1) collapse(2)
#endif
	for (int i=y1; i<y2+1; i++)
		for (int j=x1; j<x2+1; j++)
			for (int k=0; k<z_pixel; k++) // scanning through ccd
	{
		x = double(j)/x_pixel*(maxx-minx)+minx;
		y = double(i)/y_pixel*(maxy-miny)+miny;
		z = double(k)/z_pixel*(maxz-minz)+minz;
		// calculate the interpolation in the original frame of reference
		// i.e. derotate the point using angles -l and -b
		p={x*cos(b)*cos(l)+y*sin(l)+z*sin(b)*cos(l),-x*cos(b)*sin(l)+y*cos(l)-z*sin(b)*sin(l),-x*sin(b)+z*cos(b)};
		DTpointer->locate(p, lt, li, lj);

		// Only look for the nearest point and interpolate, if the point p is inside the convex hull.
		if (lt!=Delaunay_triangulation_3::OUTSIDE_CONVEX_HULL)
		{
			v=DTpointer->nearest_vertex(p);
			nearest=v->point();
/* This is how it is done in the CGAL examples		
 			pair<Coord_type,bool> tmppeak=peak(nearest);
			pair<Coord_type,bool> tmpfwhm=fwhm(nearest);
			pair<Coord_type,bool> tmplosvel=losvel(nearest);
			intpolpeak=tmppeak.first;
			intpolfwhm=tmpfwhm.first;
			intpollosvel=tmplosvel.first;*/
			intpolpeak=peakmap[nearest];
			intpolfwhm=fwhmmap[nearest];
			intpollosvel=losvelmap[nearest];
			for (int l=0; l<lambda_pixel; l++)
			{
			// lambda is made around lambda0, with a width of lambda_width 
				lambdaval=double(l)/(lambda_pixel-1)*lambda_width-lambda_width/2.;
				tempintens=intpolpeak*exp(-pow(lambdaval-intpollosvel/speedoflight*lambda0,2)/pow(intpolfwhm,2)*4.*log(2.));
				ind=((i-y1)*(x2-x1+1)+j-x1)*lambda_pixel+l;
				
				results[ind]+=tempintens;
			}
		}
		

		// print progress
		if (k == 0) ++show_progress;
	}
	if (commrank==0) cout << " Done! " << endl << flush;
}

void fillccd(cube observ, const double *results, const int x1, const int x2, const int y1, const int y2)
{
	tgrid grid=observ.readgrid();
	tphysvar intens=observ.readvar(0);
	if (intens.size() == 0) intens.resize(x_pixel*y_pixel*lambda_pixel);
	double lambda0=readgoftfromchianti(chiantifile);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=y1; i<y2+1; i++)
		for (int j=x1; j<x2+1; j++)
			for (int l=0; l<lambda_pixel; l++)
			{
				int ind=((i-y1)*(x2-x1+1)+j-x1)*lambda_pixel+l;
				double lambdaval=double(l)/(lambda_pixel-1)*lambda_width-lambda_width/2.+lambda0;
				grid[0][ind]=j;
				grid[1][ind]=i;
				grid[2][ind]=lambdaval;
				intens[ind]=results[ind];
			}

	observ.setgrid(grid);
	observ.setvar(0,intens);
	observ.settype(observcube);
}

