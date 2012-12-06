#include "../config.h"
#include <iostream>
#include <complex>

#ifdef HAVEMPI
#include <mpi.h>
#endif

using namespace std;

// fortran-subroutines in eigenfunction.f90
extern "C" {
 	void readng_(int*);
	void eigvfk_(const int*, const int*, const int*, double func[][4]);
}

// external functions in getarg.cpp
extern void getarg(int, char* array[]);

// external functions in io.cpp
extern void getfile();
extern void writefile();
extern void writeparameters(ostream&,char);
extern const int filterlength;
extern void readfilters(const int, double filter[][2]);

//external functions in calculatepollux.cpp
extern void calculatepollux();
extern void readfrequency();

//external functions in PARALLEL.f90
extern "C" {
	void calcpollux_(const double*, const double*, const double*, const double*, const double*, const double*, double*, double*);
}

//external functions in equilibrium.cpp
extern double density(const double, const double, const double);
extern double temperature(const double, const double, const double);

//external functions in writepng.cpp
extern int writepng(double * const * const image, const double);

//external functions in writemovie.cpp
extern void writemovie(char*, const double, const double);

//external functions in writearray.cpp
extern void writearray(double * const * const, const double);

//external functions in fillccd.cpp
extern const int x_pixel;
extern const int y_pixel;
// z_pixel is not really a constant of the telescope, it can be adjusted to obtain a sufficient resolution
extern const int z_pixel;
extern const double size_pixel;
extern const double size_z_pixel;
extern double findmax(double * const * const, int *, int *);
extern double findmin(double * const * const, int *, int *);
extern void fillccd(double * const * const image, const double*, const int, const int, const int, const int);
extern void mpi_getcoords(int &, int &, int &, int &);
extern void mpi_calculatemypart(double*, const int, const int, const int, const int, const double, double **);

//external functions in perturbation.cpp
extern const int ncols; // ncols determines the number of columns in the perturbed loop.
// for now the columns contain:
// 1: x coordinate
// 2: y coordinate --> in the cylindrical coordinate system aligned with the coronal loop
// 3: z coordinate
// 4: perturbed density
// 5: perturbed temperature
extern const int npointz, npointphi;
extern int ngridpoints();
extern void physicalquants(double **, const double);
extern double perturbeddensity(const double, const double, const double, const double, double **, const int);

//external functions in writetime.cpp
extern void writetime(double * const * const, const double);

//external functions in progressbar.cpp
extern void progressbar(const int, const int, const int);

// Global variables from mainprog.cpp
// reuse controls whether to compute a new problem or just to reuse the old outputfiles
extern int reuse, png, mpeg, array;
// physical parameters of the problem
extern double length, width, beta, magfield, rhoint, contrast, thickness, alpha, ampl, phase;
extern double psi, l, b;
extern double nperiods, nframes;
extern complex<double> frequency;

