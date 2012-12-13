#include "../config.h"
#include <iostream>
#include <complex>
#include <vector>

#ifdef HAVEMPI
#include <mpi.h>
#endif

using namespace std;

typedef vector<double> tcoord;
typedef tcoord * tgrid;
typedef vector<double> tphysvar;
typedef tphysvar * tvars;

enum EqType 
{
	builtineq,
	gofttable,
	emisscube,
	empty,
};

// cube class
class cube 
{
	protected:
                EqType qtype;
                int dim;
		int ng;
		int nvars;
		tgrid grid;
		tvars vars;
	public:
		cube(const int invars, const int ingrid, const int = 3);
		~cube();
		int readngrid() const;
		int readnvars() const;
		void settype(EqType intype);
		EqType readtype() const;
		void setgrid(tgrid ingrid);
		tgrid readgrid() const;
		void setvar(const int, const tphysvar);
		tphysvar readvar(const int) const;
		void fillcube();
};

//external functions in readframe.cpp
extern tphysvar log10(tphysvar const&);
extern tphysvar operator*(double const &, tphysvar const &);
extern tphysvar sqrt(tphysvar const&);

// external functions in getarg.cpp
extern void getarg(int, char* array[]);

// external functions in io.cpp
extern int nframes;
extern void getfile();
extern void writefile();
extern void writeparameters(ostream&,char);
extern const int filterlength;
extern void readfilters(const int, double filter[][2]);

//external functions in equilibrium.cpp
extern double density(const double, const double, const double);
extern double temperature(const double, const double, const double);

//external functions in writepng.cpp
extern int writepng(double * const * const image, const int);

//external functions in writemovie.cpp
extern void writemovie(char*, const double, const int);

//external functions in writearray.cpp
extern void writearray(double * const * const, const int);
extern void writearray(tphysvar const &, string);

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
extern void mpi_calculatemypart(double*, const int, const int, const int, const int, const double, cube);

//external functions in writetime.cpp
extern void writetime(double * const * const, const int);

//external functions in progressbar.cpp
extern void progressbar(const int, const int, const int);

// Global variables from mainprog.cpp
// reuse controls whether to compute a new problem or just to reuse the old outputfiles
extern int reuse, png, mpeg, array;
// physical parameters of the problem
extern double length, width, magfield, rhoint, contrast, thickness, alpha;
extern double l, b;

// external functions in equilibrium.cpp
extern const int eqx, eqy, eqz;
extern cube equilibrium();
extern void builtingrid(const int, const int, const int, tgrid);

// external functions from emissioncube.cpp
extern cube emissionfromdatacube(cube);
