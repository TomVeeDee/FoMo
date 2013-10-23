#include "../config.h"
#include <iostream>
#include <complex>
#include <vector>

#ifdef HAVEMPI
#include <mpi.h>
#endif

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>


using namespace std;

typedef vector<double> tcoord;
typedef tcoord * tgrid;
typedef vector<double> tphysvar;
typedef tphysvar * tvars;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_triangulation_3;

enum EqType 
{
	builtineq,
	gofttable,
	emisscube,
	observcube,
	patricksausage,
	stiefkink,
	empty, // this should always be at the end, because it is tested for in the reademissioncube routine in io.cpp
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
		int readdim() const;
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
extern tphysvar pow(const double, tphysvar const&);
extern tphysvar operator/(tphysvar const&, tphysvar const&);
extern tphysvar operator*(tphysvar const&, tphysvar const&);
extern tphysvar log10(tphysvar const&);
extern tphysvar operator*(double const &, tphysvar const &);
extern tphysvar sqrt(tphysvar const&);

// external functions in getarg.cpp
extern void getarg(int, char* array[]);

// external functions in io.cpp
extern void getfile();
extern void writefile();
extern void writeparameters(ostream&, char);
extern void reademissioncube(cube&, const string, Delaunay_triangulation_3* = NULL);
extern void writeemissioncube(const cube, const string, const Delaunay_triangulation_3 * = NULL);

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
extern const int lambda_pixel;
extern double lambda_width;
extern double findmax(double * const * const, int *, int *);
extern double findmin(double * const * const, int *, int *);
extern int MPE_Decomp1d(int , int , int, int *, int*);
extern void fillccd(cube, const double*, const int, const int, const int, const int);
extern void mpi_getcoords(int &, int &, int &, int &);
extern void mpi_calculatemypart(double*, const int, const int, const int, const int, const double, cube, Delaunay_triangulation_3*);

//external functions in writetime.cpp
extern void writetime(double * const * const, const int);

// Global variables from mainprog.cpp
// reuse controls whether to compute a new problem or just to reuse the old outputfiles
extern int reuse, png, mpeg, warray;
// physical parameters of the problem
extern double length, width, magfield, rhoint, contrast, thickness, alpha;
extern double l, b;

// external functions in equilibrium.cpp
extern const int eqx, eqy, eqz;
extern cube equilibrium();
extern void builtingrid(const int, const int, const int, tgrid);

// external functions from emissioncube.cpp
extern cube emissionfromdatacube(cube);
extern Delaunay_triangulation_3 triangulationfromdatacube(cube);
extern string chiantifile;
extern string abundfile;
extern string emissionsave;
extern double readgoftfromchianti(const string);
