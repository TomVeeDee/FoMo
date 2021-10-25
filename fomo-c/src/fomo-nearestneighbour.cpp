#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <set>


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <functional>
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, unsigned> value;

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
const double pi=M_PI; //pi

FoMo::RenderCube nearestneighbourinterpolation(FoMo::GoftCube goftcube, const double l, const double b, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width)
{
//
// results is an array of at least dimension (x2-x1+1)*(y2-y1+1)*lambda_pixel and must be initialized to zero
//
// determine contributions per pixel

	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	//FoMo::DataCube goftcube=object.datacube;
	FoMo::tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();
	int dim=goftcube.readdim();

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) std::cout << "Rotating coordinates to POS reference... " << std::flush;
	std::vector<double> xacc, yacc, zacc, losvel;
	xacc.resize(ng);
	yacc.resize(ng);
	zacc.resize(ng);
	losvel.resize(ng);
	// Read the physical variables
	FoMo::tphysvar peakvec=goftcube.readvar(0);//Peak intensity
	FoMo::tphysvar fwhmvec=goftcube.readvar(1);// line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);

	// initialisations for boost nearest neighbour
	point boostpoint, targetpoint;
	value boostpair;
	std::vector<value> input_values(ng),returned_values;
	box maxdistancebox;

// No openmp possible here
// Because of the insertions at the end of the loop, we get segfaults :(
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private (boostpoint, boostpair)
#endif
	for (int i=0; i<ng; i++)
	{
		std::vector<double> gridpoint;
		gridpoint.resize(dim);
		for (int j=0; j<dim; j++)	gridpoint.at(j)=grid[j][i];
		// if dim==2, then set all z-coordinates to 0.
		if (dim==2) gridpoint.push_back(0.);
		xacc.at(i)=gridpoint.at(0)*cos(b)*cos(l)-gridpoint.at(1)*cos(b)*sin(l)-gridpoint.at(2)*sin(b);// rotated grid
		yacc.at(i)=gridpoint.at(0)*sin(l)+gridpoint.at(1)*cos(l);
		zacc.at(i)=gridpoint.at(0)*sin(b)*cos(l)-gridpoint.at(1)*sin(b)*sin(l)+gridpoint.at(2)*cos(b);
		// Define the unit vector along the line-of-sight
		std::vector<double> unit = {sin(b)*cos(l), -sin(b)*sin(l), cos(b)};
		std::vector<double> velvec = {vx[i], vy[i], vz[i]};// velocity vector
		double losvelval = inner_product(unit.begin(),unit.end(),velvec.begin(),0.0);//velocity along line of sight for position [i]/[ng]
		losvel.at(i)=losvelval;

		// build r-tree from gridpoints, this part is not parallel
		boostpoint = point(gridpoint.at(0), gridpoint.at(1), gridpoint.at(2));
		boostpair=std::make_pair(boostpoint,i);
		input_values.at(i)=boostpair;
	}
	if (commrank==0) std::cout << "Done!" << std::endl;
	if (commrank==0) std::cout << "Building R-tree..." << std::flush;
	// take an rtree with the quadratic packing algorithm, it takes (slightly) more time to build, but queries are faster for large renderings
	bgi::rtree< value, bgi::quadratic<16> > rtree(input_values.begin(),input_values.end());

	// compute the bounds of the input data points, so that we can equidistantly distribute the target pixels
	double minz=*(min_element(zacc.begin(),zacc.end()));
	double maxz=*(max_element(zacc.begin(),zacc.end()));
	double minx=*(min_element(xacc.begin(),xacc.end()));
	double maxx=*(max_element(xacc.begin(),xacc.end()));
	double miny=*(min_element(yacc.begin(),yacc.end()));
	double maxy=*(max_element(yacc.begin(),yacc.end()));
	xacc.clear(); // release the memory
	yacc.clear();
	zacc.clear();
	if (commrank==0) std::cout << "Done!" << std::endl << std::flush;

	std::string chiantifile=goftcube.readchiantifile();
	double lambda0=goftcube.readlambda0();// lambda0=AIA bandpass for AIA imaging
	double lambda_width_in_A=lambda_width*lambda0/speedoflight;

	if (commrank==0) std::cout << "Building frame: " << std::flush;
	double x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens;
	int ind;

	//initialize grids
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel),yvec(x_pixel*y_pixel*lambda_pixel),lambdavec(x_pixel*y_pixel*lambda_pixel);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	if (lambda_pixel > 1) newgrid.push_back(lambdavec);
	FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel,0);

	// maxdistance is the furthest distance between a grid point and a simulation point at which the emission is interpolated
	// it is computed as the half diagonal of the rectangle around this ray, with the sides equal to the x and y distance between rays
	// i.e. it needs to be closer to this ray than to any other ray
	double maxdistance;
	maxdistance = std::sqrt(std::pow((maxx-minx)/(x_pixel-1),2)+std::pow((maxy-miny)/(y_pixel-1),2))/2.;
	// However, it is better to just take the minimum of the pixel size in either direction, because it is then used in the maxdistancebox
	maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/2.;
	// If the viewing is along one of the axis, the previous value does not work very well, and the rendering almost always shows dark stripes: make the value 6 times larger!
	maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/.3;
	if ((maxx-minx)/std::pow(ng,1./3.)>maxdistance || (maxy-miny)/std::pow(ng,1./3.)>maxdistance) std::cout << std::endl << "Warning: maximum distance to interpolated point set to " << maxdistance << "Mm. If it is too small, you have too many interpolating rays and you will have dark stripes in the image plane. Reduce x-resolution or y-resolution." << std::endl;

	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	double deltaz=(maxz-minz);
	if (z_pixel != 1) deltaz/=(z_pixel-1);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(2) shared (rtree) private (x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens,ind,returned_values,targetpoint,maxdistancebox)
#endif
	for (int i=0; i<y_pixel; i++)
		for (int j=0; j<x_pixel; j++)
		{
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel-1)*(maxx-minx)+minx;
			y = double(i)/(y_pixel-1)*(maxy-miny)+miny;

			std::vector<double> p;

			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k=0; k<z_pixel; k++) // scanning through ccd
			{
				z = double(k)*deltaz+minz;
		// calculate the interpolation in the original frame of reference
		// i.e. derotate the point using angles -l and -b
				p={x*cos(b)*cos(l)+y*sin(l)+z*sin(b)*cos(l),-x*cos(b)*sin(l)+y*cos(l)-z*sin(b)*sin(l),-x*sin(b)+z*cos(b)};

				// initialise nearestindex to point -1
				int nearestindex=-1;

				// look for nearest point to targetpoint
				targetpoint=point(p.at(0),p.at(1),p.at(2));
				returned_values.clear();
				// the second condition ensures the point is not further away than
				// - half the x-resolution in the x-direction
				// - half the y-resolution in the y-direction
				// - the maximum of both the previous numbers in the z-direction (sort of improvising a convex hull approach)
				maxdistancebox=box(point(p.at(0)-(maxx-minx)/(x_pixel-1)/2.,p.at(1)-(maxy-miny)/(y_pixel-1)/2.,p.at(2)-maxdistance),point(p.at(0)+(maxx-minx)/(x_pixel-1)/2.,p.at(1)+(maxy-miny)/(y_pixel-1)/2.,p.at(2)+maxdistance));
				// it seems the expression above is the culprit for simulations with very stretched grids producing striped emissions, let's make the box of size maxdistance
				maxdistancebox=box(point(p.at(0)-maxdistance,p.at(1)-maxdistance,p.at(2)-maxdistance),point(p.at(0)+maxdistance,p.at(1)+maxdistance,p.at(2)+maxdistance));
				//numberofpoints=
				rtree.query(bgi::nearest(targetpoint, 1) && bgi::within(maxdistancebox), std::back_inserter(returned_values));

				if (returned_values.size() >= 1)
				{
					nearestindex=returned_values.at(0).second;
					intpolpeak=peakvec.at(nearestindex);
					intpolfwhm=fwhmvec.at(nearestindex);
					intpollosvel=losvel.at(nearestindex);
				}
				else
				{
					intpolpeak=0;
				}

				if (lambda_pixel>1)// spectroscopic study
				{
					for (int il=0; il<lambda_pixel; il++) // changed index from global variable l into il [D.Y. 17 Nov 2014]
					{
						// lambda the relative wavelength around lambda0, with a width of lambda_width
						lambdaval=double(il)/(lambda_pixel-1)*lambda_width_in_A-lambda_width_in_A/2.;
						// here a ternary operator may be used
						// if intpolpeak is not zero then the correct expression is used. otherwise, the intensity is just 0
						// it remains to be tested if this is faster than just the direct computation
						tempintens=intpolpeak ? intpolpeak*exp(-pow(lambdaval-intpollosvel/speedoflight*lambda0,2)/pow(intpolfwhm,2)*4.*log(2.)) : 0; // Uncommented this line by Vaibhav pant on 22 Nov, 2018. tempintens as defined below was giving NAN values for odd wavelength bins.

						//tempintens=intpolpeak*exp(-pow(lambdaval-intpollosvel/speedoflight*lambda0,2)/pow(intpolfwhm,2)*4.*log(2.));

						ind=(i*(x_pixel)+j)*lambda_pixel+il;//
						newgrid.at(0).at(ind)=x;
						newgrid.at(1).at(ind)=y;
						newgrid.at(2).at(ind)=lambdaval+lambda0; // store the full wavelength
						// this is critical, but with tasks, the ind is unique for each task, and no collision should occur
						intens.at(ind)+=tempintens;// loop over z and lambda [D.Y 17 Nov 2014]
					}
				}
			
				if (lambda_pixel==1) // AIA imaging study
				{
					tempintens=intpolpeak;
					ind=(i*x_pixel+j);
					newgrid.at(0).at(ind)=x;
					newgrid.at(1).at(ind)=y;
					// this is critical, but with tasks, the ind is unique for each task, and no collision should occur
					intens.at(ind)+=tempintens;// loop over z and lambda [D.Y 17 Nov 2014]
				}

			// print progress
				++show_progress;
			}
		}
	if (commrank==0) std::cout << " Done! " << std::endl << std::flush;

	FoMo::RenderCube rendercube(goftcube);
	FoMo::tvars newdata;
	double pathlength=(maxz-minz)/(z_pixel-1);
	// this does not work if only one z_pixel is given (e.g. for a 2D simulation), or the maxz and minz are equal (face-on on 2D simulation)
	// assume that the thickness of the slab is 1Mm.
	if ((maxz==minz) || (z_pixel==1))
	{
		pathlength=1.;
		std::cout << "Assuming that this is a 2D simulation: setting thickness of simulation to " << pathlength << "Mm." << std::endl << std::flush;
	}
	
	// the intensity does not need to be rescaled for spectroscopic data
	double apix = 1.;
	// these are the default units if spectroscopic data
	std::vector<std::string> unitvec;
	unitvec.push_back("Mm");
	unitvec.push_back("Mm");
	if (lambda_pixel > 1) unitvec.push_back("\\AA{}");
	unitvec.push_back("erg cm^{-2} s^{-1} \\AA{}^{-1}");
	// check if the units of emissivity contain DN: then we are dealing with an instrument, and spatial units should be converted to arcsec, 
	// intensity should be rescaled to pixel size
	std::size_t found = goftcube.readunit().at(dim).find("DN");
	if (found!=std::string::npos)
	{
		unitvec.at(0)="arcsec";
		newgrid.at(0)=FoMo::operator*(1./Mmperarcsec,newgrid.at(0));
		unitvec.at(1)="arcsec";
		newgrid.at(1)=FoMo::operator*(1./Mmperarcsec,newgrid.at(1));
		float dx=(maxx-minx)/(x_pixel-1),dy=(maxy-miny)/(y_pixel-1); // are given in Mm
		apix = (dx/Mmperarcsec)*(dy/Mmperarcsec)*pow(pi/180./3600.,2); 
		std::cout << "apix" << apix << "dx" << dx << "dy" << dy;
		unitvec.back()="DN s^{-1} pixel^{-1}"; // this could be improved using Boost::units, making everything automatic, including compiler checks
	}
	
	intens=FoMo::operator*(pathlength*1e8*apix,intens); // assume that the coordinates in goftcube are given in Mm, and convert to cm
	newdata.push_back(intens);
	rendercube.setdata(newgrid,newdata,&unitvec);
	rendercube.setrendermethod("NearestNeighbour");
	rendercube.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	if (lambda_pixel == 1)
	{
		rendercube.setobservationtype(FoMo::Imaging);
	}
	else
	{
		rendercube.setobservationtype(FoMo::Spectroscopic);
	}
	return rendercube;
}

namespace FoMo
{
	FoMo::RenderCube RenderWithNearestNeighbour(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		FoMo::RenderCube rendercube(goftcube);
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit)
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit)
			{
				rendercube=nearestneighbourinterpolation(goftcube,*lit,*bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
				rendercube.setangles(*lit,*bit);
				std::stringstream ss;
				// if outfile is "", then this should not be executed.
				ss << outfile;
				ss << "l";
				ss << std::setfill('0') << std::setw(3) << std::round(*lit/pi*180.);
				ss << "b";
				ss << std::setfill('0') << std::setw(3) << std::round(*bit/pi*180.);
				ss << ".txt";
				rendercube.writegoftcube(ss.str());
				ss.str("");
			}
		return rendercube;
	}
}
