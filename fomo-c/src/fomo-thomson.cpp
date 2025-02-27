#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <sstream>
#include <iomanip>
#include <cstdlib> 
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>
#include <cmath> 
#include <vector>
#include <numeric> 
#include <algorithm> 
#include <cassert> 
#include <set>
#include <iostream> 

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
const double RSUN=695.7; // radius of the Sun in Mm
const double U=0.58; // limb darkening coefficient
const double re=2.8179e-13; // electron radii in cm
const double L=2.8e10; // the irradiance in vertical direction at 1.5Rs, if irradiance = luminosity/ 4pi distance^2 in ergs/s/cm-2

//! Calculate and return the Thomson scattering geometric functions
// from raytrace.h, but with rho dependence removed (rho/r=sin(\bar{chi}) from Inhester 2016, after Eq. 3.38)

inline void getThomsonGeomFactor(const float &r,float &totterm,float &polterm) {
    float sinomega=RSUN/r;
    float Coeff=0.5*pi*std::pow(re,2)*L; //coefficient pi*re^2*L/2 if L=Laverage in ergs/s/cm-2
    float sinsquareomega=sinomega*sinomega;
    float cossquareomega=1-sinsquareomega;
    float cosomega=sqrt(cossquareomega);   
    float logterm=log((1.+sinomega)/cosomega);
    float a, b, c, d;
    
        a = cosomega * sinsquareomega;
        b = (-1.0 / 8.0) * (1.0 - 3.0 * sinsquareomega - cossquareomega * ((1.0 + 3.0 * sinsquareomega) / sinomega) * logterm);
        c = (4.0 / 3.0) - cosomega - (cosomega * cossquareomega) / 3.0;
        d = (1.0 / 8.0) * (5.0 + sinsquareomega - cossquareomega * ((5.0 - sinsquareomega) / sinomega) * logterm);

    // Polarised components of the intensity 
    polterm=Coeff*(1./(1.-U/3.))*(a+U*(b-a)); //*sin(\bar{chi})
    // Total components of the intensity 
    totterm=Coeff*(1./(1.-U/3.))*(2.*(c+U*(d-c))); //-polterm
//    	std::cout << "polterm" << polterm << "totterm" << totterm << "coeff"<< Coeff << std::endl;
}

FoMo::GoftCube FoMo::thomsonfromdatacube(FoMo::DataCube datacube)
{
	// create goft, and manually convert with emission
	// emission is computed with raytrace.h (using losinteg and getthompsongeomfactor, probably copy paste the latter one)
	FoMo::GoftCube emission(datacube);
	FoMo::tphysvar densvar=datacube.readvar(0);
	FoMo::tphysvar totint, polint;
	FoMo::tvars exporteddata;
	// also compute r from grid positions
	float r,totterm,polterm;
	int npoints=datacube.readngrid();
	FoMo::tgrid grid=datacube.readgrid();
	// this loop can be parallelised easily with openmp
	for (int i =0; i<npoints; i++)
	{
		r=std::sqrt(std::pow(grid[0][i],2)+std::pow(grid[1][i],2)+std::pow(grid[2][i],2));
		getThomsonGeomFactor(r,totterm,polterm);
		
		// set both totterm and polterm in goftcube, and take care of integration in object.render(), 
		// intensity will be I=totterm - polterm*sin(angle radial and los)**2 (inhester 2016)
		totint.push_back(densvar.at(i)*totterm);
		polint.push_back(densvar.at(i)*polterm);
	}
	exporteddata.push_back(totint);
	exporteddata.push_back(polint);
	emission.setdata(grid,exporteddata);
	std::string emissname="thomson scattering";
	emission.setchiantifile(emissname);
	emission.setabundfile(emissname);
	emission.setlambda0(0.);
	
	return emission;
}

FoMo::RenderCube thomsoninterpolation(FoMo::GoftCube goftcube, const double l, const double b, const int x_pixel, const int y_pixel, const int z_pixel)
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
	FoMo::tphysvar totvec=goftcube.readvar(0);//unpolarised intensity 
	FoMo::tphysvar polvec=goftcube.readvar(1);//polarised intensity
	
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
		// Compute trigonometric values
    	const double epsilon = 1e-6; 
    	double cos_b = cos(b);
    	double sin_b = sin(b);
    	double cos_l = cos(l);
    	double sin_l = sin(l);
		
		// Avoid division by zero or too small numbers
		if (std::abs(cos_b) < epsilon) {
    		cos_b = epsilon;
		}
		
		if (std::abs(sin_b) < epsilon) {
        	sin_b = epsilon;  
    	}
    	
		std::vector<double> gridpoint;
		gridpoint.resize(dim);
		for (int j=0; j<dim; j++)	gridpoint.at(j)=grid[j][i];
		// if dim==2, then set all z-coordinates to 0. 
		if (dim==2) gridpoint.push_back(0.);
		
		//Default calculation of the grid rotation 
		//xacc.at(i)=gridpoint.at(0)*cos_b*cos(l)-gridpoint.at(1)*cos_b*sin(l)-gridpoint.at(2)*sin_b;// rotated grid
		//yacc.at(i)=gridpoint.at(0)*sin(l)+gridpoint.at(1)*cos(l);
		//zacc.at(i)=gridpoint.at(0)*sin_b*cos(l)-gridpoint.at(1)*sin_b*sin(l)+gridpoint.at(2)*cos_b;
		 
		// WARNING: A special treatment is applied to avoid numerical instabilities.
		// We introduce a smooth transition between the full spherical 
		// rotation and a simplified pole rotation (pure Z-axis rotation). 
		// The transition is controlled by a sigmoid blending function with configurable 
		// threshold and steepness parameters.
		// Users should be aware that this treatment ensures numerical stability 
		// but introduces a small approximation in the transition region.

		double threshold = 1.309;  // In radians, equivalent to ~85 degrees
    	double steepness = 20.0;  // Controls the sharpness of the transition
		
		double blend = 1.0 / (1.0 + exp(-steepness * (std::fabs(b) - threshold)));
		
		// Calculate the full rotation using spherical coordinates
    	double full_x = gridpoint.at(0) * cos_b * cos_l - gridpoint.at(1) * cos_b * sin_l - gridpoint.at(2) * sin_b;
    	double full_y = gridpoint.at(0) * sin_l + gridpoint.at(1) * cos_l;
    	double full_z = gridpoint.at(0) * sin_b * cos_l - gridpoint.at(1) * sin_b * sin_l + gridpoint.at(2) * cos_b;

    	// Rotation at the poles (purely around the Z-axis)
    	double pole_x = gridpoint.at(0) * cos_l - gridpoint.at(1) * sin_l;
    	double pole_y = gridpoint.at(0) * sin_l + gridpoint.at(1) * cos_l;
   		double pole_z = gridpoint.at(2);  // No change in Z at the pole

    	// Smoothly blend between the full rotation and the pole rotation
    	xacc.at(i) = full_x * (1.0 - blend) + pole_x * blend;
    	yacc.at(i) = full_y * (1.0 - blend) + pole_y * blend;
    	zacc.at(i) = full_z * (1.0 - blend) + pole_z * blend;
		
		// Define the unit vector along the line-of-sight
		std::vector<double> unit = {sin_b*cos(l), -sin_b*sin(l), cos_b};
		std::vector<double> posvec = {xacc.at(i), yacc.at(i), zacc.at(i)};// position vector
		// need to compute sin^2 \chi (angle between LOS and radial direction)
		// sin^2 chi=1-cos^2 chi = 1- (r\cdot unit)^2/r^2
		double losvelval = 1.- std::pow(inner_product(unit.begin(),unit.end(),posvec.begin(),0.0),2)/inner_product(posvec.begin(),posvec.end(),posvec.begin(),0.0);//angle with line of sight for position [i]/[ng]
		losvel.at(i)=losvelval;
		
		// build r-tree from gridpoints, this part is not parallel
		boostpoint = point(gridpoint.at(0), gridpoint.at(1), gridpoint.at(2));
		boostpair=std::make_pair(boostpoint,i);
		input_values.at(i)=boostpair;
	}
	if (commrank==0) std::cout << "Done!" << std::endl;
	if (commrank==0) std::cout << "Building R-tree... " << std::flush;
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

	if (commrank==0) std::cout << "Building frame: " << std::flush;
	double x,y,z,intpoltot,intpolpol,intpollos,tempintens1,tempintens2;
	int ind;
	
	//initialize grids
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel),yvec(x_pixel*y_pixel);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	FoMo::tphysvar intens(x_pixel*y_pixel,0);
	
	// maxdistance is the furthest distance between a grid point and a simulation point at which the emission is interpolated
	// it is computed as the half diagonal of the rectangle around this ray, with the sides equal to the x and y distance between rays
	// i.e. it needs to be closer to this ray than to any other ray
	double maxdistance; 
	//maxdistance = std::sqrt(std::pow((maxx-minx)/(x_pixel-1),2)+std::pow((maxy-miny)/(y_pixel-1),2))/2.;
	// However, it is better to just take the minimum of the pixel size in either direction, because it is then used in the maxdistancebox
	//maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/2.;
	// If the viewing is along one of the axis, the previous value does not work very well, and the rendering almost always shows dark stripes: make the value 6 times larger!
	// Default approach to calculate maxdistance
	//maxdistance = std::max((maxx-minx)/(x_pixel-1),(maxy-miny)/(y_pixel-1))/0.5; 
	//if ((maxx-minx)/std::pow(ng,1./3.)>maxdistance || (maxy-miny)/std::pow(ng,1./3.)>maxdistance) std::cout << std::endl << "Warning: maximum distance to interpolated point set to " << maxdistance << "Mm. If it is too small, you have too many interpolating rays and you will have dark stripes in the image plane. Reduce x-resolution or y-resolution." << std::endl;
	
	
	// WARNING: A special treatment is applied to avoid numerical instabilities.	
	// Calculate the latitude factor using a smooth transition
	double blending_factor = 1.0 / (1.0 + exp(-10 * (std::fabs(b) - 0.1745))); // Sigmoid-like transition around b = 10°

	// Near equator formula
	double equator_formula = std::max(0.5 * std::sqrt(std::pow((maxx - minx) / (x_pixel - 1), 2) +
                                                 std::pow((maxy - miny) / (y_pixel - 1), 2)),
                                  0.1 * (maxx - minx) / (x_pixel - 1));

	// Away from equator formula
	double away_from_equator_formula = std::max((maxx - minx) / (x_pixel - 1), (maxy - miny) / (y_pixel - 1)) / 0.005;

	// Smoothly blend the two formulas based on latitude
	maxdistance = equator_formula * (1.0 - blending_factor) + away_from_equator_formula * blending_factor;
	
	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	double deltaz=(maxz-minz);
	std::cout << "deltaz" << deltaz << std::endl;
	if (z_pixel != 1) deltaz/=(z_pixel-1);
	std::cout << "\nWarning: Numerical treatment applied near poles to improve stability.\n";
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(2) shared (rtree) private (x,y,z,intpoltot,intpolpol,intpollos,tempintens1,tempintens2,ind,returned_values,targetpoint,maxdistancebox)
#endif
	for (int i=0; i<y_pixel; i++)
		for (int j=0; j<x_pixel; j++)
		{
			// Compute trigonometric values
    			double cos_b = cos(b);
    			double sin_b = sin(b);
    			
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel-1)*(maxx-minx)+minx;
			y = double(i)/(y_pixel-1)*(maxy-miny)+miny;
						
			std::vector<double> p;
			// Initialize the tempintensMatrix to store intensity values
    		//std::vector<std::vector<double>> tempintensMatrix(y_pixel, std::vector<double>(2, 0.0)); // Two columns for two intensity types

			
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k=0; k<z_pixel; k++) // scanning through ccd
			{
				z = double(k)*deltaz+minz;
		// calculate the interpolation in the original frame of reference
		// i.e. derotate the point using angles -l and -b
				p={x*cos_b*cos(l)+y*sin(l)+z*sin_b*cos(l),-x*cos_b*sin(l)+y*cos(l)-z*sin_b*sin(l),-x*sin_b+z*cos_b};
				
				// initialise nearestindex to point -1
				int nearestindex=-1;
				
				// look for nearest point to targetpoint
				targetpoint=point(p.at(0),p.at(1),p.at(2));				returned_values.clear();
				// the second condition ensures the point is not further away than 
				// - half the x-resolution in the x-direction
				// - half the y-resolution in the y-direction
				// - the maximum of both the previous numbers in the z-direction (sort of improvising a convex hull approach)
				maxdistancebox=box(point(p.at(0)-(maxx-minx)/(x_pixel-1)/2.,p.at(1)-(maxy-miny)/(y_pixel-1)/2.,p.at(2)-maxdistance),point(p.at(0)+(maxx-minx)/(x_pixel-1)/2.,p.at(1)+(maxy-miny)/(y_pixel-1)/2.,p.at(2)+maxdistance));
				// it seems the expression above is the culprit for simulations with very stretched grids producing striped emissions, let's make the box of size maxdistance
				maxdistancebox=box(point(p.at(0)-maxdistance,p.at(1)-maxdistance,p.at(2)-maxdistance),point(p.at(0)+maxdistance,p.at(1)+maxdistance,p.at(2)+maxdistance));
				//numberofpoints=
				rtree.query(bgi::nearest(targetpoint, 5) && bgi::within(maxdistancebox), std::back_inserter(returned_values));

				if (returned_values.size() >= 1)
				{
					nearestindex=returned_values.at(0).second;
					intpoltot=totvec.at(nearestindex);
					intpolpol=polvec.at(nearestindex);
					intpollos=losvel.at(nearestindex);
				}
				else
				{
					intpoltot=0;
					intpolpol=0;
				}
			 // Calculate
            	tempintens1 = 2 * intpoltot - intpolpol * intpollos;  // Total brightness
                tempintens2 = intpolpol * intpollos;                  // Polarised brightness
				
				if (tempintens1 < 0.) tempintens1=0;
				if (tempintens2 < 0.) tempintens2=0;

            	int ind = (i * x_pixel + j);
            	
            	newgrid.at(0).at(ind) = x;
            	newgrid.at(1).at(ind) = y;
			// WARNING: Choose between Total and Polarized Brightness
            // Update intensity array on each iteration if needed
            //	intens.at(ind) +=tempintens1 ; // Compute Total Brightness //loop over z and lambda [D.Y 17 Nov 2014]
            	intens.at(ind) +=tempintens2 ; // Compute Polarized Brightness
            
			// print progress
				++show_progress;
			}
		}
	if (commrank==0) std::cout << " Done! " << std::endl << std::flush;
	// WARNING: Choose between Total and Polarized Brightness
	std::cout << "we calculate pB brightness" << std::endl << std::flush;
    //std::cout << "we calculate total brightness" << std::endl << std::flush;
    FoMo::RenderCube rendercube(goftcube);
	FoMo::tvars newdata;
	double pathlength=(maxz-minz)/(z_pixel-1);
	// this does not work if only one z_pixel is given (e.g. for a 2D simulation), or the maxz and minz are equal (face-on on 2D simulation)
	// assume that the thickness of the slab is 1Mm. 
	if ((maxz==minz) || (z_pixel==1))
	{
		pathlength=1.;
		std::cout << "Assuming that this is a 2D simulations: setting thickness of simulation to " << pathlength << "Mm." << std::endl << std::flush;
	}
	//std::cout << "pathlength" << pathlength << "Mm." << std::endl << std::flush;
	intens=FoMo::operator*(pathlength*1e8,intens); // assume that the coordinates are given in Mm, and convert to cm
	newdata.push_back(intens);
	rendercube.setdata(newgrid,newdata);
	rendercube.setrendermethod("Thomson");
	rendercube.setresolution(x_pixel,y_pixel,z_pixel,1.,100000.);
	rendercube.setobservationtype(FoMo::Imaging);

	return rendercube;
}

namespace FoMo
{
	FoMo::RenderCube RenderWithThomson(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		FoMo::RenderCube rendercube(goftcube);
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit)
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit)
			{
				rendercube=thomsoninterpolation(goftcube,*lit,*bit, x_pixel, y_pixel, z_pixel);
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

