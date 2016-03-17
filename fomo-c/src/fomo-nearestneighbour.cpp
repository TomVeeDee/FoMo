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

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
const double pi=M_PI; //pi
const unsigned int subsetgoftn=10; // include subsetgoftn * z_pixel points in the subset of points along the ray

template <typename T>
T norm(const std::vector<T> a)
{
	return std::sqrt(inner_product(a.begin(),a.end(),a.begin(),0.));
}

template <typename T> 
T distancebetweenpoints(const std::vector<T> a, const std::vector<T> b)
{
	assert(a.size() == b.size());
	std::vector<T> intermediate(a.size());
	// subtract the vectors from each other, and store in intermediate
	std::transform(a.begin(), a.end(), b.begin(), intermediate.begin(), std::minus<T>());
	// then take the norm of the difference vector
	return norm(intermediate);	
}


template <typename T>
T distancetoline(const std::vector<T> point, const std::vector<T> pointonline, const std::vector<T> directionofline)
{
	assert(point.size() == pointonline.size());
	assert(point.size() == directionofline.size());
	
	std::vector<T> pointdiff(point.size());
	std::transform(point.begin(),point.end(),pointonline.begin(),pointdiff.begin(),std::minus<T>());
	T pointdist = norm(pointdiff);
	T normunit = norm(directionofline);
	T temporary = inner_product(pointdiff.begin(),pointdiff.end(),directionofline.begin(),0.);
	// Eq. 6 on http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	return std::sqrt(std::pow(pointdist,2)*std::pow(normunit,2)-std::pow(temporary,2))/normunit;
}

template <typename T>
T coordinateonline(const std::vector<T> point, const std::vector<T> pointonline, const std::vector<T> directionofline)
{
	assert(point.size() == pointonline.size());
	assert(point.size() == directionofline.size());
	
	std::vector<T> pointdiff(point.size());
	std::transform(point.begin(),point.end(),pointonline.begin(),pointdiff.begin(),std::minus<T>());
	T normunit = norm(directionofline);
	T temporary = inner_product(pointdiff.begin(),pointdiff.end(),directionofline.begin(),0.);
	// Eq. 3 on http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	return -temporary/normunit;
}

struct point {
	double x;
	double y;
	double z;
	int index;
	double distance;
};

struct by_distance {
	bool operator()(point const &a, point const &b)
	{
		return a.distance < b.distance;
	}
	bool operator()(point const &a, double const val)
	{
		return a.distance < val;
	}
};

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
	std::vector<double> gridpoint;
	gridpoint.resize(dim);
	// Define the unit vector along the line-of-sight
	std::vector<double> unit = {sin(b)*cos(l), -sin(b)*sin(l), cos(b)};
	// Read the physical variables
	FoMo::tphysvar peakvec=goftcube.readvar(0);//Peak intensity 
	FoMo::tphysvar fwhmvec=goftcube.readvar(1);// line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);  
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);

// No openmp possible here
// Because of the insertions at the end of the loop, we get segfaults :(
// I think it is possible now, if we make losvelval private.
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
	for (int i=0; i<ng; i++)
	{
		for (int j=0; j<dim; j++)	gridpoint[j]=grid[j][i];
		xacc[i]=gridpoint[0]*cos(b)*cos(l)-gridpoint[1]*cos(b)*sin(l)-gridpoint[2]*sin(b);// rotated grid
		yacc[i]=gridpoint[0]*sin(l)+gridpoint[1]*cos(l);
		zacc[i]=gridpoint[0]*sin(b)*cos(l)-gridpoint[1]*sin(b)*sin(l)+gridpoint[2]*cos(b);
		std::vector<double> velvec = {vx[i], vy[i], vz[i]};// velocity vector
		double losvelval = inner_product(unit.begin(),unit.end(),velvec.begin(),0.0);//velocity along line of sight for position [i]/[ng]
		losvel.at(i)=losvelval;
	}
	double minz=*(min_element(zacc.begin(),zacc.end()));
	double maxz=*(max_element(zacc.begin(),zacc.end()));
	double minx=*(min_element(xacc.begin(),xacc.end()));
	double maxx=*(max_element(xacc.begin(),xacc.end()));
	double miny=*(min_element(yacc.begin(),yacc.end()));
	double maxy=*(max_element(yacc.begin(),yacc.end()));
	xacc.clear(); // release the memory
	yacc.clear();
	zacc.clear();
	if (commrank==0) std::cout << "Done!" << std::endl;

	std::string chiantifile=goftcube.readchiantifile();
	double lambda0=goftcube.readlambda0();// lambda0=AIA bandpass for AIA imaging
       	
	if (commrank==0) std::cout << "Building frame: " << std::flush;
	double x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens;
	int ind;
	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	
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
	double maxdistance = std::sqrt(std::pow((maxx-minx)/(x_pixel-1),2)+std::pow((maxy-miny)/(y_pixel-1),2))/2.;
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(2) private (x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens,ind)
#endif
	for (int i=0; i<y_pixel; i++)
		for (int j=0; j<x_pixel; j++)
		{
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel-1)*(maxx-minx)+minx;
			y = double(i)/(y_pixel-1)*(maxy-miny)+miny;
			z = minz;
			// thus the coordinates of this point in the (x,y,z) of the datacube are
			std::vector<double> startp={x*cos(b)*cos(l)+y*sin(l)+z*sin(b)*cos(l),-x*cos(b)*sin(l)+y*cos(l)-z*sin(b)*sin(l),-x*sin(b)+z*cos(b)};
			// for each ray, we select a subset of points closest to the ray
			std::vector<point> subsetofpoints;
			point tmppoint;
			
			for (int l=0; l<ng; l++)
			{
				//compute distance of each point to ray
				tmppoint.x=grid[0][l];
				tmppoint.y=grid[1][l];
				tmppoint.z=grid[2][l];
				tmppoint.index=l;
				tmppoint.distance=distancetoline({tmppoint.x,tmppoint.y,tmppoint.z},startp,unit);
				
				// if the distance is small enough, we need to add the point
				if (tmppoint.distance < maxdistance)
				{
					// find at which position the value should be inserted
					// this ensures that the vector remains sorted
					auto low = std::lower_bound(subsetofpoints.begin(), subsetofpoints.end(), tmppoint.distance, by_distance());
					subsetofpoints.insert(low,tmppoint);
					// if the subset is already too large, we need to get rid of the last element
					if (subsetofpoints.size()>subsetgoftn*z_pixel) subsetofpoints.pop_back();
				}
			}
			
			double mincoordinate, maxcoordinate, newcoordinate, pointcoordinate;
			// build in a check for the case that there are no close points, otherwise it segfaults
			if (subsetofpoints.size()==0) 
			{
				mincoordinate=0;
			}
			else
			{
				mincoordinate=coordinateonline({subsetofpoints.at(0).x,subsetofpoints.at(0).y,subsetofpoints.at(0).z},startp,unit);
			}
			maxcoordinate=mincoordinate;
			for (unsigned int l=1; l<subsetofpoints.size(); l++)
			{
				// compute the coordinate of each point on the line and store it in the struct (new member, to define)
				// this can be done with formula 3 on http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
				newcoordinate=coordinateonline({subsetofpoints.at(l).x,subsetofpoints.at(l).y,subsetofpoints.at(l).z},startp,unit);
				// go through subsetofpoints and use coordinateonline, keep the min value and max value.
				if (newcoordinate < mincoordinate) mincoordinate = newcoordinate;
				if (newcoordinate > maxcoordinate) maxcoordinate = newcoordinate;
			}
			
			std::vector<double> p;
			
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k=0; k<z_pixel; k++) // scanning through ccd
			{
				z = double(k)/(z_pixel-1)*(maxz-minz)+minz;
		// calculate the interpolation in the original frame of reference
		// i.e. derotate the point using angles -l and -b
				p={x*cos(b)*cos(l)+y*sin(l)+z*sin(b)*cos(l),-x*cos(b)*sin(l)+y*cos(l)-z*sin(b)*sin(l),-x*sin(b)+z*cos(b)};
				pointcoordinate=coordinateonline(p,startp,unit);
				
		// Only look for the nearest point and interpolate, if the point p is inside the convex hull.
		// Check if the z-coordinate is between the minimum and maximum coordinate of subsetofpoints
				if (pointcoordinate >= mincoordinate && pointcoordinate <= maxcoordinate)
				{
					for (int m=0; m<dim; m++) gridpoint.at(m)=grid[m][0];
					double nearestdistance=distancebetweenpoints(p,gridpoint);
					int nearestindex=0;
				
					for (auto subsetpointer=subsetofpoints.begin(); subsetpointer != subsetofpoints.end(); ++subsetpointer)
					{
						gridpoint.at(0)=(*subsetpointer).x;
						gridpoint.at(1)=(*subsetpointer).y;
						gridpoint.at(2)=(*subsetpointer).z;
						double tempdistance=distancebetweenpoints(p,gridpoint);
						if (tempdistance < nearestdistance)
						{
							nearestindex=(*subsetpointer).index;
							nearestdistance=tempdistance;
						}
					}

					if (nearestdistance < maxdistance)
					{
						intpolpeak=peakvec.at(nearestindex);
						intpolfwhm=fwhmvec.at(nearestindex);
						intpollosvel=losvel.at(nearestindex);
					}
					else
					{
						intpolpeak=0;
					}
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
						lambdaval=double(il)/(lambda_pixel-1)*lambda_width-lambda_width/2.;
						// here a ternary operator may be used
						// if intpolpeak is not zero then the correct expression is used. otherwise, the intensity is just 0
						// it remains to be tested if this is faster than just the direct computation
						//tempintens=intpolpeak ? intpolpeak*exp(-pow(lambdaval-intpollosvel/speedoflight*lambda0,2)/pow(intpolfwhm,2)*4.*log(2.)) : 0;
						tempintens=intpolpeak*exp(-pow(lambdaval-intpollosvel/speedoflight*lambda0,2)/pow(intpolfwhm,2)*4.*log(2.));
						ind=(i*(x_pixel)+j)*lambda_pixel+il;// 
						newgrid.at(0).at(ind)=x;
						newgrid.at(1).at(ind)=y;
						newgrid.at(2).at(ind)=lambdaval+lambda0; // store the full wavelength
						// this is critical, but with tasks, the ind is unique for each task, and no collision should occur
						intens.at(ind)+=tempintens;// loop over z and lambda [D.Y 17 Nov 2014]
					}
				}
			
				if (lambda_pixel==1) // AIA imaging study. Algorithm not verified [DY 14 Nov 2014]
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
	intens=FoMo::operator*(pathlength*1e8,intens); // assume that the coordinates are given in Mm, and convert to cm
	newdata.push_back(intens);
	rendercube.setdata(newgrid,newdata);
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
	FoMo::RenderCube RenderWithNearestNeighbour(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		/* A good speedup would be to calculate the triangulation per ray.
		 * It would be good to select only the points around the ray, make the triangulation of that.
		 * The number of points would be drastically reduced, and the triangulation would be greatly sped up.
		 */
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
				ss << std::setfill('0') << std::setw(3) << *lit/pi*180.;
				ss << "b";
				ss << std::setfill('0') << std::setw(3) << *bit/pi*180.;
				rendercube.writegoftcube(ss.str());
				ss.str("");
			}
		return rendercube;
	}
}

