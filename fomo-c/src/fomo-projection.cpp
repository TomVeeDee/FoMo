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

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
const double pi=M_PI; //pi

FoMo::RenderCube projectioninterpolation(FoMo::GoftCube goftcube, const double l, const double b, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width)
{
	//FoMo::DataCube goftcube=object.datacube;
	FoMo::tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();
	int dim=goftcube.readdim();

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	std::cout << "Rotating coordinates to POS reference... " << std::flush;
	std::vector<double> xacc, yacc, zacc, losvel;
	xacc.resize(ng);
	yacc.resize(ng);
	zacc.resize(ng);
	losvel.resize(ng);
	// Define the unit vector along the line-of-sight
	std::vector<double> unit = {sin(b)*cos(l), -sin(b)*sin(l), cos(b)};
	// Read the physical variables
	FoMo::tphysvar peakvec=goftcube.readvar(0);//Peak intensity 
	FoMo::tphysvar fwhmvec=goftcube.readvar(1);// line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);  
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);
	
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	for (int i=0; i<ng; i++)
	{
		std::vector<double> gridpoint; // declare gridpoint here so that it is private to the thread
		gridpoint.resize(dim);
		for (int j=0; j<dim; j++)	gridpoint[j]=grid[j][i];
		xacc[i]=gridpoint[0]*cos(b)*cos(l)-gridpoint[1]*cos(b)*sin(l)-gridpoint[2]*sin(b);// rotated grid
		yacc[i]=gridpoint[0]*sin(l)+gridpoint[1]*cos(l);
		zacc[i]=gridpoint[0]*sin(b)*cos(l)-gridpoint[1]*sin(b)*sin(l)+gridpoint[2]*cos(b);
		std::vector<double> velvec = {vx[i], vy[i], vz[i]};// velocity vector
		double losvelval = inner_product(unit.begin(),unit.end(),velvec.begin(),0.0);//velocity along line of sight for position [i]/[ng]
		losvel.at(i)=losvelval;
	}
	
	// compute the bounds of the input data points, so that we can equidistantly distribute the target pixels
	double minx=*(min_element(xacc.begin(),xacc.end()));
	double maxx=*(max_element(xacc.begin(),xacc.end()));
	double miny=*(min_element(yacc.begin(),yacc.end()));
	double maxy=*(max_element(yacc.begin(),yacc.end()));
	double minz=*(min_element(zacc.begin(),zacc.end()));
	double maxz=*(max_element(zacc.begin(),zacc.end()));
	std::cout << "Done!" << std::endl;

	std::string chiantifile=goftcube.readchiantifile();
	double lambda0=goftcube.readlambda0();// lambda0=AIA bandpass for AIA imaging
       	
	std::cout << "Building frame: " << std::flush;
	double x,y,lambdaval;
	int ind;
	boost::progress_display show_progress(ng);
	
	//initialize grids
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel),yvec(x_pixel*y_pixel*lambda_pixel),lambdavec(x_pixel*y_pixel*lambda_pixel);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	if (lambda_pixel > 1) newgrid.push_back(lambdavec);
	FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel,0);
	
	// initialise the rendering
#ifdef _OPENMP
#pragma omp parallel for private (x,y,lambdaval,ind)
#endif
	for (int i=0; i<y_pixel; i++)
		for (int j=0; j<x_pixel; j++)
		{
			// now we're on one ray, through point with coordinates in the image plane
			x = double(j)/(x_pixel-1)*(maxx-minx)+minx;
			y = double(i)/(y_pixel-1)*(maxy-miny)+miny;
						
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
					ind=(i*(x_pixel)+j)*lambda_pixel+il;// 
					newgrid.at(0).at(ind)=x;
					newgrid.at(1).at(ind)=y;
					newgrid.at(2).at(ind)=lambdaval+lambda0; // store the full wavelength
				}
			}
			
			if (lambda_pixel==1) // AIA imaging study. Algorithm not verified [DY 14 Nov 2014]
			{
				ind=(i*x_pixel+j); 
				newgrid.at(0).at(ind)=x;
				newgrid.at(1).at(ind)=y;
			}
		}

	int i,j;
	double tempintens;
	// we step through the data points, and add their emissivity to the correct pixel
#ifdef _OPENMP
#pragma omp parallel for private(lambdaval,ind,i,j,tempintens) shared(xacc,yacc,peakvec,losvel,fwhmvec,intens)
#endif
	for (int k=0; k<ng; k++)
	{
		// xacc.at(k) contains x coordinate of pixel to be added
		j=std::round((xacc.at(k)-minx)*(x_pixel-1)/(maxx-minx));
		// yacc.at(k) contains y coordinate of pixel to be added
		i=std::round((yacc.at(k)-miny)*(y_pixel-1)/(maxy-miny));
		
		if (lambda_pixel>1)// spectroscopic study
		{
			for (int il=0; il<lambda_pixel; il++) // changed index from global variable l into il [D.Y. 17 Nov 2014]
			{
				// lambda the relative wavelength around lambda0, with a width of lambda_width
				lambdaval=static_cast<double>(il)/(lambda_pixel-1)*lambda_width-lambda_width/2.;
				tempintens=peakvec.at(k)*exp(-pow(lambdaval-losvel.at(k)/speedoflight*lambda0,2)/pow(fwhmvec.at(k),2)*4.*log(2.));
				ind=(i*(x_pixel)+j)*lambda_pixel+il;// 
#ifdef _OPENMP
#pragma omp atomic
#endif
				intens.at(ind)+=tempintens;// loop over z and lambda [D.Y 17 Nov 2014]
			}
		}
		
		if (lambda_pixel==1) // AIA imaging study. Algorithm not verified [DY 14 Nov 2014]
		{
			ind=(i*x_pixel+j); 
#ifdef _OPENMP
#pragma omp atomic
#endif
			intens.at(ind)+=peakvec.at(k);
		}
		// print progress
		++show_progress;
	}
	std::cout << " Done! " << std::endl << std::flush;
	
	FoMo::RenderCube rendercube(goftcube);
	FoMo::tvars newdata;
	double pathlength=(maxz-minz)/(z_pixel-1); // we set the pathlength to this value, corresponding to other integration methods. However, here it is not correct, and relies on the user to set an appropriate value for the number of pixels along the LOS.
	std::cout << "******Warning! Multiplying the emission with a guessed path length in the voxels of " << pathlength << "Mm" <<std::endl;
	std::cout << "******This is because you specified that there are approximately z_pixel " << z_pixel << " voxels along the largest dimension of the datacube." << std::endl;
	std::cout << "******Take care to correct this value if you need absolute intensities!" << std::endl << std::flush;
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
	FoMo::RenderCube RenderWithProjection(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		FoMo::RenderCube rendercube(goftcube);
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit)
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit)
			{
				rendercube=projectioninterpolation(goftcube,*lit,*bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
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

