#include "../config.h"
#include "FoMo.h"
#include "FoMo-internal.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <gsl/gsl_const_mksa.h>
#include <boost/progress.hpp>

#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#if defined(CGAL_LINKED_WITH_TBB) && CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,4,0) 
typedef CGAL::Triangulation_data_structure_3< 
    CGAL::Triangulation_vertex_base_3<K>, 
    CGAL::Triangulation_cell_base_3<K>, 
    CGAL::Parallel_tag>                                       Tds;
// no fast_location if openmp is used :(
//typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>		    Delaunay_triangulation_3;
typedef CGAL::Delaunay_triangulation_3<K, Tds>		    Delaunay_triangulation_3;
#else
// typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_triangulation_3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay_triangulation_3;
#endif

const double speedoflight=GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light
const double pi=M_PI; //pi

Delaunay_triangulation_3 triangulationfromdatacube(FoMo::DataCube goftcube)
{
	typedef K::Point_3                                    Point;
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	FoMo::tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();
	std::vector<Point> delaunaygrid;
	delaunaygrid.resize(ng);

	for (int i=0; i<ng; i++)
	{
		delaunaygrid[i]=Point(grid[0][i],grid[1][i],grid[2][i]);
	}

	// compute the Delaunay triangulation
	if (commrank==0) std::cout << "Doing Delaunay triangulation for interpolation onto rays... " << std::flush;
	Delaunay_triangulation_3 DT;
	// The triangulation should go quicker if it is sorted
	// CGAL::spatial_sort(delaunaygrid.begin(),delaunaygrid.end());
	// but I don't know how this affects the values in the maps.
	// Apparently, this is already done internally.
#if defined(CGAL_LINKED_WITH_TBB) && CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,4,0) 
	double minz=*(min_element(grid[2].begin(),grid[2].end()));
	double maxz=*(max_element(grid[2].begin(),grid[2].end()));
	double minx=*(min_element(grid[0].begin(),grid[0].end()));
	double maxx=*(max_element(grid[0].begin(),grid[0].end()));
	double miny=*(min_element(grid[1].begin(),grid[1].end()));
	double maxy=*(max_element(grid[1].begin(),grid[1].end()));
	Delaunay_triangulation_3::Lock_data_structure locking_ds(CGAL::Bbox_3(minx, miny, minz, maxx, maxy, maxz), 50);
	DT.insert(delaunaygrid.begin(),delaunaygrid.end());
#else	
	DT.insert(delaunaygrid.begin(),delaunaygrid.end());
#endif
	if (commrank==0) std::cout << "Done!" << std::endl << std::flush;
	return DT;
}

FoMo::RenderCube CGALinterpolation(FoMo::GoftCube goftcube, Delaunay_triangulation_3* DTpointer, const double l, const double b, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width)
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
	//FoMo::DataCube goftcube=object.datacube;
	FoMo::tgrid grid = goftcube.readgrid();
	int ng=goftcube.readngrid();

	// We will calculate the maximum image coordinates by projecting the grid onto the image plane
	// Rotate the grid over an angle -l (around z-axis), and -b (around y-axis)
	// Take the min and max of the resulting coordinates, those are coordinates in the image plane
	if (commrank==0) std::cout << "Rotating coordinates to POS reference... " << std::flush;
	std::vector<double> xacc, yacc, zacc;
	xacc.resize(ng);
	yacc.resize(ng);
	zacc.resize(ng);
	std::vector<double> gridpoint;
	gridpoint.resize(3);
	Point temporarygridpoint;
	// Define the unit vector along the line-of-sight
	std::vector<double> unit = {sin(b)*cos(l), -sin(b)*sin(l), cos(b)};
	// Read the physical variables
	FoMo::tphysvar peakvec=goftcube.readvar(0);//Peak intensity 
	FoMo::tphysvar fwhmvec=goftcube.readvar(1);// line width, =1 for AIA imaging
	FoMo::tphysvar vx=goftcube.readvar(2);  
	FoMo::tphysvar vy=goftcube.readvar(3);
	FoMo::tphysvar vz=goftcube.readvar(4);
	double losvelval;

// No openmp possible here
// Because of the insertions at the end of the loop, we get segfaults :(
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
	for (int i=0; i<ng; i++)
	{
		for (int j=0; j<3; j++)	gridpoint[j]=grid[j][i];
		xacc[i]=gridpoint[0]*cos(b)*cos(l)-gridpoint[1]*cos(b)*sin(l)-gridpoint[2]*sin(b);// rotated grid
		yacc[i]=gridpoint[0]*sin(l)+gridpoint[1]*cos(l);
		zacc[i]=gridpoint[0]*sin(b)*cos(l)-gridpoint[1]*sin(b)*sin(l)+gridpoint[2]*cos(b);
		temporarygridpoint=Point(grid[0][i],grid[1][i],grid[2][i]); //position vector
		// also create the map function_values here
		std::vector<double> velvec = {vx[i], vy[i], vz[i]};// velocity vector
		losvelval = inner_product(unit.begin(),unit.end(),velvec.begin(),0.0);//velocity along line of sight for position [i]/[ng]
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
	xacc.clear(); // release the memory
	yacc.clear();
	zacc.clear();
	if (commrank==0) std::cout << "Done!" << std::endl;

	std::string chiantifile=goftcube.readchiantifile();
	double lambda0=goftcube.readlambda0();// lambda0=AIA bandpass for AIA imaging
	double lambda_width_in_A=lambda_width*lambda0/speedoflight;
       	
	if (commrank==0) std::cout << "Building frame: " << std::flush;
	double x,y,z,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens;
	int li,lj,ind;
	Point p,nearest;
	Delaunay_triangulation_3::Vertex_handle v;
	Delaunay_triangulation_3::Locate_type lt;
	Delaunay_triangulation_3::Cell_handle c_old,c_new;
	bool could_lock_zone=false;
	boost::progress_display show_progress(x_pixel*y_pixel*z_pixel);
	
	//initialize grids
	FoMo::tgrid newgrid;
	FoMo::tcoord xvec(x_pixel*y_pixel*lambda_pixel),yvec(x_pixel*y_pixel*lambda_pixel),lambdavec(x_pixel*y_pixel*lambda_pixel);
	newgrid.push_back(xvec);
	newgrid.push_back(yvec);
	if (lambda_pixel > 1) newgrid.push_back(lambdavec);
	FoMo::tphysvar intens(x_pixel*y_pixel*lambda_pixel,0);
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(2) private (x,y,z,p,lt,li,lj,v,nearest,intpolpeak,intpolfwhm,intpollosvel,lambdaval,tempintens,ind,could_lock_zone,c_old,c_new)
#endif
	for (int i=0; i<y_pixel; i++)
		for (int j=0; j<x_pixel; j++)
		{
			#ifdef _OPENMP
			#pragma omp task
			#endif
			for (int k=0; k<z_pixel; k++) // scanning through ccd
			{
				x = double(j)/(x_pixel-1)*(maxx-minx)+minx;
				y = double(i)/(y_pixel-1)*(maxy-miny)+miny;
				z = double(k)/(z_pixel-1)*(maxz-minz)+minz;
		// calculate the interpolation in the original frame of reference
		// i.e. derotate the point using angles -l and -b
				p={x*cos(b)*cos(l)+y*sin(l)+z*sin(b)*cos(l),-x*cos(b)*sin(l)+y*cos(l)-z*sin(b)*sin(l),-x*sin(b)+z*cos(b)};
				
				#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,4,0) 
				// This is the proper way of doing a parallel location, in order to avoid collisions
				// Should we release the locks somehow?
				while (!could_lock_zone) 
				{
					c_new=DTpointer->locate(p, lt, li, lj, c_old, &could_lock_zone);
				}
				c_old=c_new;
				could_lock_zone=false;
				#else
				c_new=DTpointer->locate(p, lt, li, lj, c_old);
				c_old=c_new;
				#endif
				
		// Only look for the nearest point and interpolate, if the point p is inside the convex hull.
				if (lt!=Delaunay_triangulation_3::OUTSIDE_CONVEX_HULL)
				{
					v=DTpointer->nearest_vertex_in_cell(p,c_new);
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
						newgrid[0][ind]=x;
						newgrid[1][ind]=y;
						intens[ind]+=tempintens; // loop over z [D.Y 17 Nov 2014]
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
	rendercube.setrendermethod("CGAL");
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
	FoMo::RenderCube RenderWithCGAL(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, std::string outfile)
	{
		/* A good speedup would be to calculate the triangulation per ray.
		 * It would be good to select only the points around the ray, make the triangulation of that.
		 * The number of points would be drastically reduced, and the triangulation would be greatly sped up.
		 */
		Delaunay_triangulation_3 DT=triangulationfromdatacube(goftcube);
		FoMo::RenderCube rendercube(goftcube);
		for (std::vector<double>::iterator lit=lvec.begin(); lit!=lvec.end(); ++lit)
			for (std::vector<double>::iterator bit=bvec.begin(); bit!=bvec.end(); ++bit)
			{
				rendercube=CGALinterpolation(goftcube,&DT,*lit,*bit, x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
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
#endif
