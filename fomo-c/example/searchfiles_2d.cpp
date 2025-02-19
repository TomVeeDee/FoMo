#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <map>
#include <cmath>


using namespace std;

static const string rootdir="/users/cpa/tomvd/data/idl/FoMo/chiantitables/";

static const map<string,string>::value_type AIAfilters[]=
{
	map<string,string>::value_type("193",rootdir+"goft_table_aia193_abco.dat"),
	map<string,string>::value_type("094",rootdir+"goft_table_aia094_abco.dat"),
	map<string,string>::value_type("131",rootdir+"goft_table_aia131_abco.dat"),
	map<string,string>::value_type("171",rootdir+"goft_table_aia171_abco.dat"),
	map<string,string>::value_type("211",rootdir+"goft_table_aia211_abco.dat"),
	//map<string,string>::value_type("304",rootdir+"goft_table_aia304_abco.dat"),
	map<string,string>::value_type("335",rootdir+"goft_table_aia335_abco.dat"),
	map<string,string>::value_type("novalue","novalue")
};

static const map<string,string> AIAmap{ &AIAfilters[0], &AIAfilters[6]};

int main(int argc, char* argv[])
{
	// search the current path for files that contain string compstring
	// boost::filesystem::path full_path( boost::filesystem::initial_path<boost::filesystem::path>() );
	string full_path = "./";
	std::filesystem::directory_iterator end_itr; // default construction yields past-the-end
	vector<string> filelist;
	vector<string> complist;
	int ng=512;
	string compstring="datcubes_oscp005la8d512v0_2dzx_0000";
	complist.push_back(compstring);
	compstring=".dat";
	complist.push_back(compstring);
	for ( auto & itr : std::filesystem::directory_iterator{full_path} )
	{
		if ( std::filesystem::is_regular_file( itr ) )
		{
			string filename=itr.path();
			bool matches=true;
			for (unsigned int i=0; i<complist.size(); i++)
			{
				matches = ( matches && ( filename.find(complist[i]) != std::string::npos ) );
			}
			if ( ( matches == true) && ( filename.find(".fomo") == std::string::npos ) )
			{
				filelist.push_back( filename );
				// cout << "frame " << filelist.size() << " read from " << filename << endl;
			}
		}
	}
	int nframes=filelist.size();
	// Let's first get the code working for a single frame. Later on, we can still change to more frames, although that 
	//nframes=1;

	double pi=4*atan(1.);
	vector<double> langles={0.,pi/4,pi/2};
	vector<double> bangles={0.};
	int dim=2;	
	
	// Initialize the FoMo object
	FoMo::FoMoObject Object(dim);
	
	for (int t=0; t<nframes; t++)
	{
		string filename=filelist[t];
		cout << "Doing timestep " << t+1 << " of " << nframes << ": read from " << filename << endl << flush;
		
		// read the file
		ifstream in(filename,ios::binary);
		if (in.is_open())
		{
			string test;
			double tmpvar;
			vector<double> xvec,yvec;
			in >> test;
			for (int i=0; i<ng; i++)
			{
				in >> tmpvar;
				xvec.push_back(tmpvar);
			}
			
			for (int i=0; i<ng; i++)
			{
				in >> tmpvar;
				yvec.push_back(tmpvar);
			}
			
			// build x vector and y vector for loading into FoMo
			FoMo::tcoord x,y;
			double xval;
			int nz=1;
			for (int k=0; k<nz; k++)
			for (int i=0; i<ng; i++)
			{
				xval=xvec[i];
				for (int j=0; j<ng; j++)
				{
					x.push_back(xval);
					y.push_back(yvec[j]);
					//z.push_back(0.);
				}
			}
			FoMo::tgrid grid;
			grid.push_back(x);
			// it's better to push back y=0 as the second coordinate, 
			// so that we're modelling observations in the x,z plane
			// l is thus to be taken as 0
			// -b is the angle with the z-axis in the horizontal plane
			// grid.push_back(z);
			grid.push_back(y);
			
			FoMo::tvars vars;
			FoMo::tphysvar variable;
			int nvars=4;
			// read in rho, T, vx, vy
			for (int k=0; k<nvars; k++)
			{
				for (int i=0; i<ng; i++)
					for (int j=0; j<ng; j++)
					{
						in >> tmpvar;
						variable.push_back(tmpvar);
					}
				vars.push_back(variable);
				variable.clear();
			}
			
			Object.setdata(grid,vars);
			// data is in structure, now start the rendering
			int x_pixel=500;
			int y_pixel=500;
			int z_pixel=1;
			int lambda_pixel=1;
			double lambda_width=200000;
			Object.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
			Object.setrendermethod("CGAL2D");
			Object.setobservationtype(FoMo::Imaging);
			for (map<string,string>::const_iterator it=AIAmap.begin(); it!=AIAmap.end(); ++it)
			{
				Object.setchiantifile(it->second);
				cout << "Now rendering with the " << it->first << "A filter." << endl << flush;
				stringstream ss;
				ss << filename << ".fomo_" << it->first;
				Object.setoutfile(ss.str());
				Object.render(langles,bangles);
			}
		}
	}
}
