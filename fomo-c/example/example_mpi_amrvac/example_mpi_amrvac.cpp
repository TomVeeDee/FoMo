#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>

using namespace std;

double L_unit = 1e8; // length normalisation value
double L_unitFoMo = 1e6; // (express length unit in Mm for FoMo)

double rho_unit = 1e-12;
double Teunit = 1e+6;
double kboltz = 1.38065e-23;   //                 # J/K = kg/K m_^2/s_^2 (boltzmann const)
double mhydro = 1.6726e-27;    //                # kg  (proton mass for the amrvac_nonAMR.par)
double R_spec = kboltz/(0.5e0*mhydro);
double n_unit = rho_unit * 1.204 * 1.e21;
double V_unit = sqrt(R_spec*Teunit); //         velocity
double p_unit = rho_unit*pow(V_unit,2); //              pressure
// p = R_spec * rho * T -> T = p/R_spec/rho 

void loopthroughleafs(vector<bool> forest, int & forestposition, int & level, const int ndim, vector<vector<int>> & block_info, vector<int> index)
{
	forestposition++;
	if (forestposition<50) cout << forest.at(forestposition) << " " << forestposition << " " << level << endl;
	if (forest.at(forestposition))
	{
		// now we are a leaf
		vector<int> tempblockinfo(4);
		tempblockinfo.at(0)=level;
		for (int i=0; i<ndim; i++)
			tempblockinfo.at(i+1)=index.at(i);
		block_info.push_back(tempblockinfo);
	}
	else
	{
		// we are a parent, loop through the children
		level++;
		vector<int> newindex(ndim);
		for (int i=0; i<pow(2,ndim); i++)
		{
			if (ndim==2) newindex = {2*(index.at(0)-1) + 1 + i%2,2*(index.at(1)-1) + 1 + i/2,1};
			if (ndim==3) newindex = {2*(index.at(0)-1) + 1 + i%2,2*(index.at(1)-1) + 1 + (i/2)%2,2*(index.at(2)-1) + 1 + i/4};
			loopthroughleafs(forest, forestposition, level, ndim, block_info, newindex);
		}
	}
}

int main(int argc, char* argv[])
{
	// First check to see if two filenames are given as argument.
	if (argc != 3)
	{
		cout << "This program should be run with two file names as argument." << endl;
		cout << "The first file is the AMRVAC .dat file to be rendered, the second one is the amrvac.par file." << endl;
		exit(EXIT_FAILURE);
	}
	
	// Initialize the FoMo object
		FoMo::FoMoObject Object;
	
	// Open the file in the argument as ifstream
	// it is a binary file, and we open it at the end (because we first read the end)
	ifstream file(argv[1],ios::in|ios::binary|ios::ate);
	
	if (file.is_open())
	{
		// the end of the AMRVAC dat file has nleafs, levmax, ndim, ndir, nw, neqpar+nspecialpar, it, t
		int nleafs, levmax, ndim, ndir, nw, neqpar, it;
		double t;
		file.seekg(-7*sizeof(nleafs)-sizeof(t),ios_base::cur);
		file.read(reinterpret_cast<char*>(&nleafs), sizeof(int));
		file.read(reinterpret_cast<char*>(&levmax), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndim), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndir), sizeof(int));
		file.read(reinterpret_cast<char*>(&nw), sizeof(int));
		cout << "The simulation contains " << ndim << "D data, with " << nw << " variables. It has " << nleafs << " top-level blocks (=leafs)." << endl;
		file.read(reinterpret_cast<char*>(&neqpar), sizeof(int));
		file.read(reinterpret_cast<char*>(&it), sizeof(int));
		file.read(reinterpret_cast<char*>(&t), sizeof(double));
		//go back to where we started reading
		file.seekg(-7*sizeof(int)-sizeof(double),ios_base::cur);
		// before this information, a vector of length neqpar+nspecialpar was written, and the number of points in each dimension 
		// we can just read neqpar times a double, since we do not need to store it
		double tmpeqpar;
		vector<int> nx(ndim);
		// nglev1 contains the number of cells per block
		int nglev1=1;
		file.seekg(-neqpar*sizeof(double)-ndim*sizeof(int),ios_base::cur);
		// this position is the end of the forest
		auto endofforest=file.tellg();
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&nx.at(i)), sizeof(int));
			nglev1*=nx.at(i);
		}
		for (int i=0; i<neqpar; i++)
		{
			file.read(reinterpret_cast<char*>(&tmpeqpar), sizeof(double));
		}
		
		// let's read the data blocks, they are at the beginning of the file
		file.seekg(0,ios_base::beg);
		double temp;
		vector<double> block(nw*nglev1);
		vector<vector<double>> leafs;
		for (int i=0; i<nleafs; i++)
		{
			// read nleaf blocks
			// blocks are of size nw*nglev1*sizeof(double)
			for (int j=0; j<nw; j++)
				for (int k=0; k<nglev1; k++)
				{
					file.read(reinterpret_cast<char*>(&temp), sizeof(double));
					block.at(j+k*nw)=temp;
				}
			leafs.push_back(block);
		}
		// this position is the start of the forest
		auto startofforest=file.tellg();
		
		// now we read the so-called forest, which indicates which blocks are leafs or parents in AMRVAC's block structure
		// bool is stored as int in the unformatted fortran file
		int tempint;
		vector<bool> forest((endofforest-startofforest)/sizeof(int));
		for (unsigned int i=0; i<forest.size(); i++)
		{
			file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
			if (i<50) cout << tempint << endl;
			forest.at(i)=tempint;
		}
		file.close();
		
		// to convert the forest to physical units, we need some data from the amrvac.par file
		vector<int> nxlone(ndim);
		vector<double> xprobmin(ndim),xprobmax(ndim);
		ifstream vacpar(argv[2],ios_base::ate);
		if (vacpar.is_open())
		{
			// read the amrvac.par file in a buffer as a string
			streamsize vacpar_size = vacpar.tellg();
			vacpar.seekg(0,ios_base::beg);
			vector<char> buffer(vacpar_size);
			vacpar.read(buffer.data(),vacpar_size);
			string stringbuffer(buffer.begin(),buffer.end());
			
			// iterate over keywords to be found: "nxlone", "xprobmin", "xprobmax"
			vector<string> keywords({"nxlone","xprobmin","xprobmax"});
			
			// search for the values of keyword1, keyword2, keyword3 (if 3D)
			for (unsigned int j=0; j<keywords.size(); j++)
			for (int i=0; i<ndim; i++)
			{
				stringstream key;
				key << keywords.at(j) << i+1;
				// find the first occurrence of keyword1
				size_t result=stringbuffer.find(key.str());
				if (result == string::npos)
				{
					cerr << argv[2] << " does not contain the key " << key.str();
					exit (EXIT_FAILURE);
				}
				// go to the start of the numeric value after the equal sign (i.e. skip blanks or =)
				// this will fail if no blank or = is found! should we check for this or assume that the amrvac.par file is well formed?
				size_t start=stringbuffer.find_first_of("\t\n= ",result+1);
				// convert the string to integer or double, and store it in the appropriate variable
				switch(j){
					case 0: 
					{
						nxlone.at(i)=stoi(stringbuffer.substr(start+1));
						cout << argv[2] << ": " << key.str() << "=" << nxlone.at(i) << endl;
						break;
					}
					case 1:
					{
						xprobmin.at(i)=stod(stringbuffer.substr(start+1));
						cout << argv[2] << ": " << key.str() << "=" << xprobmin.at(i) << endl;
						break;
					}
					case 2:
					{
						xprobmax.at(i)=stod(stringbuffer.substr(start+1));
						cout << argv[2] << ": " << key.str() << "=" << xprobmax.at(i) << endl;
						break;
					}
				}
			}
			vacpar.close();
		}
		else 
		{
			cerr << "amrvac.par file not found at " << argv[2] << endl;
			exit (EXIT_FAILURE);
		}
		
		//compute the number of blocks in the simulation
		vector<int> nblocks(ndim);
		//compute the dx in each direction
		vector<double> cellsize(ndim);
		for (int i=0; i<ndim; i++)
		{
			nblocks.at(i)=nxlone.at(i)/nx.at(i);
			cellsize.at(i)=(xprobmax.at(i)-xprobmin.at(i))/nxlone.at(i);
			// if the number of grid points is not an integer multiple of nx, then there is a problem
			if (nxlone.at(i)%nx.at(i)!=0)
			{
				cerr << "Non-integer number of blocks. Are you sure " << argv[2] << " contains the correct amrvac.par file?";
				exit (EXIT_FAILURE);
			}
		}
		
		// now go through the forest and put the correct leaf blocks in the FoMoObject
		vector<vector<int>> block_info;
		int forestposition = -1;
		for (int k=0; k<nblocks.at(2); k++)
		for (int j=0; j<nblocks.at(1); j++)
		for (int i=0; i<nblocks.at(0); i++)
		{
			int level = 1;
			loopthroughleafs(forest, forestposition, level, ndim, block_info, {i+1,j+1,k+1});
		}
		if (block_info.size() != nleafs)
		{
			cerr << "block_info.size() is " << block_info.size() << " and does not match nleafs=" << nleafs << endl;
			exit (EXIT_FAILURE);
		}
		
		// loop through nleafs blocks and load into FoMoObject
		for (int i=0; i<nleafs; i++)
		{
			// calculate coordinates of bottom left of block
			vector<double> bottomleft;
			// local cell size
			vector<double> localcellsize;
			for (int j=0; j<ndim; j++)
			{
				bottomleft.push_back(xprobmin.at(j)+(block_info.at(i).at(1+j) - 1)*cellsize.at(j)*nx.at(j)/pow(2,(block_info.at(i).at(0)-1)));
				localcellsize.push_back(cellsize.at(j)/pow(2,(block_info.at(i).at(0)-1)));
				//if (i<200 && j==0) cout << block_info.at(i).at(0) << " " << block_info.at(i).at(j+1) << " " << bottomleft.at(j) << " " << cellsize.at(j)*nx.at(j)/pow(2,(block_info.at(i).at(0)-1)) << endl;
			}
			
			// loop over nglev1 cells in block
			for (int k=0; k<nglev1; k++)
			{
				// calculate coordinates of cell in the block
				vector<int> localcoord;
				localcoord.push_back(1+k%nx.at(0));
				localcoord.push_back(1+k/nx.at(0)%nx.at(0));
				if (ndim == 3) localcoord.push_back(1+k/nx.at(0)/nx.at(1));
				vector<double> coordinates;
				for (int j=0; j<ndim; j++)
				{
					coordinates.push_back(bottomleft.at(j)+(localcoord.at(j)-.5)*localcellsize.at(j));
					coordinates.at(j)*=L_unit/L_unitFoMo; // convert all lengths to Mm 
					//if (k==0 && i==0) cout << localcellsize.at(j) << " " << bottomleft.at(j) << " " ;
					//cout << coordinates.at(j) << " ";
					//if (j==ndim-1) cout << endl;
				}
				
				vector<double> variables;
				// rho
				variables.push_back(leafs.at(i).at(0+k*nw)); // 0 here corresponds to which variable needs to be loaded into FoMo
				variables.at(0)*=n_unit;
				// T
				variables.push_back(leafs.at(i).at(4+k*nw)/leafs.at(i).at(0+k*nw));
				//variables.at(1)*=p_unit/R_spec/rho_unit; 
				variables.at(1)*=1e4;
				// vx
				variables.push_back(leafs.at(i).at(1+k*nw));
				variables.at(2)*=V_unit;
				// vy
				variables.push_back(leafs.at(i).at(2+k*nw));
				variables.at(3)*=V_unit;
				// vz
				variables.push_back(leafs.at(i).at(3+k*nw));
				variables.at(4)*=V_unit;
				//if (k==0 && i==0) cout << R_spec << " " <<  p_unit/R_spec/rho_unit << " " << variables.at(1) << endl;
				// cout << "rho " << variables.at(0) << " T " << variables.at(1) << " v " << variables.at(2) << " " << variables.at(3) << " " << variables.at(4) << endl; 
				
				// load data into FoMo
				Object.push_back_datapoint(coordinates, variables);
			}
		}
	}
	else
	{
		cerr << ".dat file not found at " << argv[1] << endl;
		exit (EXIT_FAILURE);
	}
	
	// data is in structure, now start the rendering
	
	/// [Set rendering options]
	Object.setchiantifile("../../../chiantitables/goft_table_fe_12_0194_abco.dat"); // the default value is "../chiantitables/goft_table_fe_12_0194small_abco.dat"
	Object.setabundfile("/empty"); //use "/empty" or do not set it at all for the default sun_coronal_2012_schmelz.abund file
	Object.setrendermethod("NearestNeighbour"); // NearestNeighbour is the default rendermethod
	Object.setobservationtype(FoMo::Spectroscopic);
	// adjust the resolution with these parameters
	int x_pixel=149;
	int y_pixel=148;
	int z_pixel=500;
	int lambda_pixel=100;
	double lambda_width=200000; // in m/s
	Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
	// determine where the output will be written
	Object.setoutfile("fomo-example-out.");
	/// [Set rendering options]
	
	/// [Render]
	Object.render(2*atan(1.),2*atan(1.));
	// alternatively, you could add different angles in radians as argument, e.g. Object.render(1.57/2.,1.57/2.) to obtain some nice Doppler shifts.
	/// [Render]
	
	//FoMo::GoftCube goftcube=Object.readgoftcube();
	//string outfile("goftcube.txt");
	//goftcube.writegoftcube(outfile);
}
