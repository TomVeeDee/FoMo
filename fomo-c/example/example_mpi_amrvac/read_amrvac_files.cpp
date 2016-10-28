#include "FoMo.h"
#include "FoMo-amrvac.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <gsl/gsl_const_mksa.h>

using namespace std;

double kboltz = GSL_CONST_MKSA_BOLTZMANN;   //                 # J/K = kg/K m_^2/s_^2 (boltzmann const)
double mhydro = GSL_CONST_MKSA_MASS_PROTON;    //                # kg  (proton mass for the amrvac_nonAMR.par)
double R_spec = kboltz/(0.5e0*mhydro);

// this recursive procedure loops through the forest and calculates the indices of each leaf block
// the implementation is a carbon copy of the python code of Oliver Porth (included in AMRVAC)
void loopthroughleafs(vector<bool> forest, int & forestposition, int & level, const int ndim, vector<vector<int>> & block_info, vector<int> index)
{
	forestposition++;
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
		vector<int> newindex(ndim);
		for (int i=0; i<pow(2,ndim); i++)
		{
			if (ndim==2) newindex = {2*(index.at(0)-1) + 1 + i%2,2*(index.at(1)-1) + 1 + i/2,1};
			if (ndim==3) newindex = {2*(index.at(0)-1) + 1 + i%2,2*(index.at(1)-1) + 1 + (i/2)%2,2*(index.at(2)-1) + 1 + i/4};
			int tempint=level+1;
			loopthroughleafs(forest, forestposition, tempint, ndim, block_info, newindex);
		}
	}
}

// to convert the forest to physical units, we need some data from the amrvac.par file
// this is the routine to read it in
// the relevant information is stored in nxlone, xprobmin, xprobmax, vectors of length ndim
void read_amrvac_par_file(const char* amrvacpar, const int ndim, vector<int> & nxlone, vector<double> & xprobmin, vector<double> & xprobmax)
{
		nxlone.resize(ndim);
		xprobmin.resize(ndim);
		xprobmax.resize(ndim);
		ifstream vacpar(amrvacpar,ios_base::ate);
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
					cerr << amrvacpar << " does not contain the key " << key.str();
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
						cout << amrvacpar << ": " << key.str() << "=" << nxlone.at(i) << endl;
						break;
					}
					case 1:
					{
						xprobmin.at(i)=stod(stringbuffer.substr(start+1));
						cout << amrvacpar << ": " << key.str() << "=" << xprobmin.at(i) << endl;
						break;
					}
					case 2:
					{
						xprobmax.at(i)=stod(stringbuffer.substr(start+1));
						cout << amrvacpar << ": " << key.str() << "=" << xprobmax.at(i) << endl;
						break;
					}
				}
			}
			vacpar.close();
		}
		else 
		{
			cerr << "amrvac.par file not found at " << amrvacpar << endl;
			exit (EXIT_FAILURE);
		}
}

FoMo::FoMoObject read_amrvac_dat_file(const char* datfile, const char* amrvacpar, string amrvac_version, const int gamma_eqparposition, const double n_unit, const double Teunit, const double L_unit)
{
	// Initialize the FoMo object
	FoMo::FoMoObject Object;
	
	// Open the file in the argument as ifstream
	// it is a binary file, and we open it at the end (because we first read the end)
	ifstream file(datfile,ios::in|ios::binary|ios::ate);
	
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
		double tmpeqpar, gamma;
		vector<int> nx(ndim);
		vector<double> eqpar(neqpar);
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
			eqpar.at(i)=tmpeqpar;
		}
		gamma = eqpar.at(gamma_eqparposition);
		
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
			forest.at(i)=tempint;
		}
		file.close();
		
		vector<int> nxlone;
		vector<double> xprobmin, xprobmax;
		read_amrvac_par_file(amrvacpar, ndim, nxlone, xprobmin, xprobmax);
		
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
				cerr << "Non-integer number of blocks. Are you sure " << amrvacpar << " contains the correct amrvac.par file?";
				exit (EXIT_FAILURE);
			}
		}
		
		// now go through the forest and put the correct leaf blocks in the FoMoObject
		vector<vector<int>> block_info;
		vector<vector<int>> mortoncurve;
		// begin the loop through the leafs at forestposition 0 (forestposition++ at start of loopthroughleafs!)
		int forestposition = -1;
		// make sure the code is also working for 2D data!
		if (ndim<3) nblocks.push_back(1);
		if (amrvac_version.compare("old")==0)
		{
			for (int k=0; k<nblocks.at(2); k++)
			for (int j=0; j<nblocks.at(1); j++)
			for (int i=0; i<nblocks.at(0); i++)
			{
				int level = 1;
				loopthroughleafs(forest, forestposition, level, ndim, block_info, {i+1,j+1,k+1});
			}
		}
		else if (amrvac_version.compare("gitlab")==0)
		{
			// compute the maximum length of the morton curve
			int maxgridlength = *(max_element(nblocks.begin(),nblocks.end()));
			// resize the morton curve to the length of a hypothetic cube with sides maxgridlength
			// later on, we will prune simulation points from the mortoncurve
			// can this be optimized for 2D data?
			mortoncurve.resize(pow(maxgridlength,3));
			for (int k=0; k<maxgridlength; k++)
			for (int j=0; j<maxgridlength; j++)
			for (int i=0; i<maxgridlength; i++)
			{
				mortoncurve.at(mortonEncode_for(i,j,k))={i+1,j+1,k+1};
			}
			// loop through mortoncurve
			// at each block within the domain, store the leafs from that parent into block_info
			for (unsigned int i=0; i<mortoncurve.size(); i++)
			{
				// check if the block is within the domain. If not, just skip it and go to the next block
				if (mortoncurve.at(i).at(0)<=nblocks.at(0) && mortoncurve.at(i).at(1)<=nblocks.at(1) && mortoncurve.at(i).at(2)<=nblocks.at(2))
				{
					int level = 1;
					loopthroughleafs(forest, forestposition, level, ndim, block_info, mortoncurve.at(i));
				}
			}
		}
		else
		{
			cerr << "The stated AMRVAC implementation " << amrvac_version << " is not know." << endl;
			exit(EXIT_FAILURE);
		}
		// check if the resulting blocks correspond to the expected number of leafs
		if (block_info.size() != nleafs)
		{
			cerr << "block_info.size() is " << block_info.size() << " and does not match nleafs=" << nleafs << endl;
			exit (EXIT_FAILURE);
		}
		
		double V_unit = sqrt(R_spec*Teunit); //         velocity
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
					coordinates.at(j)*=L_unit; // convert all lengths to Mm 
				}
				
				vector<double> variables;
				// rho
				variables.push_back(leafs.at(i).at(0+k*nw)); // 0 here corresponds to which variable needs to be loaded into FoMo
				variables.at(0)*=n_unit;
				// T
				// first we calculate the pressure from the internal energy 
				// p = (5/3 -1)*(e- K-B), with K = 0.5*rho *(vi^2), and B = 0.5*(bi^2)
				double kineticenergy = inner_product(&leafs.at(i).at(1+k*nw), &leafs.at(i).at(3+k*nw), &leafs.at(i).at(1+k*nw), 0.)/leafs.at(i).at(0+k*nw)/2.;
				double magneticenergy = inner_product(&leafs.at(i).at(5+k*nw), &leafs.at(i).at(7+k*nw), &leafs.at(i).at(5+k*nw), 0.)/2.;
				double p = (gamma-1)*(leafs.at(i).at(4+k*nw)-kineticenergy - magneticenergy);
				// then T = p/rho*Teunit
				variables.push_back(p/leafs.at(i).at(0+k*nw));
				variables.at(1)*=Teunit; 
				// vx = m1/rho
				variables.push_back(leafs.at(i).at(1+k*nw)/leafs.at(i).at(0+k*nw));
				variables.at(2)*=V_unit;
				// vy = m2/rho
				variables.push_back(leafs.at(i).at(2+k*nw)/leafs.at(i).at(0+k*nw));
				variables.at(3)*=V_unit;
				// vz = m3/rho
				variables.push_back(leafs.at(i).at(3+k*nw)/leafs.at(i).at(0+k*nw));
				variables.at(4)*=V_unit;
				
				// load data into FoMo
				Object.push_back_datapoint(coordinates, variables);
			}
		}
	}
	else
	{
		cerr << ".dat file not found at " << datfile << endl;
		exit (EXIT_FAILURE);
	}

	return Object;
}
