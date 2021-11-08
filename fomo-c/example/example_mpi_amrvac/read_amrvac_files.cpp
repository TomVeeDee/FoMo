#include "FoMo.h"
#include "FoMo-amrvac.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <gsl/gsl_const_mksa.h>

using namespace std;

double kboltz = GSL_CONST_MKSA_BOLTZMANN;   //                 # J/K = kg/K m_^2/s_^2 (boltzmann const)
double mhydro = GSL_CONST_MKSA_MASS_PROTON;    //                # kg  (proton mass for the amrvac_nonAMR.par)
double R_spec = kboltz/(0.5e0*mhydro);

//Version 5 (AMRVAC 2.2) compatibility added by Krishna
//Version 4 may also be integrated here after testing
std::vector<int> compatible_versions={3, 5};  

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

void blocksandsize(std::vector<int> nxlone, std::vector<int> nx, std::vector<double> xprobmin,std::vector<double> xprobmax, std::vector<int> & nblocks, std::vector<double> & cellsize)
{
	int ndim=nxlone.size();
	for (int i=0; i<ndim; i++)
	{
		nblocks.at(i)=nxlone.at(i)/nx.at(i);
		cellsize.at(i)=(xprobmax.at(i)-xprobmin.at(i))/nxlone.at(i);
		// if the number of grid points is not an integer multiple of nx, then there is a problem
		if (nxlone.at(i)%nx.at(i)!=0)
		{
			cerr << "Non-integer number of blocks. Are you sure you passed the correct amrvac.par file?";
			exit (EXIT_FAILURE);
		}
	}
}

std::vector<bool> read_forest(istream & file, int forestsize)
{
	// now we read the so-called forest, which indicates which blocks are leafs or parents in AMRVAC's block structure
	// bool is stored as int in the unformatted fortran file
	int tempint;
	std::vector<bool> forest(forestsize);
	for (unsigned int i=0; i<forest.size(); i++)
	{
		file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
		forest.at(i)=tempint;
	}
	return forest;
}

std::vector<std::vector<double>> read_blocks(istream & file, const int ndim, const int nw, const int nglev1, const int nleafs, const int version)
{
	double tempdouble;
	int tempint;
	std::vector<double> block(nw*nglev1);
	std::vector<std::vector<double>> leafs;
	for (int i=0; i<nleafs; i++)
	{
		// read nleaf blocks
		// blocks are of size nw*nglev1*sizeof(double)
		// and are preceded by two ghost cell integers per dimension (if version is 3)
		if (version==3)
		{
			for (int j=0; j<ndim; j++)
			{
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
			}
		}

		//////////////////////////////////////////////////////////
		// added by vaibhav for amrvac datfiles version 4. Test version.
		//////////////////////////////////////////////////////////
		if (version==4)
		{
			for (int j=0; j<ndim; j++)
			{
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
			}
		}
		///////////////////////////////////////////////////////////////
		// added by Krishna for compatibility with AMRVAC 2.2 (version 5)
		///////////////////////////////////////////////////////////////
		if (version==5)
		{
			for (int j=0; j<ndim; j++)
			{
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
				file.read(reinterpret_cast<char*>(&tempint), sizeof(int));
			}
		}
		///////////////////////////////////////////////////////////////
		// then start reading the block
		// in principle, if the above number of ghost cells is non-zero, then we need to take a different block size, and select different numbers
		// implement this later, if needed
		// we should also differentiate between version numbers for this aspect, though
		for (int j=0; j<nw; j++)
			for (int k=0; k<nglev1; k++)
			{
				file.read(reinterpret_cast<char*>(&tempdouble), sizeof(double));
				block.at(j+k*nw)=tempdouble;
			}
		leafs.push_back(block);
	}

	return leafs;
}

std::vector<std::vector<int>> build_block_info_morton(std::vector<int> nblocks, std::vector<bool> forest, int ndim, int nleafs)
{
	std::vector<std::vector<int>> block_info;
	std::vector<std::vector<int>> mortoncurve;
	// begin the loop through the leafs at forestposition 0 (forestposition++ at start of loopthroughleafs!)
	int forestposition = -1;

	// compute the maximum length of the morton curve
	int maxgridlength = *(max_element(nblocks.begin(),nblocks.end()));
	// resize the morton curve to the length of a hypothetic cube with sides maxgridlength
	// later on, we will prune simulation points from the mortoncurve
	// can this be optimized for 2D data?
	// If number of blocks is odd, make cube with +1 size
	//if (maxgridlength % 2 !=  0) maxgridlength++; //removed by Krishna
	//it seems mortonEncode_for will produce indices of the same size as the number of grid cells only
	//when the dimensions of the grid is in powers of 2. So added the below line to accommodate this. 
	//perhaps not an efficient approach when maxgridlength is large. Should be optimized in the future.
        maxgridlength = pow(2, ceil(log(maxgridlength)/log(2))); //added by Krishna
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

	// check if the resulting blocks correspond to the expected number of leafs
	if (block_info.size() != nleafs)
	{
		cerr << "block_info.size() is " << block_info.size() << " and does not match nleafs=" << nleafs << endl;
		exit (EXIT_FAILURE);
	}

	return block_info;
}

FoMo::FoMoObject fomo_from_amrvac(const int ndim,const int nleafs,const int nglev1,const int nw, std::vector<double> xprobmin, std::vector<int> nx, std::vector<double> cellsize, std::vector<std::vector<int>> block_info, std::vector<std::vector<double>> leafs, const double n_unit, const double Teunit, const double L_unit, const double gamma)
{
	// could take a vector of necessary variables for data into fomo (from varnames in read_amrvac_new_dat_file, or fixed in read_amrvac_old_dat_file)

	// Initialize the FoMo object, and load AMRVAC data into it
	FoMo::FoMoObject Object;

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
			double kineticenergy = std::inner_product(&leafs.at(i).at(1+k*nw), &leafs.at(i).at(3+k*nw), &leafs.at(i).at(1+k*nw), 0.)/leafs.at(i).at(0+k*nw)/2.;
			double magneticenergy = std::inner_product(&leafs.at(i).at(5+k*nw), &leafs.at(i).at(7+k*nw), &leafs.at(i).at(5+k*nw), 0.)/2.;
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

	return Object;
}

FoMo::FoMoObject read_amrvac_old_dat_file(const char* datfile, const char* amrvacpar, string amrvac_version, const int gamma_eqparposition, const double n_unit, const double Teunit, const double L_unit)
{
	FoMo::FoMoObject Object;

	// Open the file in the argument as ifstream
	// it is a binary file
	ifstream filehandle(datfile,ios::in|ios::binary);

	if (filehandle.is_open())
	{
		// read the filehandle into a buffer, this saves around 10s for each 2GB file
	        stringstream file;
        	file << filehandle.rdbuf();
		// go to the end of the file to start reading
        	file.seekg(0,std::ios::end);

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
		std::vector<std::vector<double>> leafs=read_blocks(file,ndim,nw,nglev1,nleafs,0);
		// this position is the start of the forest
		auto startofforest=file.tellg();

		// read the forest
		std::vector<bool> forest=read_forest(file,(endofforest-startofforest)/sizeof(int));

		filehandle.close();

		vector<int> nxlone;
		vector<double> xprobmin, xprobmax;
		read_amrvac_par_file(amrvacpar, ndim, nxlone, xprobmin, xprobmax);

		//compute the number of blocks in the simulation
		vector<int> nblocks(ndim);
		//compute the dx in each direction
		vector<double> cellsize(ndim);
		blocksandsize(nxlone,nx,xprobmin,xprobmax,nblocks,cellsize);

		// now go through the forest and put the correct leaf blocks in the FoMoObject
		vector<vector<int>> block_info;
		// make sure the code is also working for 2D data!
		if (ndim<3) nblocks.push_back(1);
		if (amrvac_version.compare("old")==0)
		{
			// begin the loop through the leafs at forestposition 0 (forestposition++ at start of loopthroughleafs!)
			int forestposition = -1;

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
			block_info=build_block_info_morton(nblocks,forest,ndim,nleafs);
		}
		else
		{
			cerr << "The stated AMRVAC implementation " << amrvac_version << " is not know." << endl;
			exit(EXIT_FAILURE);
		}

		// Initialize the FoMo object, and load AMRVAC data into it
		Object=fomo_from_amrvac(ndim,nleafs,nglev1,nw,xprobmin,nx,cellsize,block_info,leafs,n_unit,Teunit,L_unit,gamma);
	}
	else
	{
		cerr << ".dat file not found at " << datfile << endl;
		exit (EXIT_FAILURE);
	}

	return Object;
}

int read_dat_version(const char* datfile)
{
	ifstream filehandle(datfile,ios::in|ios::binary);
	int version_number;

	if (filehandle.is_open())
	{
		filehandle.read(reinterpret_cast<char*>(&version_number), sizeof(int));
		filehandle.close();
		std::cout<<"This dat version is "<<version_number<<std::endl; //added by vaibhav.
	}
	return version_number;
}

FoMo::FoMoObject read_amrvac_new_dat_file(const char* datfile, int version_number, const double n_unit, const double Teunit, const double L_unit)
{
	// Initialize the FoMo object
	FoMo::FoMoObject Object;

	// Open the file in the argument as ifstream
	// it is a binary file
	ifstream filehandle(datfile,ios::in|ios::binary);

	if (filehandle.is_open())
	{
		// read the filehandle into a buffer, this saves around 10s for each 2GB file
		stringstream file;
		file << filehandle.rdbuf();
		filehandle.close();

		int nleafs, levmax, ndim, ndir, nw, neqpar, it;
		int version, treeoffset, blockoffset, nparents;
		double t;

		// follow specification at http://amrvac.org/md_doc_snapshot_format.html
		file.read(reinterpret_cast<char*>(&version), sizeof(int));
		file.read(reinterpret_cast<char*>(&treeoffset), sizeof(int));
		file.read(reinterpret_cast<char*>(&blockoffset), sizeof(int));
		file.read(reinterpret_cast<char*>(&nw), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndir), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndim), sizeof(int));
		file.read(reinterpret_cast<char*>(&levmax), sizeof(int));
		file.read(reinterpret_cast<char*>(&nleafs), sizeof(int));
		file.read(reinterpret_cast<char*>(&nparents), sizeof(int));
		cout << "The simulation contains " << ndim << "D data, with " << nw << " variables. It has " << nleafs << " top-level blocks (=leafs)." << endl;
		file.read(reinterpret_cast<char*>(&it), sizeof(int));
		file.read(reinterpret_cast<char*>(&t), sizeof(double));

		vector<int> nxlone(ndim),nx(ndim);
		vector<double> xprobmin(ndim), xprobmax(ndim);
		// block_nx is earlier nxblock
		// domain_nx is earlier nxlone
		// amrvac.par is not needed anymore!
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&xprobmin.at(i)), sizeof(double));
			std::cout << "xprobmin(" << i << "):" << xprobmin.at(i) << std::endl;
		}
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&xprobmax.at(i)), sizeof(double));
			std::cout << "xprobmax(" << i << "):" << xprobmax.at(i) << std::endl;
		}
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&nxlone.at(i)), sizeof(int));
			std::cout << "domain_nx(" << i << "):" << nxlone.at(i) << std::endl;
		}
		int nglev1=1;
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&nx.at(i)), sizeof(int));
			nglev1*=nx.at(i);
			std::cout << "block_nx(" << i << "):" << nx.at(i) << std::endl;
		}		
		///////////////////////////////////////////////////////////////
		// added by Krishna for compatibility with AMRVAC 2.2 (version 5)
		///////////////////////////////////////////////////////////////
		if (version == 5)
		{		
			std::vector<int> periodic(ndim);
		        for (int i=0; i<ndim; i++)
		        {
		                file.read(reinterpret_cast<char*>(&periodic.at(i)), sizeof(int));
		                std::cout << "Boundary Periodicity ["<< i <<"]:" << periodic.at(i)  << std::endl;
		        }
		        char geom[16];
		        file.read(geom,16);
		        std::cout << "Geometry:" << geom <<std::endl;
			int stag;
		        file.read(reinterpret_cast<char*>(&stag),sizeof(int));
		        std::cout << "Staggered:" << stag <<std::endl;
		}
		///////////////////////////////////////////////////////////////
		char temp[16];
		std::vector<std::string> varnames(nw);
		for (int i=0; i<nw; i++)
		{
			// read the names of variables. This could be used for error checking.
			// Or we could look for the variable number that we need, by matching the expected char[16] to the output here.
			file.read(temp,16);
			varnames.at(i)=temp;
		}
		// read physics module name, should be 'mhd'
		file.read(temp,16);
		std::string physics_module=temp;

		// read number of added parameters, should 1 for physics_module='mhd'
		// iterate through all and find gamma by the name of the variable
		file.read(reinterpret_cast<char*>(&neqpar), sizeof(int));
		double tmpeqpar;
		double gamma=0;
		std::vector<double> eqpar(neqpar);
		std::vector<std::string> eqparnames(neqpar);
		std::string gammaname="gamma"; //removed spaces by Krishna for compatibility with versions 3 & 5
		for (int i=0; i<neqpar; i++)
		{
			file.read(reinterpret_cast<char*>(&tmpeqpar), sizeof(double));
			eqpar.at(i)=tmpeqpar;
		}
		for (int i=0; i<neqpar; i++)
		{
			file.read(temp,16);
			eqparnames.at(i)=temp;
			// if the variable name matches gammaname, then take that value for gamma
			// this has as consequence that the last of such values will eventually be used
			if (eqparnames.at(i).substr(0,5).compare(gammaname)==0) gamma=eqpar.at(i); //added substr by Krishna
			std::cout << "Gamma:" << gamma << std::endl; //added by Krishna
		}
		if (gamma==0.)
		{
			std::cerr << "No value for gamma found" << std::endl;
			exit(EXIT_FAILURE);
		}

		// read the tree structure (previously called the "forest")
		file.seekg(treeoffset,ios_base::beg);
		std::vector<bool> forest=read_forest(file,nparents+nleafs);

		// read the block structure
		file.seekg(blockoffset,ios_base::beg);
		vector<vector<double>> leafs=read_blocks(file,ndim,nw,nglev1,nleafs,version);
		// we are now done with reading

		//compute the number of blocks in the simulation
		std::vector<int> nblocks(ndim);
		//compute the dx in each direction
		std::vector<double> cellsize(ndim);
		blocksandsize(nxlone,nx,xprobmin,xprobmax,nblocks,cellsize);
		for (int i=0; i<ndim; i++) std::cout << "cell size ["<< i << "]:" << cellsize[i] << std::endl; //added by Krishna

		// make sure the code is also working for 2D data!
		if (ndim<3) nblocks.push_back(1);
		std::vector<std::vector<int>> block_info=build_block_info_morton(nblocks,forest,ndim,nleafs);

		// Initialize the FoMo object, and load AMRVAC data into it
		Object=fomo_from_amrvac(ndim,nleafs,nglev1,nw,xprobmin,nx,cellsize,block_info,leafs,n_unit,Teunit,L_unit,gamma);
	}
	else
	{
		std::cerr << "could not open " << datfile;
	}

	return Object;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Added by Vaibhav Pant for AMRVAC 2 dat files (version 4)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FoMo::FoMoObject read_amrvac_dat_file_v4(const char* datfile, int version_number, const double n_unit, const double Teunit, const double L_unit)
{
	// Initialize the FoMo object
	FoMo::FoMoObject Object;

	// Open the file in the argument as ifstream
	// it is a binary file
	ifstream filehandle(datfile,ios::in|ios::binary);

	if (filehandle.is_open())
	{
		// read the filehandle into a buffer, this saves around 10s for each 2GB file
		stringstream file;
		file << filehandle.rdbuf();
		filehandle.close();

		int nleafs, levmax, ndim, ndir, nw, neqpar, it,stop;
		int version, treeoffset, blockoffset, nparents;
		double t;

		// follow specification at http://amrvac.org/md_doc_snapshot_format.html
		file.read(reinterpret_cast<char*>(&version), sizeof(int));
		file.read(reinterpret_cast<char*>(&treeoffset), sizeof(int));
		file.read(reinterpret_cast<char*>(&blockoffset), sizeof(int));
		file.read(reinterpret_cast<char*>(&nw), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndir), sizeof(int));
		file.read(reinterpret_cast<char*>(&ndim), sizeof(int));
		file.read(reinterpret_cast<char*>(&levmax), sizeof(int));
		file.read(reinterpret_cast<char*>(&nleafs), sizeof(int));
		file.read(reinterpret_cast<char*>(&nparents), sizeof(int));
		cout << "The simulation contains " << ndim << "D data, with " << nw << " variables. It has " << nleafs << " top-level blocks (=leafs)." << endl;
		file.read(reinterpret_cast<char*>(&it), sizeof(int));
		file.read(reinterpret_cast<char*>(&t), sizeof(double));
		//file.read(reinterpret_cast<char*>(&gamma), sizeof(double));
		std::cout<<"time= "<< t<< std::endl;
		vector<int> nxlone(ndim),nx(ndim);
		vector<double> xprobmin(ndim), xprobmax(ndim);
		// block_nx is earlier nxblock
		// domain_nx is earlier nxlone
		// amrvac.par is not needed anymore!
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&xprobmin.at(i)), sizeof(double));
			std::cout << "xprobmin(" << i << "):" << xprobmin.at(i) << std::endl;
		}
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&xprobmax.at(i)), sizeof(double));
			std::cout << "xprobmax(" << i << "):" << xprobmax.at(i) << std::endl;
		}
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&nxlone.at(i)), sizeof(int));
			std::cout << "domain_nx(" << i << "):" << nxlone.at(i) << std::endl;
		}
		int nglev1=1;
		for (int i=0; i<ndim; i++)
		{
			file.read(reinterpret_cast<char*>(&nx.at(i)), sizeof(int));
			nglev1*=nx.at(i);
			std::cout << "block_nx(" << i << "):" << nx.at(i) << std::endl;
		}
		std::cout << "nglev= " << nglev1 << std::endl;
		char temp[16];
		std::vector<std::string> varnames(nw);
		for (int i=0; i<nw; i++)
		{
			// read the names of variables. This could be used for error checking.
			// Or we could look for the variable number that we need, by matching the expected char[16] to the output here.
			file.read(temp,16);
			std::string vname=temp;
			varnames.at(i)=vname.substr(0,3);
			std::cout<<"varnames = "<<varnames.at(i)<<std::endl;
		}
		//std::cout<<nleafs<<std::endl;
		// read physics module name, should be 'mhd'
		file.read(temp,16);
		std::string physics_module=temp;
		physics_module = physics_module.substr(0,3);
		std::cout<<"physica module = "<<physics_module<<std::endl;

		// read number of added parameters, should 1 for physics_module='mhd'
		// iterate through all and find gamma by the name of the variable
		file.read(reinterpret_cast<char*>(&neqpar), sizeof(int));
		std::cout<<"neqpar ="<<neqpar<<std::endl;
		double tmpeqpar;
		double gamma=0;
		std::vector<double> eqpar(neqpar);
		std::vector<std::string> eqparnames(neqpar);
		std::string gammaname="gamma";
		for (int i=0; i<neqpar; i++)
		{
			file.read(reinterpret_cast<char*>(&tmpeqpar), sizeof(double));
			eqpar.at(i)=tmpeqpar;
			std::cout<<"tmpeqpar ="<<tmpeqpar<<std::endl;
		}
		for (int i=0; i<neqpar; i++)
		{
			file.read(temp,16);
			std::string name=temp;
			eqparnames.at(i)=name.substr(0,5);
			//std::cout<<"eqparnames ="<<eqparnames.at(i)<<"done."<<std::endl;
			// if the variable name matches gammaname, then take that value for gamma
			// this has as consequence that the last of such values will eventually be used
			if (eqparnames.at(i).compare(gammaname)==0) gamma=eqpar.at(i);
			//gamma=eqpar.at(i);
			std::cout<<"gamma is "<<gamma<<std::endl;
		}
		if (gamma==0.)
		{
			std::cerr << "No value for gamma found" << std::endl;
			exit(EXIT_FAILURE);
		}

		// read the tree structure (previously called the "forest")
		file.seekg(treeoffset,ios_base::beg);
		std::vector<bool> forest=read_forest(file,nparents+nleafs);

		// read the block structure
		file.seekg(blockoffset,ios_base::beg);
		vector<vector<double>> leafs=read_blocks(file,ndim,nw,nglev1,nleafs,version);
		//for (int i=0; 7; i++)
		//{
		//	std::cout<<leafs.at(0).at(0)<<std::endl;
		//	std::cout<<leafs.at(0).at(1)<<std::endl;
		//	std::cout<<leafs.at(0).at(2)<<std::endl;
		//	std::cout<<leafs.at(0).at(3)<<std::endl;
		//	std::cout<<leafs.at(0).at(4)<<std::endl;
		//	std::cout<<leafs.at(0).at(5)<<std::endl;
		//	std::cout<<leafs.at(0).at(6)<<std::endl;
		//	std::cout<<leafs.at(0).at(7)<<std::endl;
		//}
		//stop=0;
		//if (stop==0.)
		//{
		//	std::cerr << "Stop is encountered" << std::endl;
		//	exit(EXIT_FAILURE);
		//}
		// we are now done with reading

		//compute the number of blocks in the simulation
		std::vector<int> nblocks(ndim);
		//compute the dx in each direction
		std::vector<double> cellsize(ndim);
		blocksandsize(nxlone,nx,xprobmin,xprobmax,nblocks,cellsize);

		// make sure the code is also working for 2D data!
		if (ndim<3) nblocks.push_back(1);
		std::vector<std::vector<int>> block_info=build_block_info_morton(nblocks,forest,ndim,nleafs);

		// Initialize the FoMo object, and load AMRVAC data into it
		Object=fomo_from_amrvac(ndim,nleafs,nglev1,nw,xprobmin,nx,cellsize,block_info,leafs,n_unit,Teunit,L_unit,gamma);
		//std::cout<<Object.at(0).at(0)<<std::endl;
	}
	else
	{
		std::cerr << "could not open " << datfile;
	}

	return Object;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FoMo::FoMoObject read_amrvac_dat_file(const char* datfile, const char* amrvacpar, string amrvac_version, const int gamma_eqparposition, const double n_unit, const double Teunit, const double L_unit)
{
	int version_number=read_dat_version(datfile);

	FoMo::FoMoObject Object;

	if(std::find(compatible_versions.begin(), compatible_versions.end(), version_number) != compatible_versions.end()) {
		/* This is found in the compatible versions, and thus we use the new snapshot reader */
		std::cout << "Using new datfile reader." << std::endl;
		Object=read_amrvac_new_dat_file(datfile,version_number,n_unit,Teunit,L_unit);
	}
	//////////// Added by vaibhav pant for DATfiles version 4///////////////////////////////////////////////////
	else if(version_number == 4) {   //modified (to "else if") by Krishna
		/* We use the new snapshot reader, which is similar to version 3 except few changes in the definition of Gamma */
		std::cout << "Using new datfile reader (version 4) for AMRVAC 2.0" << std::endl;
		Object=read_amrvac_dat_file_v4(datfile,version_number,n_unit,Teunit,L_unit);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else {
		/* Not a compatible version, use the old reader */
		std::cout << "Using old datfile reader." << std::endl;
		Object=read_amrvac_old_dat_file(datfile,amrvacpar,amrvac_version,gamma_eqparposition,n_unit,Teunit,L_unit);
	}

	return Object;
}
