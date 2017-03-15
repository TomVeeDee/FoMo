#include "FoMo.h"
#include "FoMo-internal.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

/*#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_triangulation_3;

void FoMo::GoftCube::writegoftcube(const std::string filename, const Delaunay_triangulation_3 *DTpointer)
{
	// write out goftcube to file "filename"
	int commrank;
	string space=" ";
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		ofstream out(filename,ios::binary|ios::ate);
		if (out.is_open())
		{
			int nvars = this->readnvars();
			int dim = this->readdim();
			int ng = this->readngrid();
			out << dim << std::endl;
			out << ng << std::endl;
			out << nvars << std::endl;
			out << this->readchiantifile() << std::endl;
			out << this->readabundfile() << std::endl;
			FoMo::tgrid grid=this->readgrid();
			FoMo::tvars allvar;
			allvar.resize(nvars);
			for (int i=0; i<nvars; i++)
			{
				allvar[i]=this->readvar(i);
			}
			for (int j=0; j<ng; j++)
			{
				for (int i=0; i<dim; i++)
				{ 
					out << grid[i][j] << space;
				}
				for (int i=0; i<nvars; i++)
				{
					out << allvar[i][j] << space;
				}
				out << endl;
			}
			if (DTpointer!=NULL) out << *DTpointer;
		}
		else cout << "Unable to write to " << filename << endl;
		out.close();
	}
}

void FoMo::GoftCube::readgoftcube(const std::string emissionsave, Delaunay_triangulation_3 * DTpointer)
{
	// Read GoftCube datacube from file "emissionsave"
	// Otherwise read in data cube (simulation snapshots) from directory specified in main function

	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		ifstream in(emissionsave,ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> ng;
			in >> nvars;
			in >> chiantifile;
			in >> abundfile;
			
			grid.resize(dim);
			vars.resize(nvars);
			lambda0=readgoftfromchianti(chiantifile);
			
			for (int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
			}
			for (int i=0; i<nvars; i++)
			{
				vars[i].resize(ng);
			}
			for (int j=0; j<ng; j++)
			{
				for (int i=0; i<dim; i++) in >> grid[i][j];
				for (int i=0; i<nvars; i++) in >> allvar[i][j];
			}

			if (!in.eof() && (DTpointer!=NULL)) 
			{
				in >> *DTpointer;
				// Check is this is a valid Delaunay triangulation
				// assert(DTpointer->is_valid());
			}
		}
		else cout << "Unable to read " << emissionsave << endl;
		in.close();
	}
}*/

/**
 * @brief This writes out the contents of the GoftCube.
 * 
 * This member function writes out the contents of the GoftCube to the hard disk, in the filename given by the argument. 
 * The file can then be read into IDL with the commands provide under the idl/ directory, here in particular readgoftcube.pro.
 * Normally, one would use it to post-process the forward modelling results (when using with a RenderCube), but it is also very useful for debugging purposes when using it with a GoftCube before the rendering.
 * @param filename This parameter specifies which filename the data needs to be written to.
 */
void FoMo::GoftCube::writegoftcube(const std::string filename)
{
	// write out goftcube to file "filename"
	int commrank;
	std::string space=" ";
	size_t chiantisize=chiantifile.size();
	size_t abundsize=abundfile.size();
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		std::ofstream out(filename,std::ios::binary|std::ios::ate);
		if (out.is_open())
		{
			int nvars = this->readnvars();
			int dim = this->readdim();
			int ng = this->readngrid();
			out.write(reinterpret_cast<const char*>(&dim),sizeof(dim));
			out.write(reinterpret_cast<const char*>(&ng),sizeof(ng));
			out.write(reinterpret_cast<const char*>(&nvars),sizeof(nvars));
			out.write(reinterpret_cast<const char*>(&chiantisize),sizeof(chiantisize));
			out.write(chiantifile.c_str(),chiantisize);
			out.write(reinterpret_cast<const char*>(&abundsize),sizeof(abundsize));
			out.write(abundfile.c_str(),abundsize);
			
			for (int j=0; j<ng; j++)
			{
				for (int i=0; i<dim; i++)
				{ 
					out.write(reinterpret_cast<const char*>(&grid[i][j]),sizeof(grid[i][j]));
				}
				for (int i=0; i<nvars; i++)
				{
					out.write(reinterpret_cast<const char*>(&vars[i][j]),sizeof(vars[i][j]));
				}
			}
		}
		else std::cerr << "Unable to write to " << filename << std::endl;
		out.close();
	}
}

/**
 * @brief This is a routine for reading a previously saved GoftCube.
 * 
 * This member reads the information from file emissionsave, which was previously written by GoftCube::writegoftcube, and stores it in the GoftCube.
 * @param emissionsave The filename from which the GoftCube should be read.
 */
void FoMo::GoftCube::readgoftcube(const std::string emissionsave)
{
	// Read GoftCube datacube from file "emissionsave"
	// Otherwise read in data cube (simulation snapshots) from directory specified in main function

	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		std::ifstream in(emissionsave,std::ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> ng;
			in >> nvars;
			in >> chiantifile;
			in >> abundfile;
			
			grid.resize(dim);
			vars.resize(nvars);
			lambda0=readgoftfromchianti(chiantifile);
			
			for (unsigned int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
			}
			for (unsigned int i=0; i<nvars; i++)
			{
				vars[i].resize(ng);
			}
			for (unsigned int j=0; j<ng; j++)
			{
				for (unsigned int i=0; i<dim; i++) in >> grid[i][j];
				for (unsigned int i=0; i<nvars; i++) in >> vars[i][j];
			}
		}
		else std::cerr << "Unable to read " << emissionsave << std::endl;
		in.close();
	}
}