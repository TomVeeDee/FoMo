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

void FoMo::GoftCube::writegoftcube(const std::string filename)
{
	// write out goftcube to file "filename"
	int commrank;
	std::string space=" ";
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
			out << dim << std::endl;
			out << ng << std::endl;
			out << nvars << std::endl;
			out << this->readchiantifile() << std::endl;
			out << this->readabundfile() << std::endl;
			for (int j=0; j<ng; j++)
			{
				for (int i=0; i<dim; i++)
				{ 
					out << grid[i][j] << space;
				}
				for (int i=0; i<nvars; i++)
				{
					out << vars[i][j] << space;
				}
				out << std::endl;
			}
		}
		else std::cerr << "Unable to write to " << filename << std::endl;
		out.close();
	}
}

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