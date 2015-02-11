#include "header.h"
#include <fstream>
#include <cstdlib>

string configfile="fomo-c.conf";

/*void getfile()
{
	ifstream in(configfile);
	if (!in) {
		cout << "No " << configfile << ": please rerun the program without the --reuse option";
		exit(EXIT_FAILURE);
	}
	in >> length
	   >> width
	   >> magfield
	   >> rhoint
	   >> contrast
	   >> thickness
	   >> alpha
	   >> l
	   >> b;
	   
}
*/
void writeparameters(ostream& s, char v = '0')
{
	int commrank;
#ifdef HAVEMPI
	MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	switch (v) {
		case 'v':
			s << "length " << length << endl
			  << "width " << width << endl
			  << "magnetic field " << magfield << endl
			  << "interior density " << rhoint << endl
			  << "density contrast " << contrast << endl
			  << "thickness " << thickness << endl
			  << "alpha " << alpha << endl
			  << "longitude " << l << endl
			  << "lattitude " << b << endl;
			break;
		case '0':
			s << length << endl
			  << width << endl
			  << magfield << endl
			  << rhoint << endl
			  << contrast << endl
			  << thickness << endl
			  << alpha << endl
			  << l << endl
			  << b << endl;
			break;
		default:
			cout << "Warning: incorrect use of writeparameters";
			break;
	}
}
/*
void writefile()
{
	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		ofstream out(configfile);
		writeparameters(out);
	}
} */

void writeemissioncube(const cube goftcube, const string filename, const Delaunay_triangulation_3 *DTpointer)
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
			int nvars = goftcube.readnvars();
			int dim = goftcube.readdim();
			int intqtype=int(goftcube.readtype());
			int ng = goftcube.readngrid();
			out << dim << endl;
			out << intqtype << endl;
			out << ng << endl;
			out << nvars << endl;
			tgrid grid=goftcube.readgrid();
			vector<tphysvar> allvar;
			allvar.resize(nvars);
			for (int i=0; i<nvars; i++)
			{
				allvar[i]=goftcube.readvar(i);
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

void reademissioncube(cube &resultcube, const string emissionsave, Delaunay_triangulation_3 * DTpointer)
{
	// Read emission datacube from file "emissionsave" when the reuse parameter is switched on
	// Otherwise read in data cube (simulation snapshots) from directory specified in main function

	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		int dim, ng, nvars, intqtype;
		EqType resulttype=resultcube.readtype();

		ifstream in(emissionsave,ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> intqtype;
			in >> ng;
			in >> nvars;
			cube goftcube(nvars,ng,dim);
			if ((EqType(intqtype) <= empty) && (resulttype == empty))
			{
				goftcube.settype(EqType(intqtype));
			}
			else
			{
				goftcube.settype(resulttype);
			}
			tgrid grid= new tcoord[dim];
			vector<tphysvar> allvar;
			allvar.resize(nvars);
			for (int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
			}
			for (int i=0; i<nvars; i++)
			{
				allvar[i].resize(ng);
			}
			for (int j=0; j<ng; j++)
			{
				for (int i=0; i<dim; i++) in >> grid[i][j];
				for (int i=0; i<nvars; i++) in >> allvar[i][j];
			}
			goftcube.setgrid(grid);
			for (int i=0; i<nvars; i++)
			{
				goftcube.setvar(i,allvar[i]);
			}
			resultcube=goftcube;

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
}
/*
void reademissioncube_byt(cube &resultcube, const string emissionsave, Delaunay_triangulation_3 * DTpointer)
{
 // read binary data added by D.Y 10 Feb 2015
	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
                int dim, nvars, intqtype;
                long int ng;
		EqType resulttype=resultcube.readtype();

		ifstream infile(emissionsave,ios_base::binary | ios_base::in);
		if (infile.is_open())
		{
                        infile.read(&dim,sizeof(dim));
                        infile.read(&intqtype,sizeof(intqtype));
                        infile.read(&ng,sizeof(ng));
                        infile.read(&nvars,sizeof(nvars));
                        float array[dim+nvars][ng];
                        infile.read(&array,sizeof(array));
			cube goftcube(nvars,ng,dim);
			if ((EqType(intqtype) <= empty) && (resulttype == empty))
			{
				goftcube.settype(EqType(intqtype));
			}
			else
			{
				goftcube.settype(resulttype);
			}
			tgrid grid= new tcoord[dim];
			vector<tphysvar> allvar;
			allvar.resize(nvars);
			for (int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
			}
			for (int i=0; i<nvars; i++)
			{
				allvar[i].resize(ng);
			}
                     
			for (int j=0; j<ng; j++)
	        	{
				for (int i=0; i<dim; i++) grid[i][j]=array[i][j];
				for (int i=0; i<nvars; i++) allvar[i][j]==array[dim+i][j];
			}
                        
			goftcube.setgrid(grid);
			for (int i=0; i<nvars; i++)
			{
				goftcube.setvar(i,allvar[i]);
			}
			resultcube=goftcube;

			if (!infile.eof() && (DTpointer!=NULL)) 
			{
				infile >> *DTpointer;
				// Check is this is a valid Delaunay triangulation
				// assert(DTpointer->is_valid());
			}
		}
		else cout << "Unable to read " << emissionsave << endl;
		infile.close();
	}
}*/
