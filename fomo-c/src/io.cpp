#include "header.h"
#include <fstream>
#include <cstdlib>

char* configfile="fomo-c.conf";

void getfile()
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
}

void writeemissioncube(const cube goftcube)
{
	// write out goftcube to file "emissionsave"
	int commrank;
	char* space=" ";
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		ofstream out(emissionsave,ios::binary|ios::ate);
		if (out.is_open())
		{
			int nvars = goftcube.readnvars();
			int dim = goftcube.readdim();
			int intqtype=int(goftcube.readtype());
			int ng = goftcube.readngrid();
			out << dim << space;
			out << intqtype << space;
			out << ng << space;
			out << nvars << space;
			cout << dim << " " << intqtype << " " << ng << " " << nvars;
			tgrid grid=goftcube.readgrid();
			for (int i=0; i<dim; i++)
				for (int j=0; j<ng; j++)
				{ 
					out << grid[i][j] << space;
				}
			for (int i=0; i<nvars; i++)
			{
				tphysvar var=goftcube.readvar(i);
				for (int j=0; j<ng; j++)
				{
					out << var[j] << space;
				}
			}
		}
		else cout << "Unable to write to " << emissionsave << endl;
		out.close();
	}
}

cube reademissioncube()
{
	// read emission datacube from file "emissionsave" when the reuse parameter is switched on
	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		int dim, ng, nvars, intqtype;

		ifstream in(emissionsave,ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> intqtype;
			in >> ng;
			in >> nvars;
			cout << "test" << dim << " " << intqtype << " " << ng << " " << nvars;
			cube goftcube(nvars,ng,dim);
			goftcube.settype(EqType(intqtype));
			return goftcube;
			tgrid grid= new tcoord[dim];
			for (int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
				for (int j=0; j<ng; j++)
				{
					in >> grid[i][j];
				}
			}
			goftcube.setgrid(grid);
			tphysvar var;
			var.resize(ng);
			for (int i=0; i<nvars; i++)
			{
				for (int j=0; j<ng; j++)
				{
					in >> var[j];
				}
				goftcube.setvar(i,var);
			}
			return goftcube;
		}
		else cout << "Unable to read " << emissionsave << endl;
		in.close();
	}
	cube badcube(1,1,1);
	return badcube;
}
