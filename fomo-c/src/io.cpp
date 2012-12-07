#include "header.h"
#include <fstream>
#include <cstdlib>

int nframes=1;

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

const int filterlength = 2149;

void readfilters(const int i, double array[][2])
// i is the selected filter
// i=1: 195/171
// i=2: 284/195
// i=3: 171
// i=4: 195
// i=5: 284
{
	char* file = "filter.dat";
	ifstream in(file);
	if (!in) {
		cout << "No " << file << ": unable to determine filtervalues of TRACE";
		exit(EXIT_FAILURE);
	}
	double temp[6];
	int j=0;
	in >> temp[0] >> temp[1] >> temp[2] >> temp[3] >> temp[4] >> temp[5];
	while (!in.eof() )
	{
		array[j][0]=temp[0];
		array[j][1]=temp[i];
		in >> temp[0] >> temp[1] >> temp[2] >> temp[3] >> temp[4] >> temp[5];
		j++;
	}
}
