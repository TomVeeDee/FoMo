#include "header.h"
#include <getopt.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>

// Initialize global variables
int reuse = 0, png = 0, mpeg = 0, warray = 0;
double length = 200, width = 5, magfield = 15, rhoint = 1.4, contrast = 3, thickness = 2, alpha = 0.78, ampl = 0.2, phase=0.;
double l=M_PI/6., b=M_PI/3.;
string chiantifile="chiantitables/goft_table_fe_12_0194.dat";
string abundfile="/empty";
string emissionsave="fomo-c.emissionsave";

void printusage(const char* programname){
// print the usage of the program
	printf("Usage: %s [--reuse] [--help] [--parameter value] \n\n",programname);
	printf("  --reuse       Reuse the previous results. Setting the physical parameters does not have any effect!\n");
	printf("  --emissionsave    When running without --reuse, this is the filename where the G(T) interpolated datacube will be stored (filename 'none' skips this step). When running with --reuse, this is the file that will be read in\n");
	printf("  --help        Print this message\n\n");
	printf("CHIANTI parameters:\n");
	printf("  --goftfile    Filename containing the table with the G(T) written by CHIANTI\n");
	printf("  --abundfile   Filename containing the table with the abundances (standard sun_coronal.abund is used)\n\n");
	printf("Equilibrium parameters:\n");
	printf("  --length      Length (in Mm)\n");
	printf("  --width       Width (in Mm)\n");
	printf("  --contrast    Density contrast = exterior density / interior density\n");
	printf("  --thickness   Thickness of the inhomogeneous layer (l/R)\n");
	printf("  --alpha       Alpha (between 0 and 1)\n\n");
	printf("Rescaling parameters:\n");
	printf("  --magfield    Magnetic field strength (in Gauss) (still unused)\n");
	printf("  --rhoint      Interior density (in 10^9 cm-3) (still unused)\n");
	printf("Orientation parameters:\n");
	printf("  --l           Longitude of the symmetry centre of the loop (in radians)\n");
	printf("  --b           Lattitude of the symmetry centre of the loop (in radians)\n\n");
	printf("Output parameters:\n");
#ifdef HAVEPNG
	printf("  --png		Enable writing of png's\n");
#endif
#ifdef HAVEMPEG
	printf("  --mpeg	Enable writing of movies\n");
#endif
	printf("  --array	Enable writing of array\n");
	printf("\n");
}

void getarg(int argc, char* argv[])
{
// Get arguments and set parameters
// If reuse is different from 0, everything else is unused
	char* optstr="C:A:L:w:t:c:a:B:i:l:q:GMy:rs:?";
	struct option longopts[]={
// name, require argument, flag, value
		{"goftfile",1,0,'C'},
		{"abundfile",1,0,'A'},
		{"length",1,0,'L'},
		{"width",1,0,'w'},
		{"thickness",1,0,'t'},
		{"contrast",1,0,'c'},
		{"alpha",1,0,'a'},
		{"magfield",1,0,'B'},
		{"rhoint",1,0,'i'},
		{"l",1,0,'l'},
		{"b",1,0,'q'},
// mark the q above, not an evident abbreviation
		{"png",0,0,'G'},
		{"mpeg",0,0,'M'},
		{"array",0,0,'y'},
		{"reuse",0,0,'r'},
		{"emissionsave",1,0,'s'},
		{"help",0,0,'h'},
		{NULL,0,NULL,0}
		};
	int longindex = 0;
	if (argc == 1) {
		cout << "Error: no arguments given\n\n";
		printusage(argv[0]);
		exit(EXIT_FAILURE);
	}
	int opt = getopt_long(argc, argv, optstr, longopts, &longindex);
	while (opt != -1){
		switch (opt){
			// CHIANTI parameters
			case 'C':
				chiantifile=optarg;
				break;
			case 'A':
				abundfile=optarg;
				break;
			// physical parameters
			case 'L':
				if (reuse!=1) length = atof(optarg);
				break;
			case 'w':
				if (reuse!=1) width = atof(optarg);
				break;
			case 't':
				if (reuse!=1) thickness = atof(optarg);
				break;
			case 'c':
				if (reuse!=1) contrast = atof(optarg);
				break;
			case 'a':
				if (reuse!=1) alpha = atof(optarg);
				break;
			// rescaling parameters
			case 'B':
				magfield = atof(optarg);
				break;
			case 'i':
				rhoint = atof(optarg);
				break;
			// orientation parameters
			case 'l':
				l = atof(optarg);
				break;
			case 'q':
				b = atof(optarg);
				break;
			// output parameters
			case 'G':
#ifdef HAVEPNG
				png = 1;
#endif
				break;
			case 'M':
#ifdef HAVEMPEG
				mpeg = 1;
#endif
				break;
			case 'y':
				warray = 1;
				break;
			case 'r':
				reuse = 1;
				printf("Warning: reusing old parameters\n");
				getfile();
				break;
			case 's':
				emissionsave = optarg;
				break;
			case 'h':
				printusage(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			default:
				printf("Error: unknown option\n\n");
				printusage(argv[0]);
				exit(EXIT_FAILURE);
				break;
		}
		opt = getopt_long(argc, argv, optstr, longopts, &longindex);
	}
	if ((!mpeg)&&(!png)) warray = 1;
}
