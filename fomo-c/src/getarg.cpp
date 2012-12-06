#include "header.h"
#include <getopt.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>

// Initialize global variables
int reuse = 0, png = 0, mpeg = 0, array = 0;
double length = 200, width = 5, beta = 0, magfield = 15, rhoint = 1.4, /*rhoext = 0.4,*/ contrast = 3, thickness = 2, alpha = 0.78, ampl = 0.2, phase=0.;
double psi = M_PI/3., l=M_PI/6., b=M_PI/3.;
double nperiods = 1., nframes = 5.;

void printusage(const char* programname){
// print the usage of the program
	printf("Usage: %s [--reuse] [--help] [--parameter value] \n\n",programname);
	printf("  --reuse       Reuse the previous results. Setting the physical parameters does not have any effect!\n");
	printf("  --help        Print this message\n\n");
	printf("Physical parameters:\n");
	printf("  --length      Length (in Mm)\n");
	printf("  --width       Width (in Mm)\n");
	printf("  --beta        Plasma-beta\n");
	printf("  --contrast    Density contrast = exterior density / interior density\n");
	printf("  --thickness   Thickness of the inhomogeneous layer\n");
	printf("  --alpha       Alpha (between 0 and 1)\n\n");
	printf("Rescaling parameters:\n");
	printf("  --magfield    Magnetic field strength (in Gauss) (still unused)\n");
	printf("  --rhoint      Interior density (in 10^9 cm-3) (still unused)\n");
	printf("  --ampl        Amplitude of oscillation (in Mm) (still unused)\n");
	printf("  --phase       Phase of the wave (exp(I*m*phi+phase*I))\n\n");
	printf("Orientation parameters:\n");
	printf("  --psi         Angle of the plane of the loop with the circle parallel to the equator (in radians)\n");
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
	printf("  --nP          Number of periods observed\n");
	printf("  --nf          Number of frames per period\n");
	printf("\n");
}

void getarg(int argc, char* argv[])
{
// Get arguments and set parameters
// If reuse is different from 0, everything else is unused
	char* optstr="L:w:b:c:t:a:B:i:A:z:p:l:q:GMyP:f:r?";
	struct option longopts[]={
// name, require argument, flag, value
		{"length",1,0,'L'},
		{"width",1,0,'w'},
		{"beta",1,0,'b'},
		{"contrast",1,0,'c'},
		{"thickness",1,0,'t'},
		{"alpha",1,0,'a'},
		{"magfield",1,0,'B'},
		{"rhoint",1,0,'i'},
		{"ampl",1,0,'A'},
		{"phase",1,0,'z'},
		{"psi",1,0,'p'},
		{"l",1,0,'l'},
		{"b",1,0,'q'},
// mark the q above, not an evident abbreviation
		{"png",0,0,'G'},
		{"mpeg",0,0,'M'},
		{"array",0,0,'y'},
		{"nP",1,0,'P'},
		{"nf",1,0,'f'},
		{"reuse",0,0,'r'},
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
			// physical parameters
			case 'L':
				if (reuse!=1) length = atof(optarg);
				break;
			case 'w':
				if (reuse!=1) width = atof(optarg);
				break;
			case 'b':
				if (reuse!=1) beta = atof(optarg);
				break;
			case 'c':
				if (reuse!=1) contrast = atof(optarg);
				break;
			case 't':
				if (reuse!=1) thickness = atof(optarg);
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
			case 'A':
				ampl = atof(optarg);
				break;
			case 'z':
				phase = atof(optarg);
				break;
			// orientation parameters
			case 'p':
				psi = atof(optarg);
				break;
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
				array = 1;
				break;
			case 'P':
				nperiods = atof(optarg);
				break;
			case 'f':
				nframes = atof(optarg);
				break;
			case 'r':
				reuse = 1;
				printf("Warning: reusing old parameters\n");
				getfile();
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
	if ((!mpeg)&&(!png)) array = 1;
}
