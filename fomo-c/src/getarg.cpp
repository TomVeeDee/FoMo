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
// added by D.Y. 28 Oct 2014
string infileini="infileini";
string outfileini="outfileini";
int tstart=0, tend=0,tstep=1;


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
        printf("  --infileini    input data cube initials\n\n");
        printf("  --outfileini   output data cube initials\n\n");
        printf("  --tstart   the first frame to process, default=0\n\n");
        printf("  --tend     the final frame to process, default=1\n\n");
        printf("  --tstep     the time steps to process input files, default=1\n\n");
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
//	char* optstr="C:A:L:w:t:c:a:B:i:l:q:GMy:rs:f:u:T:E:P:?";
 char const *optstr="C:A:L:w:t:c:a:B:i:l:q:GMy:rs:f:u:T:E:P:?";// D.Y. 28 oct 2014
	struct option longopts[]={
// name, require argument, flag, value
		{"goftfile",     required_argument,0,'C'},
		{"abundfile",    required_argument,0,'A'},
		{"length",       required_argument,0,'L'},
		{"width",        required_argument,0,'w'},
		{"thickness",    required_argument,0,'t'},
		{"contrast",     required_argument,0,'c'},
		{"alpha",        required_argument,0,'a'},
		{"magfield",     required_argument,0,'B'},
		{"rhoint",       required_argument,0,'i'},
		{"l",            required_argument,0,'l'},
		{"b",            required_argument,0,'q'},
// mark the q above, not an evident abbreviation
		{"png",          no_argument,      0,'G'},
		{"mpeg",         no_argument,      0,'M'},
		{"array",        no_argument,      0,'y'},
		{"reuse",        no_argument,      0,'r'},
		{"emissionsave", required_argument,0,'s'},
                {"infileini",    required_argument,   0,'f'}, // D.Y. 28 Oct 2014
                {"outfileini",   required_argument,  0,'u'}, // D.Y. 28 Oct 2014
                {"tstart",       required_argument,  0,'T'},
                {"tend",         required_argument,  0,'E'},
                {"tstep",        required_argument,  0,'P'}, 
		{"help",         no_argument,      0,'h'},
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
                        case 'f':
                                infileini = optarg;
                                break;
                        case 'u':
                                outfileini = optarg;
                                 break;
                        case 'T':
                                 tstart = atoi(optarg);
                                  break;
                        case 'E':
                                 tend = atoi(optarg);
                                  break;
                        case 'P':
                                  tstep = atoi(optarg);
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
