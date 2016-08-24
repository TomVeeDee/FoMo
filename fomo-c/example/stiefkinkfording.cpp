#include "FoMo.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <time.h> 

using namespace std;

void printusage(const char* programname){
// print the usage of the program
	printf("Usage: %s [--help] [--parameter value] \n\n",programname);
	printf("  --help        Print this message\n\n");
	printf("CHIANTI parameters:\n");
	printf("  --goftfile    Filename containing the table with the G(T) written by CHIANTI\n");
	printf("  --abundfile   Filename containing the table with the abundances (standard sun_coronal.abund is used)\n\n");
	printf("IO parameters:\n");
        printf("  --infileini   input data cube initials\n");
        printf("  --outfileini  output data cube initials\n");
        printf("  --tstart      the first frame to process, default=0\n");
        printf("  --tend        the final frame to process, default=1\n");
        printf("  --tstep       the time steps to process input files, default=1\n");
	printf("\n");
}

string chiantifile="chiantitables/goft_table_fe_12_0194.dat";
string abundfile="/empty";
string emissionsave="fomo-c.emissionsave";
string infileini="infileini";
string outfileini="outfileini";
int tstart=0, tend=0,tstep=1;

void getarg(int argc, char* argv[])
{
// Get arguments and set parameters
// If reuse is different from 0, everything else is unused
//	char* optstr="C:A:L:w:t:c:a:B:i:l:q:GMy:rs:f:u:T:E:P:?";
 char const *optstr="C:A:f:u:T:E:P:?";// D.Y. 28 oct 2014
	struct option longopts[]={
// name, require argument, flag, value
		{"goftfile",     required_argument,0,'C'},
		{"abundfile",    required_argument,0,'A'},
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
		std::cout << "Error: no arguments given\n\n";
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
                        // added by DY 31 Oct 2014
                        case 'f':    
                                infileini = optarg;
                                break;
                        case 'u':
                                outfileini = optarg;
                                 break;
                        case 'T':
                                 tstart = atoi(optarg);
                                 if (tend == 0) {tend=tstart;}  
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
}

int main(int argc, char* argv[])
{
        clock_t t0;
	getarg(argc,argv);  // Read in arguments from call to main function. For an overview of different input parameters, run program with --help argument.

	double globalmax, globalmin;
	globalmax = 0.; globalmin = 0.;
	// Create G(T) interpolated cube or artificial images (depending on --reuse option) 
        //const int nframes=tend-tstart+1; //Number of time steps = number of simulation snapshots
	// double pi=4*atan(1.);
	vector<double> langles={M_PI/6.}; // M_PI/2.,M_PI/3.,M_PI/4. Rotation angles around y-axis
	vector<double> bangles={M_PI/6.};
	stringstream ss;
       // if not dynamicDT the triangulation can be done only once, however, it is saved for every snapshot for restart!
        bool dynamicDT=false;
        bool doDT=true;
        bool readDT=true;
	for (int t=tstart; t<=tend; t=t+tstep)
	{
		cout << endl << "Doing timestep " << t << endl << flush;
	        ss << infileini; // argument input
         	ss << setfill('0') << setw(3) << t;
		ss << ".dat";
		string filename=ss.str(); // input file
		
		// create a FoMoObject of dimension 3
		FoMo::FoMoObject Object(3);
		Object.setchiantifile(chiantifile);
		Object.setabundfile(abundfile);
		Object.setrendermethod("CGAL"); // CGAL is the default rendermethod
		Object.setobservationtype(FoMo::Spectroscopic);
		// adjust the resolution with these parameters
		int x_pixel=101;
		int y_pixel=102;
		int z_pixel=300;
		int lambda_pixel=30;
		double lambda_width=200000;
		Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);

		if (t==tstart)
                {
                	doDT=true;
                  	readDT=true;
                }
                
                t0=clock();
            	cout << "Reading in snapshot " << filename << "... " << flush;
            	int dim, ng, nvars, intqtype;
		ifstream in(filename,ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> intqtype;
			in >> ng;
			in >> nvars;
			double tmpvar;
			while ( !in.eof() )
			{
				vector<double> coordinates;
				// The points should be three dimensional, and thus we iterate until 3.
				for (unsigned int i=0; i<3; i++)
				{
					// Read one value of a coordinate at a time.
					in >> tmpvar;
					// push it back into the coordinate vector.
					coordinates.push_back(tmpvar);
				}
				vector<double> variables;
				// There should be 5 variables (n, T, vx, vy, vz), and iterate until 5.
				for (unsigned int i=0; i<5; i++)
				{
					// Read one value of a variable at a time.
					in >> tmpvar;
					if (i > 1) tmpvar*= 1e6; // convert speeds to m/s
					// push it back into the variable vector.
					variables.push_back(tmpvar);
				}
				
				// If the file has not reached the end, this is probably a valid data point. Push it back 
				// into the FoMoObject.
				if ( !in.eof() ) Object.push_back_datapoint(coordinates, variables);
			}
		}
		else cout << "Unable to read " << filename << endl;
		in.close();
		cout << "Done!" << endl << flush;
                cout<<"Reading time used "<<(float(clock()-t0)/CLOCKS_PER_SEC)<<"[cpu second]"<<endl;
                t0=clock();
		
		// write out observ cube
		ss.str("");
		// ss << "/users/cpa/dyuan/fomodata/testemislos";
                ss << outfileini; // output file initials
		ss << "t";
		ss << setfill('0') << setw(3) << t;
		string outfile = ss.str();
		ss.str("");
		Object.render(langles,bangles);
	}			// end t
	std::cout << "Gelukt\n" ;
	return EXIT_SUCCESS;
}
