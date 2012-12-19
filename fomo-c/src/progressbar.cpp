#include <iostream>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

const int nstars=56;
const double nperc=8.;
const char* point=".";

void progressbar(const int i, const int mini, const int maxi)
{
	int range=maxi-mini+1;
	int threadrank=0;
	int nthreads=1;
#ifdef _OPENMP
	nthreads=omp_get_num_threads();
	range/=nthreads;
	threadrank=omp_get_thread_num();
#endif
	for (int l=1; l<nperc; l++)
		if (((i-mini-1)<l/nperc*range)&&(i-mini>=l/nperc*range)&&(threadrank==0))
		{
			cout << floor(l/nperc*1000+0.5)/10 << "%";
			cout.flush();
		}
	if (((i-mini)%(range/nstars)==0)&&(threadrank==0)) 
	{
		cout << point;
		// need a flush, otherwise everything is kept until the next endl
		cout.flush();
	}
}
