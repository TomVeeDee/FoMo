#include <iostream>
#include <cmath>

using namespace std;

const int nstars=56;
const double nperc=8.;
const char* point=".";

void progressbar(const int i, const int mini, const int maxi)
{
	for (int l=1; l<nperc; l++)
		if (((i-mini-1)<l/nperc*(maxi-mini+1))&&(i-mini>=l/nperc*(maxi-mini+1)))
		{
			cout << floor(l/nperc*1000+0.5)/10 << "%";
			cout.flush();
		}
	if ((i-mini)%((maxi-mini+1)/nstars)==0) 
	{
		cout << point;
		// need a flush, otherwise everything is kept until the next endl
		cout.flush();
	}
}
