#include "header.h"
#include <sstream>
#include <fstream>
#include <cmath>

void writearray(double * const * const image, const int t)
{
	int seqnum = t;
	stringstream ss;
	ss << "outputarray.t" << seqnum << ".txt";
	string graphic_out=ss.str();
	ofstream s(graphic_out.c_str());
	for (int j=0; j<y_pixel; j++)
		for (int k=0; k<x_pixel; k++)
			s << image[j][k] << endl;
}

