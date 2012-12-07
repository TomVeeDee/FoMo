#include "header.h"
#include <cmath>

void writebar(double * const * const image, const int pos, const int digit)
{
	// char color = findmax(image,i,j);
	int i,j;
	double color = findmin(image,&i,&j);
	int spacing = 22;
	switch (digit){
		case 1:
			for (int y=y_pixel-42; y<y_pixel-39; y++)
			for (int x=x_pixel-24-(pos-1)*spacing; x<x_pixel-12-(pos-1)*spacing; x++)
			{
				image[y][x] = color;
			}
			break;
		case 2:
			for (int y=y_pixel-39 ; y<y_pixel-27; y++)
			for (int x=x_pixel-27-(pos-1)*spacing ; x<x_pixel-24-(pos-1)*spacing; x++)
				image[y][x] = color;
			break;
		case 3:
			for (int y=y_pixel-39 ; y<y_pixel-27; y++)
			for (int x=x_pixel-12-(pos-1)*spacing ; x<x_pixel-9-(pos-1)*spacing; x++)
				image[y][x] = color;
			break;
		case 4:
			for (int y=y_pixel-27 ; y<y_pixel-24; y++)
			for (int x=x_pixel-24-(pos-1)*spacing ; x<x_pixel-12-(pos-1)*spacing; x++)
				image[y][x] = color;
			break;
		case 5:
			for (int y=y_pixel-24 ; y<y_pixel-12; y++)
			for (int x=x_pixel-27-(pos-1)*spacing ; x<x_pixel-24-(pos-1)*spacing; x++)
				image[y][x] = color;
			break;
		case 6:
			for (int y=y_pixel-24 ; y<y_pixel-12; y++)
			for (int x=x_pixel-12-(pos-1)*spacing ; x<x_pixel-9-(pos-1)*spacing; x++)
				image[y][x] = color;
			break;
		case 7:
			for (int y=y_pixel-12 ; y<y_pixel-9; y++)
			for (int x=x_pixel-24-(pos-1)*spacing ; x<x_pixel-12-(pos-1)*spacing; x++)
			{
				image[y][x] = color;
			}
			break;
		default:
			break;
	}
}

void writedigit(double * const * const image, const int pos, const int i)
{
	switch (i){
	case 0:
		writebar(image,pos,1);
		writebar(image,pos,2);
		writebar(image,pos,3);
		writebar(image,pos,5);
		writebar(image,pos,6);
		writebar(image,pos,7);
		break;
	case 1:
		writebar(image,pos,3);
		writebar(image,pos,6);
		break;
	case 2:
		writebar(image,pos,1);
		writebar(image,pos,3);
		writebar(image,pos,4);
		writebar(image,pos,5);
		writebar(image,pos,7);
		break;
	case 3: 
		writebar(image,pos,1);
		writebar(image,pos,3);
		writebar(image,pos,4);
		writebar(image,pos,6);
		writebar(image,pos,7);
		break;
	case 4:
		writebar(image,pos,2);
		writebar(image,pos,3);
		writebar(image,pos,4);
		writebar(image,pos,6);
		break;
	case 5:
		writebar(image,pos,1);
		writebar(image,pos,2);
		writebar(image,pos,4);
		writebar(image,pos,6);
		writebar(image,pos,7);
		break;
	case 6:
		writebar(image,pos,2);
		writebar(image,pos,4);
		writebar(image,pos,5);
		writebar(image,pos,6);
		writebar(image,pos,7);
		break;
	case 7: 
		writebar(image,pos,1);
		writebar(image,pos,3);
		writebar(image,pos,6);
		break;
	case 8:
		writebar(image,pos,1);
		writebar(image,pos,2);
		writebar(image,pos,3);
		writebar(image,pos,4);
		writebar(image,pos,5);
		writebar(image,pos,6);
		writebar(image,pos,7);
		break;
	case 9:
		writebar(image,pos,1);
		writebar(image,pos,2);
		writebar(image,pos,3);
		writebar(image,pos,4);
		writebar(image,pos,6);
		break;
	}
}

void writetime(double * const * const image, const int t)
{
        int seqnum = t;
	int a,b,c;
	a = seqnum/100;
	b = (seqnum-a*100)/10;
	c = seqnum-a*100-b*10;
	if (a!=0) writedigit(image,3,a);
	if ((a!=0) || (b!=0))writedigit(image,2,b);
	writedigit(image,1,c);
}
