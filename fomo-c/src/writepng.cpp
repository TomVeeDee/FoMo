#include "header.h"
#include <cstdlib>
#include <cmath>

#ifdef HAVEPNG
#include <png.h>
#endif

#include <sstream>
#include <string>

#ifdef HAVEPNG
void user_error_fn(png_structp, png_const_charp error_msg)
{
	cerr << error_msg <<endl;
}

void user_warning_fn(png_structp, png_const_charp warning_msg)
{
	cerr << warning_msg << endl;
}
#endif

int writepng(double * const * const ccd, const double t)
{
#ifdef HAVEPNG
	double period = 2*M_PI/frequency.imag();
	int seqnum = (int)floor(t*nframes/period+0.5);
	stringstream ss;
	ss << "outputfile.t" << seqnum << ".png";
	string graphic_out=ss.str();
	FILE *fp = fopen((char *)graphic_out.c_str(), "wb");
	char user_error_ptr[] = "PngEncoder";
	if (!fp)
	{
		return (EXIT_FAILURE);
	}
	png_structp png_ptr = png_create_write_struct
	(PNG_LIBPNG_VER_STRING, (png_voidp)user_error_ptr,
	user_error_fn, user_warning_fn);
	if (!png_ptr)
	{
		return (EXIT_FAILURE);
	}
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_write_struct(&png_ptr,
		(png_infopp)NULL);
		return (EXIT_FAILURE);
	}
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fp);
		return (EXIT_FAILURE);
	}
	png_init_io(png_ptr, fp);
//
//	Setting the png_info
//
// for trace ccdwidth = ccdheight = 1024
	const png_uint_32 ccdwidth = x_pixel;
	const png_uint_32 ccdheight = y_pixel;
	const png_byte bit_depth = 8;
	const double maxval = 255;
	png_set_IHDR(png_ptr, info_ptr, ccdwidth, ccdheight,
	bit_depth, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
	PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
// convert func to image
	png_byte **image;
// allocate image
	image = (png_byte **)malloc(ccdheight*sizeof(png_byte *));
	for (int row = 0; row < ccdheight; row++)
	{
		image[row] = (png_byte *)malloc(ccdwidth*sizeof(png_byte));
	}
// fill image with relevant data
	int i_max, i_min, j_max, j_min;
	double max = findmax(ccd, &i_max, &j_max);
	double min = findmin(ccd, &i_min, &j_min);
	png_byte temp;
	for (int y=0; y<y_pixel; y++)
		for (int x=0; x<x_pixel; x++)
	{
// x is r coordinate
// y is z coordinate
		temp = (png_byte)floor((ccd[y][x]-min)/(max-min)*maxval+0.5);
		image[y][x] = temp;
	}
// write out image
	png_write_image(png_ptr, image);
// finish image
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fp);
#endif
	return EXIT_SUCCESS;
}
