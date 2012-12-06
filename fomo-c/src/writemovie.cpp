#include "header.h"

#ifdef HAVEMPEG
#include <avcodec.h>
#endif

#include <cmath>
#include <fstream>

char * moviefile="movie.mpeg";

void writemovie(char* imagefile, const double globalmin, const double globalmax)
{
#ifdef HAVEMPEG
	// initialize encoding
	avcodec_init();
	avcodec_register_all();
	AVCodec *codec;
	AVCodecContext *c= NULL;
	int out_size, size, x, y, outbuf_size;
	FILE *f;
	AVFrame *picture;
	uint8_t *outbuf, *picture_buf;
	double t;
// allocate image
	double **image;
	image = (double **)malloc(y_pixel*sizeof(double *));
        for (int row = 0; row < y_pixel; row++)
        {
                image[row] = (double *)malloc(x_pixel*sizeof(double));
// initialize image to black
		for (int j=0; j<x_pixel; j++)
			image[row][j]=0.;
        }
	
	printf("Video encoding\n");
	/* find the mpeg1 video encoder */
	codec = avcodec_find_encoder(CODEC_ID_MPEG2VIDEO);
	if (!codec) {
		fprintf(stderr, "codec not found\n");
		exit(EXIT_FAILURE);
	}
	
	c= avcodec_alloc_context();
	picture= avcodec_alloc_frame();
	
	/* put sample parameters */
	c->bit_rate = 400000;
	/* resolution must be a multiple of two */
	c->width = x_pixel;  
	c->height = y_pixel;
	/* frames per second */
	c->frame_rate = 25;  // this is stupid, why do we need to encode 25 images per second?
	// we only want 1 picture per second, so write 25 times the same
	c->frame_rate_base= 1;
	c->gop_size = 1; /* emit one intra frame every ten frames */
	// increase gop_size to have a smaller movie, but of worse quality
	c->max_b_frames=0;

	/* open it */
	if (avcodec_open(c, codec) < 0) {
		fprintf(stderr, "could not open codec\n");
		exit(EXIT_FAILURE);
	}
	
	/* the codec gives us the frame size, in samples */

	f = fopen(moviefile, "wb");
	if (!f) {
		fprintf(stderr, "could not open %s\n", moviefile);
		exit(EXIT_FAILURE);
	}
	
	/* alloc image and output buffer */
	outbuf_size = 100000;
	outbuf = (uint8_t *)malloc(outbuf_size);
	size = c->width * c->height;
	picture_buf = (uint8_t *)malloc((size * 3) / 2); /* size for YUV 420 */
	
	picture->data[0] = picture_buf;
	picture->data[1] = picture->data[0] + size;
	picture->data[2] = picture->data[1] + size / 4;
	picture->linesize[0] = c->width;
	picture->linesize[1] = c->width / 2;
	picture->linesize[2] = c->width / 2;
	
	// initialize file
	ifstream s(imagefile);
	
	double period = 2*M_PI/frequency.imag();
	for (t=0.; t<nperiods*period; t+=period/nframes)
	{
		// read image from disc
		for (int i=0; i<y_pixel; i++)
			for (int j=0; j<x_pixel; j++)
				s >> image[i][j];
		//writepng(image,t);
		// encode frame
		/* Y */
		for(y=0;y<c->height;y++) 
			for(x=0;x<c->width;x++) 
				picture->data[0][y * picture->linesize[0] + x] = floor((image[y][x]-globalmin)/(globalmax-globalmin)*255+0.5);
		/* Cb and Cr */
		for(y=0;y<c->height/2;y++)
			for(x=0;x<c->width/2;x++) {
				picture->data[1][y * picture->linesize[1] + x] = 255*0.5;
				picture->data[2][y * picture->linesize[2] + x] = 255*0.5;
			}
	
		/* encode the image */
		out_size = avcodec_encode_video(c, outbuf, outbuf_size, picture);
		cout << "encoding frame: t=" << t << " (size=" << out_size << ")"<< endl;
		// write 25 times the same -> artificially 1 frame per second
		for (int i=0; i<25; i++)
			fwrite(outbuf, 1, out_size, f);
	}
	// finish encoding
	/* get the delayed frames */
	for(; out_size; t++) {
		out_size = avcodec_encode_video(c, outbuf, outbuf_size, NULL);
		fwrite(outbuf, 1, out_size, f);
	}

	/* add sequence end code to have a real mpeg file */
	outbuf[0] = 0x00;
	outbuf[1] = 0x00;
	outbuf[2] = 0x01;
	outbuf[3] = 0xb7;
	fwrite(outbuf, 1, 4, f);
	fclose(f);
	free(picture_buf);
	free(outbuf);
	
	avcodec_close(c);
	free(c);
	free(picture);
	printf("\n");
#endif
}
