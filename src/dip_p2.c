#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Size 723

int write_pgm_image(char* filename, int x_dim, int y_dim, unsigned char* image)
{
	unsigned char* y = image;
	FILE* filehandle = NULL;
	filehandle = fopen(filename, "wb");
	if (filehandle) 
	{
		fprintf(filehandle, "P5\n\n%d %d 255\n", x_dim, y_dim);
		fwrite(y, 1, x_dim * y_dim, filehandle);
		fclose(filehandle);
		return 0;
	} 
	else
	  return 1;
}

void build_map(int width, int height, unsigned char *image, unsigned char *image_r, double dir)
{
	double x=0, y=0, z=0;
	double theta=0;
	int radius=width/2;
	int circle = 2*M_PI*radius;
	double mag=0;
	int x_dst=0, y_dst=0;

	for(y=0; y<radius; y++)
	{
		for(x=0; x<circle; x++)
		{

			theta=dir*x/radius;
			mag=radius-y;

			x_dst=lrint(radius+mag*cos(theta));
			y_dst=lrint(radius+mag*sin(theta));
			if(y_dst<0 || y_dst>=height || x_dst<0 || x_dst>=width)
				continue;
			else
				image_r[(int)y*width+(int)x]=image[(int)y_dst*width+(int)x_dst];

		}
	}
}

int main(int argc, char** argv)
{
	FILE *file = NULL;
	unsigned char Imagedata[Size*Size] = {};

	char fname[1024]={};
	if(argv[1] != NULL && strlen(argv[1])>0)
		strcpy(fname, argv[1]);
	else
	{
		fprintf(stderr, "please specify filename of raw input image.\n");
		exit(-1);
	}

	if (!(file=fopen(fname,"rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(Imagedata, sizeof(unsigned char), Size*Size, file);
	fclose(file);

	/* save the original image for comparision */
	write_pgm_image("sample2.pgm", Size, Size, Imagedata);

	unsigned char pan2[Size*Size*2]={};
	build_map(Size, Size, Imagedata, pan2, 3.14);
	write_pgm_image("panorama2.pgm", Size*2, Size/2, pan2);

	exit(0);
	return 0;
}


