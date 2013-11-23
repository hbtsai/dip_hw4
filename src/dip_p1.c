#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#define Size 256

typedef struct complex {
	double r,i;			/* Real and imaginery parts */
} COMPLEX;

/* Some functions for complex arithmetics */

COMPLEX complex_plus(c1, c2)
COMPLEX c1, c2;
{
  COMPLEX c;

  c.r = c1.r + c2.r;
  c.i = c1.i + c2.i;

  return( c );
}

COMPLEX complex_minus(c1, c2)
COMPLEX c1, c2;
{
  COMPLEX c;

  c.r = c1.r - c2.r;
  c.i = c1.i - c2.i;

  return( c );
}

COMPLEX complex_times(c1, c2)
COMPLEX c1, c2;
{
  COMPLEX c;

  c.r = c1.r*c2.r - c1.i*c2.i;
  c.i = c1.r*c2.i + c2.r*c1.i;

  return( c );
}

int complex_print(COMPLEX c) {
   printf("[%f,%fi]\n", c.r, c.i);
   return (1);
}

double complex_mag(COMPLEX c) {
   return (sqrt( pow(c.r,2) + pow(c.i,2)));
}


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
	{
	  return 1;
	}
}


int paint_histogram(int width, int height, unsigned char* image, char* filename)
{

	unsigned char *histoimg = calloc(1, 256*256);
	double histo[256]={};
	int i=0, j=0;
	for(i=0; i<256; i++)
		histo[i]=0;

	for(i=0; i<width*height; i++)
		histo[image[i]]++;

	for(i=0; i<256; i++)
	{
		histo[i]/=65536;
		histo[i]*=255*40;
	}

	for(i=0; i<256; i++)
	{
		for(j=0; j<histo[i]; j++)
		{
			if(j>255)
				break;
			histoimg[i*256+j]=238;
		}

		/*
		if(histo[i]>200)
		{
			fprintf(stderr, " pixel value %d is major\n", i);
		}
		*/

	}

	write_pgm_image(filename, 256, 256, histoimg);
	free(histoimg);
	histoimg=NULL;
	return 0;
}

void DFT( int N, COMPLEX* image_clx,  COMPLEX* dft_clx, int inv)
{
    int u, v, j, k; 

    COMPLEX sum;
    COMPLEX tmp;
    double theta;
    double p;
    theta = 2*M_PI/N;

    if(inv)
        theta=-theta;

    /*
    for(u=0; u<Size; u++)
        for(v=0; v<Size; v++)
            image_clx[u*Size+v].r=image_clx[u*Size+v].r*pow(-1, (u+v));
    */

    for(u=-N/2; u<N/2; u++)
    {
        for(v=-N/2; v<N/2; v++)
        {
            sum.r=0;
            sum.i=0;
            for(j=-N/2; j<N/2; j++) 
            {
                tmp.r = 0;
                tmp.i = 0;
                for(k=-N/2; k<N/2; k++)
                {
                    p = (theta*u*j+theta*v*k);
                    tmp.r = cos(p);
                    tmp.i = -1*sin(p);
                    sum = complex_plus(sum, complex_times( image_clx[(j+N/2)*N+(k+N/2)], tmp ));
                }
            }
            sum.i=sum.i;
            sum.r=sum.r;
            dft_clx[(u+N/2)*N+(v+N/2)] = sum;
        }
    }
}

int main(int argc, char** argv)
{
	FILE *file = NULL;
	// image data array
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
	write_pgm_image("sample1.pgm", Size, Size, Imagedata);

    COMPLEX dft_r[Size*Size]={};
    COMPLEX image_c[Size*Size]={};
    int i=0;
    for(i=0; i<Size*Size; i++)
    {
        image_c[i].r=(double)Imagedata[i];
        image_c[i].i=0; 
    }

    DFT( Size, image_c, dft_r, 0);
    /* return dft raw data*/

    double mag_d[Size*Size]={};
    unsigned char mag[Size*Size]={};

    double mmax=DBL_MIN;
    double mmin=DBL_MAX;
    for(i=0; i<Size*Size; i++)
    {
        mag_d[i]=complex_mag(dft_r[i]);
        if(mag_d[i]>mmax)
            mmax=mag_d[i];
        if(mag_d[i]<mmin)
            mmin=mag_d[i];
    }

    fprintf(stderr, "mmax=%f min=%f\n", mmax, mmin);

    for(i=0; i<Size*Size; i++)
        mag[i]=(unsigned char)(complex_mag(dft_r[i])*254/mmax);

	write_pgm_image("dft_mag.pgm", Size, Size, mag);

    double c = 255/log(1+mmax);
    unsigned char mag_log[Size*Size]={};
    for(i=0; i<Size*Size;i++)
        mag_log[i]=(unsigned char)(c*complex_mag(dft_r[i]));

	write_pgm_image("dft_mag_log.pgm", Size, Size, mag_log);


	exit(0);
	return 0;
}




	
