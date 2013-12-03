#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <omp.h>

#define Size 256
#define SKIP_DFT 1
#define SKIP_SQUARE 1

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


int write_raw_image(char* filename, int x_dim, int y_dim, unsigned char* image)
{
	unsigned char* y = image;
	FILE* filehandle = NULL;
	filehandle = fopen(filename, "wb");
	if (filehandle) 
	{
//		fprintf(filehandle, "P5\n\n%d %d 255\n", x_dim, y_dim);
		fwrite(y, 1, x_dim * y_dim, filehandle);
		fclose(filehandle);
		return 0;
	} 
	else
	{
	  return 1;
	}
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


void image_centering(int N, COMPLEX* clx_in, COMPLEX* clx_out)
{
    int u=0, v=0;
    for(u=0; u<N; u++)
        for(v=0; v<N; v++)
            clx_out[u*N+v].r=clx_in[u*N+v].r*pow((-1), (u+v));
    return;
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

    COMPLEX sum[Size*Size]={};
    COMPLEX tmp, tmp2;
    double theta;
    double p;
    theta = 2*M_PI/N;

    if(inv)
        theta=-theta;

    for(u=0; u<N; u++)
    {
        for(v=0; v<N; v++)
        {
            sum[(u)*N+(v)].r=0;
            sum[(u)*N+(v)].i=0;
            for(j=0; j<N; j++) 
            {
                for(k=0; k<N; k++)
                {
                    p = (theta*u*j+theta*v*k);
                    tmp.r = cos(p);
                    if(inv)
                        tmp.i = sin(p);
                    else
                        tmp.i = -1*sin(p);
                    tmp2 = complex_times( image_clx[(j)*N+(k)], tmp );
                    sum[(u)*N+(v)] = complex_plus(sum[(u)*N+(v)], tmp2);
                }
            }
            sum[(u)*N+(v)].r=sum[(u)*N+(v)].r/N;
            sum[(u)*N+(v)].i=sum[(u)*N+(v)].i/N;
            dft_clx[(u)*N+(v)] = sum[(u)*N+(v)];
        }
    }

}



int main(int argc, char** argv)
{
	FILE *file = NULL;
	// image data array
	unsigned char Imagedata[Size*Size] = {};
    int i=0;
    double mmax=DBL_MIN;
    double mmin=DBL_MAX;

    struct timeval start, check;
    gettimeofday(&start, NULL);

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

    //COMPLEX image_c[Size*Size]={};
    COMPLEX* image_c=NULL;
    image_c = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);
    COMPLEX* image_d=NULL;
    image_d = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);
    COMPLEX* dft_r=NULL;
    dft_r = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);
    COMPLEX* dft_amp=NULL;
    dft_amp = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);

#if 1
    for(i=0; i<Size*Size; i++)
    {
        image_c[i].r=(double)Imagedata[i];
        image_c[i].i=0; 
    }

    image_centering(Size, image_c, image_d);

#if !defined(SKIP_DFT)

    DFT( Size, image_d, dft_r, 0);
    /* return dft raw data*/

    /* save raw DFT result to a file*/
    FILE* dft_fp = NULL;
    dft_fp = fopen("dft_r.raw", "wb");
    fwrite(dft_r, sizeof(COMPLEX), Size*Size, dft_fp);
    fclose(dft_fp);
    /*------------------------------*/

#else

    /* read raw DFT from a file */
fprintf(stderr, "%s:%d checkpoint \n", __FILE__, __LINE__);
    FILE* dft_fp = NULL;
    dft_fp = fopen("dft_r.raw", "rb");
    fread(dft_r, sizeof(COMPLEX), Size*Size, dft_fp);
    fclose(dft_fp);

fprintf(stderr, "%s:%d checkpoint \n", __FILE__, __LINE__);

    /* -----------------------------*/

#endif


    double mag_d[Size*Size]={};
    unsigned char mag[Size*Size]={};

    for(i=0; i<Size*Size; i++)
    {
        mag_d[i]=complex_mag(dft_r[i]);
        if(mag_d[i]>mmax)
            mmax=mag_d[i];
        if(mag_d[i]<mmin)
            mmin=mag_d[i];
    }

    fprintf(stderr, "mmax=%f min=%f\n", mmax, mmin);

    unsigned char mag_log[Size*Size]={};
    for(i=0; i<Size*Size; i++)
        mag[i]=(unsigned char)(mag_d[i]*254/mmax);

	write_pgm_image("dft_mag.pgm", Size, Size, mag);
	write_raw_image("dft_sample1_mag.raw", Size, Size, mag_log);

    for(i=0; i<Size*Size;i++)
        mag_log[i]=(unsigned char)((log(mag_d[i]+1)*255/log(1+mmax)));
	write_pgm_image("dft_mag_log.pgm", Size, Size, mag_log);



    mmax=DBL_MIN;
    mmin=DBL_MAX;
    /* create phase here */
    double phase[Size*Size]={};
    for(i=0; i<Size*Size; i++)
    {
        phase[i] = atan2(dft_r[i].i, dft_r[i].r);
        if(mmax<phase[i])
            mmax = phase[i];
        if(mmin>phase[i])
            mmin = phase[i];
    }

    FILE* phase_fp=NULL;
    phase_fp=fopen("phase.raw", "wb");
    fwrite(phase, sizeof(double), Size*Size, phase_fp);
    fclose(phase_fp);

    fprintf(stderr, "out of phase,  mmax =%f mmin=%f \n", mmax, mmin);

    for(i=0; i<Size*Size; i++)
    {
        phase[i]+= (M_PI+1);
        phase[i]*= 180/(2*M_PI);
        if(phase[i]<0)
            fprintf(stderr, "WARNING! phase %f < 0\n", phase[i]);
        if(phase[i]>255)
            fprintf(stderr, "WARNING! phase %f > 255\n", phase[i]);
    }

    unsigned char phase_pgm[Size*Size]={};
    for(i=0; i<Size*Size; i++)
        phase_pgm[i] = (unsigned char)phase[i];

//	write_raw_image("dft_sample1_phase.raw", Size, Size, phase_pgm);
	write_pgm_image("dft_sample1_phase.pgm", Size, Size, phase_pgm);

    unsigned char mag_log2[Size*Size]={};
    for(i=0; i<Size*Size;i++)
        mag_log2[i]=(unsigned char)((log2(mag_d[i]+1)*255/log2(1+mmax)));
	write_pgm_image("dft_mag_log2.pgm", Size, Size, mag_log2);

    unsigned char mag_log10[Size*Size]={};
    for(i=0; i<Size*Size;i++)
        mag_log10[i]=(unsigned char)((log10(mag_d[i]+1)*255/log10(1+mmax)));
	write_pgm_image("dft_mag_log10.pgm", Size, Size, mag_log10);

    gettimeofday(&check, NULL);
      printf(" The DFT operation takes %f seconds\n", (double)((check.tv_sec*1000000+check.tv_usec)- (start.tv_sec * 1000000 + start.tv_usec))/1000000);
    gettimeofday(&start, NULL);



#if !defined(SKIP_SQUARE)

    /* create not centering dft */
    memset(dft_r, 0, sizeof(dft_r));
    DFT( Size, image_c, dft_r, 0);


    memset(mag_d, 0, sizeof(mag_d));
    memset(mag, 0, sizeof(mag));
    mmax=DBL_MIN;
    mmin=DBL_MAX;
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
        mag[i]=(unsigned char)(mag_d[i]*254/mmax);

	write_pgm_image("dft_rec_mag.pgm", Size, Size, mag);

    memset(mag_log, 0, sizeof(mag_log));
    for(i=0; i<Size*Size;i++)
        mag_log[i]=(unsigned char)((255/log(1+mmax))*log(mag_d[i]+1));
	write_pgm_image("dft_rec_mag_log.pgm", Size, Size, mag_log);

    memset(mag_log2, 0, sizeof(mag_log2));
    for(i=0; i<Size*Size;i++)
        mag_log2[i]=(unsigned char)((255/log2(1+mmax))*log2(mag_d[i]+1));
	write_pgm_image("dft_rec_mag_log2.pgm", Size, Size, mag_log2);

    memset(mag_log10, 0, sizeof(mag_log10));
    for(i=0; i<Size*Size;i++)
        mag_log10[i]=(unsigned char)((255/log10(1+mmax))*log10(mag_d[i]+1));
	write_pgm_image("dft_rec_mag_log10.pgm", Size, Size, mag_log10);

    gettimeofday(&check, NULL);
      printf(" create non-centering DFT takes %f seconds\n",(double) ((check.tv_sec * 1000000 + check.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec))/1000000);
    gettimeofday(&start, NULL);
#endif

#endif

    unsigned char dft_mag_data[Size*Size]={};
    //unsigned char dft_phase_data[Size*Size]={};
    double dft_phase_data[Size*Size]={};
	if (!(file=fopen("dft_sample1.raw","rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(dft_mag_data, sizeof(unsigned char), Size*Size, file);
	fclose(file);
	if (!(file=fopen("phase.raw","rb")))
	{
		fprintf(stderr, "Cannot open file!\n");
		exit(1);
	}
	fread(dft_phase_data, sizeof(double), Size*Size, file);
	fclose(file);

    mmax=DBL_MIN;
    mmin=DBL_MAX;
    //double phase_dft[Size*Size]={};
    for(i=0; i<Size*Size; i++)
    {
     //   phase_dft[i] = ((double)dft_phase_data[i]*2*M_PI/180)-1-M_PI;
        if(mmax<dft_phase_data[i])
            mmax = dft_phase_data[i];
        if(mmin>dft_phase_data[i])
            mmin = dft_phase_data[i];
    }

    fprintf(stderr, "out of recovered phase,  mmax =%f mmin=%f \n", mmax, mmin);

    COMPLEX dft_recover[Size*Size]={};
    for(i=0; i<Size*Size; i++)
    {
        dft_recover[i].r =  dft_mag_data[i]*cos(dft_phase_data[i]);
        dft_recover[i].i =  dft_mag_data[i]*sin(dft_phase_data[i]);
    }

    COMPLEX* f2_clx=NULL;
    f2_clx = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);

    //COMPLEX* f2_clx_r=NULL;
    //f2_clx_r = (COMPLEX*) calloc(sizeof(COMPLEX), Size*Size);

    //image_centering(Size, dft_amp, f2_clx_r);
    DFT(Size, dft_recover, f2_clx, 1);

    mmax=DBL_MIN;
    mmin=DBL_MAX;
    for(i=0; i<Size*Size; i++)
    {
        if(mmax<f2_clx[i].r)
            mmax = f2_clx[i].r;
        if(mmin>f2_clx[i].r)
            mmin = f2_clx[i].r;
    }

    fprintf(stderr, "out of recovered data,  mmax =%f mmin=%f \n", mmax, mmin);

    unsigned char f2[Size*Size]={};
    for(i=0; i<Size*Size; i++)
        f2[i]=(unsigned char)f2_clx[i].r;

	write_pgm_image("sample_from_dft.pgm", Size, Size, f2);

    gettimeofday(&check, NULL);
      printf(" convert DFT back takes %f seconds\n",(double) ((check.tv_sec * 1000000 + check.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec))/1000000);

   // gettimeofday(&start, NULL);
   
   free(image_c);
   free(image_d);
   free(dft_r);
   free(dft_amp);
   free(f2_clx);
   //free(f2_clx_r);

	exit(0);
	return 0;
}




	
