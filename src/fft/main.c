/*
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/********************************************
** This program determines the Fast Four-  **
** ier Transform (FFT) of a supplied image.**
** The program does this by calling the    **
** routines referenced in the fft.h file.  **
**                                         **
** Author:    Daryle Niedermayer           **
** Date:      Oct 5, 1999                  **
**                                         **
** Credits:   Derived from the template    **
**            supplied by Xue Dong Yang    **
**            Oct 2, 1999                  **
**                                         **
** To see the FFT run though its procedure **
** define _DEBUG_                          **
********************************************/
#undef _DEBUG_

/********************************************
** ROWS and COLS must be defined for the   **
** correct image size.                     **
********************************************/
#define ROWS 64
#define COLS 64

/********************************************
** This same procedure can be used to gen  **
** erate either the FFT or DFT algorithms. **
** Define either _FFT_OPT_ or _DFT_OPT_    **
** here.                                   **
** _FFT_OPT_ is the default.               **
********************************************/
#define _FFT_OPT_

#ifdef _DFT_OPT_
#undef _FFT_OPT_
#else
#define _FFT_OPT_
#undef _DFT_OPT_
#endif

/********************************************
** Other files and definitions to include  **
********************************************/
#include "img.h"		/* Assume it is at the current directory */
#include "fft.h"		/* Fast Fourier Transform routines       */


unsigned char in_buf[ROWS][COLS];	/* Buffer for the input image */
unsigned char out_buf[ROWS][COLS];	/* Buffer for the output Fourier
					   spectrum image */

main(argc, argv)
int argc;
char **argv;
{
  FILE *fin, *fout;
  int  i, j, n;                     /* i and j are loop counters */
  									/* n is the log(N) of the    */
									/*   sample size.            */
  float max, scale_factor;			/* used for scaling the im-  */
  									/* ages brightness.          */
  COMPLEX *ft;						/* pointer to the array con- */
  									/* taining the results of    */
									/* the ft.                   */
  COMPLEX F_buffer[ROWS];			/* Buffer of complex numbers */
  									/* to hand off to FT function*/
  int half_width, half_height;
  unsigned char swap_buf[ROWS][COLS]; 
  									/* Shifted image for output  */


  /* Check the number of arguments in the command line. */
  if (argc != 3) {
    fprintf(stderr, "Usage: %s in.img out.img\n", argv[0]);
    exit(1);
  }

  /* Open the input image file */
  if ((fin = fopen(argv[1], "rb")) == NULL) {
    fprintf(stderr, "ERROR: Cann't open input image file %s\n", argv[1]);
    exit(1);
  }

  /* Open the output image file */
  if ((fout = fopen(argv[2], "wb")) == NULL) {
    fprintf(stderr, "ERROR: Cann't open output image file %s\n", argv[2]);
    exit(1);
  }

  /* Load the input image */
  if ((n=fread(in_buf, sizeof(char), ROWS*COLS, fin)) < ROWS*COLS*sizeof(char)){
    fprintf(stderr, "ERROR: Read input image file %s error)\n", argv[1]);
    exit(1);
  }

  /* Convert the real image into complex image first. */
  for (i=0; i<ROWS; i++)
    for (j=0; j<COLS; j++) {
      F[i][j].r = (float)in_buf[i][j]; /* The real part = input image */
      F[i][j].i = 0.0;			/* The imaginery part = 0 */
    }

  /* ========================= PASS 1 ==============================
	Applying the 1D FFT function to each columun of the input
	image, and save the intermediate results into a temparory
	array F1.
     =============================================================== */

  for (j=0; j<COLS; j++) {

    /* Copy a column from array F into the temporary vector TF*/
	for (i=0; i<ROWS; i++) {
	   F_buffer[i].r = F[i][j].r;
	   F_buffer[i].i = F[i][j].r;
    }

    /* Call the corresponding function to compute the Fourier transform 
	   depending on whether the _DFT_OPT_ or _FFT_OPT_ flag is set   */
#ifdef _DFT_OPT_
	ft = DFT(F_buffer,ROWS);
#else
	ft = FFT(F_buffer,ROWS);
#endif

    /* Copy the returned result into a temporary array F1 */
	for (i=0; i<ROWS; i++) {
       F1[i][j].r = ft[i].r;
	   F1[i][j].i = ft[i].i;
	}
  }

  /* ========================= PASS 2 ==============================
	Applying the 1D FFT function to each row of the intermediate
	results, and save the result back into array F.
     =============================================================== */

  for (i=0; i<ROWS; i++) {
    /* Copy a row from array F1 into the temporary vector TF*/

	for (j=0; j<COLS; j++) {
	   F_buffer[j].r = F1[i][j].r;
	   F_buffer[j].i = F1[i][j].r;
    }

    /* Call the corresponding function to compute the Fourier transform 
	   depending on whether the _DFT_OPT_ or _FFT_OPT_ flag is set   */
#ifdef _DFT_OPT_
	ft = DFT(F_buffer,ROWS);
#else
	ft = FFT(F_buffer,ROWS);
#endif

    /* Copy the returned result back into array F */
	for (j=0; j<COLS; j++) {
       F[i][j].r = ft[j].r;
	   F[i][j].i = ft[j].i;
	}
  }

  /* Compute the Fourier spectrum |F(u,v)| and save it into the
     output image buffer out_buf.
     Note that proper scaling for the values is needed in order
     to obtain a reasonably bright image. Refer to the discussion
     on page 92, the last paragraph.
  */

  max = 0;

  for (i=0; i<ROWS; i++) {
    for (j=0; j<COLS; j++) {
      if (complex_mag(F[i][j]) > max) max = complex_mag(F[i][j]);
    }
  }

#ifdef _DEBUG_
  printf ("The min and maximum values are: %f and %f \n",min,max);
#endif

  scale_factor = 255/(log(1 + max));
  for (i=0; i < ROWS; i++) {
    for (j=0; j < COLS; j++) {
      out_buf[i][j] = (unsigned char) scale_factor*(log(1 + complex_mag(F[i][j])));
    }
  }

  /* Swap sectors of the image to reposition the image to the center of the
     screen.
  */

  half_height = (int)ROWS/2;
  half_width  = (int)COLS/2;

  for (i=0; i<half_height; i++) {
    for (j=0; j<half_width; j++) {
	  swap_buf[half_height+i][half_width+j] = out_buf[i][j];
	  swap_buf[i][j]             = out_buf[(half_height-1)-i][(half_height-1)-j];
	  swap_buf[i][j+half_width]  = out_buf[(half_height-1)-i][j];
	  swap_buf[half_height+i][j] = out_buf[i][(half_width-1)-j];
    }
  }

  /* Save the output Fourier spectrum image */
  if ((n=fwrite(swap_buf, sizeof(char), ROWS*COLS, fout)) < ROWS*COLS*sizeof(char)) {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[1]);
    exit(1);
  }

  fclose(fin);
  fclose(fout);
}

