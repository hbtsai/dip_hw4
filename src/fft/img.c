/* =======================================
 *
 *	IMAGE FORMAT DEFINITIONS
 *
 * Dr. Xue Dong Yang
 * =======================================
 */

#include <stdio.h>
#include <math.h>
#ifndef ROWS
#define ROWS 256
#endif

#ifndef COLS
#define COLS 256
#endif

unsigned char in_img[ROWS][COLS];	/* Input raw raster image */
unsigned char out_img[ROWS][COLS];	/* Output raster image */

typedef struct complex {
	float r,i;			/* Real and imaginery parts */
	} COMPLEX;

COMPLEX F[ROWS][COLS];			/* Arrays used FFT computation */
COMPLEX F1[ROWS][COLS];
COMPLEX TF[COLS];

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
   printf("[%f,%f]", c.r, c.i);
   return (1);
}

double complex_mag(COMPLEX c) {
   return (sqrt(pow(c.r,2) + pow(c.i,2)));
}

