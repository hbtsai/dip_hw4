/********************************************
** Function prototypes and header infor-   **
** mation for the "fft.c" file.            **
**                                         **
** Author:    Daryle Niedermayer           **
** Date:      Oct 5, 1999                  **
********************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#ifndef _FFT_
#define _FFT_


/***********************************************
** unsigned char: the integer to convert
**      using bitwise operators.
** int: the size of the passed integer in bits
***********************************************/
unsigned char reverse(unsigned char, int);

/***********************************************
** DFT is a one dimensional standard Fourier
** Transform based of the FT definition.
***********************************************/
COMPLEX *DFT(COMPLEX[], int);

/***********************************************
** FFT receives parameters:
** int: the dimension of the given array
**      where array is of size int*int
** COMPLEX**: a pointer to a two dimensional
**      array of complex numbers containing
**      the initial image.
***********************************************/
COMPLEX *FFT(COMPLEX[], int);

#include "fft.c"

/***********************************************
** int: the integer to print out in bits
***********************************************/
void bit_print(int);

#endif

