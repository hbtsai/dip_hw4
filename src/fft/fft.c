/********************************************
** FFT: the subroutine to calculate and    **
** return the 1D-FF transform of a input-  **
** ted image.                              **
**                                         **
** reverse: a subroutine to take an        **
** integer of n bits and reverse it from   **
** high-endian to low-endian order.        **
**                                         **
** Author:    Daryle Niedermayer           **
** Date:      Oct 5, 1999                  **
** Credits:   Thanks to Paul Bourke for    **
**            his contribution and explan- **
**            ations                       **
**            (http://www.swin.edu.au/     **
**             austronomy/pbourke)         **
********************************************/

/********************************************
** REVERSE: Strategy...                    **
** using a single pass, take the most sig- **
** nificant and least significant bits and **
** switch them. For this purpose, leftbit  **
** and rightbit are temporary values of    **
** the bits requiring switching, and left- **
** place and rightplace are the positions  **
** of these bits in the integer.           **
********************************************/
unsigned char reverse(unsigned char input, int n) {
   
   unsigned char leftbit, rightbit, output;
   int i, leftplace, rightplace;

   leftplace = n-1, rightplace = 0;
   output    = 0;

   while (leftplace >= rightplace) {
    
      /* Create bit masks for left and right places */
      leftbit  = 1 << leftplace;
      rightbit = 1 << rightplace;

      /* Use masks to get the value of the bits     */
      leftbit  = input & leftbit;
      rightbit = input & rightbit;

      /* Switch the places for these bits           */
      leftbit  = leftbit >> (leftplace - rightplace);
      rightbit = rightbit << (leftplace - rightplace);

      /* Use masks to insert values into input      */
      output = output | leftbit;
      output = output | rightbit;

      /* Increment/Decrement placement holders      */
      leftplace--;
      rightplace++;

   } /* while (leftplace > rightplace) */
   
   return output;
}

/********************************************
** DFT calculates the 1D-Discrete Fourier  **
** Transform using the definition of a     **
** Fourier Transform:                      **
** F(u)=(1/N)Sum[from x=0 to N-1]f(x)exp^  **
**      (-j2*pi*ux/N)                      **
********************************************/
COMPLEX *DFT(COMPLEX F[], int N) {
   
   int u, x;
   double cosarg = 0;
   double sinarg = 0;
   double W;
   COMPLEX *ft;
   COMPLEX tmp;

   for (u=0; u<N; u++) {
      
      /* Zero out the temporary array  as we go */
      TF[u].r = 0;
      TF[u].i = 0;
      W = 2*PI*u/N;
	  
      for (x=0; x<N; x++) {
         tmp.r = cos(W*x);
         tmp.i = -1*sin(W*x);
   
         /****************************************************
         ** To find the sums in complex form, we must
         ** cross-multiply f(x) with the cos and sin functions
         ** from Euler's Formula:
         ** TF.r*(cosarg - sinarg) + TF.i*(cosarg - sinarg) =
         ** TF.r*cosarg - TF.i*(sinarg) = TF.r
         ** TF.i*cosarg + TF.r*(sinarg) = TF.i 
         ****************************************************/
         TF[u] = complex_plus(TF[u],(complex_times(F[x],tmp))); 
      }
      TF[u].r = TF[u].r/N;
      TF[u].i = TF[u].i/N;
   }
   ft = TF;
   return ft;
}

/********************************************
** FFT calculates the 1D-FFT of a complex  **
** array passed into the function. The ar- **
** ray is expected to be of size N.        **
********************************************/
COMPLEX *FFT(COMPLEX F[], int N) {

   int i, M, j, k, u, i1, i2;
   COMPLEX F2[N];
   COMPLEX Half, W2M;
   int new_index[N];
   int n=(int) log(N) / log(2);

#ifdef _DEBUG_
   printf("The value of n is: %d\n",n);
#endif
 
   Half.r = 0.5;
   Half.i = 0.0;

   /*****************************************
   ** Reverse the array indexes
   *****************************************/
   for (i=0; i<N; i++) {
      new_index[i] = reverse(i,n);
   }

   /*****************************************
   ** Rearrange array elements
   *****************************************/
   for (i=0; i<N; i++) {
      TF[new_index[i]] = F[i];
   }

   /*****************************************
   ** Begin decomposing the array into 
   ** subgroups.
   *****************************************/
   M = 1;			/* The initial length of subgroups  */
   j = N/2;			/* The number of pairs of subgroups */

   /** Begin successive merges for n levels             */
   for (i=0; i<n; i++) {

      /** Merge pairs at the current level              */
      for (k=0; k<j; k++) {

         i1 = k*2*M;        /* Start of first group     */
         i2 = (k*2+1)*M;    /* Start of second group    */

#ifdef _DEBUG_
         printf("Cycling through k=%d... i1 and i2 = %d and %d...\n",k,i1, i2);
#endif
   
         for (u=0; u<M; u++) {
            W2M.r = cos(PI*u/M);
			W2M.i = -sin(PI*u/M);

#ifdef  _DEBUG_
            printf("M is currently: %d\n",M);
	        printf("Cycling through u=%d...\n",u);
            printf("Cycling through i=%d... W2M=%f\n",i,W2M);
#endif
		    
            /* Calculate the values for the first set of numbers */
			F2[u] = complex_times(Half,(complex_plus(TF[i1+u], \
			        complex_times(TF[i2+u],W2M))));


            /*Calculate the values for the second set of numbers */
			F2[u+M] = complex_times(Half,(complex_minus(TF[i1+u], \
                      complex_times(TF[i2+u],W2M))));

	        TF[i1+u] = F2[u];
	        TF[i1+u+M] = F2[u+M];

#ifdef  _DEBUG_
	        printf("F2[%d] is: ",u);
	        complex_print(F2[u]);
	        printf("\n");
	        printf("F2[%d] is: ",u+M);
	        complex_print(F2[u+M]);
	        printf("\n");
	        printf("TF[%d] is: ",i1+u);
	        complex_print(TF[i1+u]);
	        printf("\n");
	        printf("TF[%d] is: ",i1+u+M);
	        complex_print(TF[i1+u+M]);
	        printf("\n");
#endif

	     }
      }
      M = 2*M;
      j = j/2;
   }
   
   return TF;
}


/********************************************
** bit_print is a little utility to print  **
** out an integer bitwise for debugging    **
** purposes.                               **
********************************************/
void bit_print(int a) {
   int i;
   int n = sizeof(int) * CHAR_BIT;
   int mask = 1 << (n - 1);

   for (i=1; i<=n; ++i) {
      putchar(((a & mask) == 0) ? '0' : '1');
      a <<= 1;
      if (i % CHAR_BIT == 0 && i<n)
          putchar(' ');
   }
   printf("\n");
}
