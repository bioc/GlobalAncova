/*******************************************************************************

  ludcmp.c
  Project: GlobalAncova 
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/
/*
  Quelle: http://jrfonseca.planetaclix.pt/projects/old/xlmat.html
  download: 12.03.2007
  Book: "Numerical Recipes in C"
*/

#include <stdlib.h>
#include "R.h"
#include "ludcmp.h"

#define TINY 1.0e-20	// A small number.

/* ludcmp
 *
 * Given a matrix a[n * n], this routine replaces it by the LU decomposition of a rowwise
 * permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
 * indx[n] is an output vector that records the row permutation elected by the partial
 * pivoting; d is output as  1 depending on whether the number of row interchanges was even
 * or odd, respectively. This routine is used in combination with lubksb to solve linear equations
 * or invert a matrix.
 */
int ludcmp(double *a, int n, int *indx, double *d)
{
	int i, imax=0, j, k;
	double big, dum, sum, temp;
	double *vv;	// vv stores the implicit scaling of each row.

	vv = (double *) malloc(n * sizeof(double));

	*d = 1.0;	// No row interchanges yet.

	// Loop over rows to get the implicit scaling information.
	for (i = 0; i < n; i++)
	{
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i * n + j])) > big)
				big = temp;

		if (big == 0.0)
		{
			// Singular matrix in routine ludcmp
			free(vv);
			return 0;
		}

		// No nonzero largest element.
		vv[i] = 1.0 / big;       // Save the scaling.
	}

	// This is the loop over columns of Crout's method.
	for (j = 0; j < n; j++)
	{
		// This is equation (2.3.12) except for i = j.
		for (i = 0; i < j; i++)
		{
			sum = a[i * n + j];
			for (k = 0; k < i; k++)
				sum -= a[i * n + k] * a[k * n + j];
			a[i * n + j] = sum;
		}

		// Initialize for the search for largest pivot element.
		big = 0.0;

		for (i = j; i < n; i++)
		{
			// This is i = j of equation (2.3.12) and i = j + 1 : : : N of equation (2.3.13).
			sum = a[i * n + j];
			for (k = 0; k < j; k++)
				sum -= a[i * n + k] * a[k * n + j];
			a[i * n + j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big)
			{
				// Is the figure of merit for the pivot better than the best so far?
				big = dum;
				imax = i;
			}
		}
		// Do we need to interchange rows?
		if (j != imax)
		{
			// Yes, do so...
			for (k = 0; k < n; k++)
			{
				dum = a[imax * n + k];
				a[imax * n + k] = a[j * n + k];
				a[j * n + k] = dum;
			}

			// ...and change the parity of d.
			*d = -(*d);

			// Also interchange the scale factor.
			vv[imax] = vv[j];

		}

		// If the pivot element is zero the matrix is singular (at least to the precision of the
		// algorithm). For some applications on singular matrices, it is desirable to substitute
		// TINY for zero.
		indx[j] = imax;  
		if (a[j * n + j] == 0.0) 
			a[j * n + j] = TINY;

		if (j != n - 1)
		{
			// Now, finally, divide by the pivot element.
			dum = 1.0 / (a[j * n + j]);
			for (i = j + 1; i < n; i++)
				a[i * n + j] *= dum;
		}
		
	} // Go back for the next column in the reduction.

	free(vv);

	return 1;
}

/* lubksb
 *
 * Solves the set of n linear equations A . X = B. Here a[n * n] is input, not as the matrix
 * A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
 * as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
 * B, and returns with the solution vector X. a, n, and indx are not modified by this routine
 * and can be left in place for successive calls with diferent right-hand sides b. This routine takes
 * into account the possibility that b will begin with many zero elements, so it is elcient for use
 * in matrix inversion.
 */
void lubksb(double *a, int n, int *indx, double b[])
{
	int i, ii = -1, ip, j;
	double sum;

	for (i = 0; i < n; i++)
	{
		// When ii is set to a positive value, it will become the
		// index of the first nonvanishing element of b. We now
		// do the forward substitution, equation (2.3.6). The
		// only new wrinkle is to unscramble the permutation
		// as we go.

		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];

		if (ii >= 0)
			for (j = ii; j < i; j++)
				sum -= a[i * n + j] * b[j];
		else
			if (sum)
				// A nonzero element was encountered, so from now on we
				// will have to do the sums in the loop above.
				ii = i;

		b[i] = sum;
	}

	for (i = n - 1; i >= 0; i--)
	{
		// Now we do the backsubstitution, equation (2.3.7).
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum -= a[i * n + j] * b[j];

		// Store a component of the solution vector X.
		b[i] = sum / a[i * n + i];
	}
}  /****END****/




double * matinv ( double *a , int *n , double *y , double * det ) {
//#define N
  double  d, *col ;
  int i,j,*indx ;

  if ( ( col = (double *) malloc((*n) * sizeof(double)) ) == NULL )
  { error("No memory allocation."); }
  
  if ( ( indx = (int *) malloc((*n) * sizeof(int)) ) == NULL )
  { error("No memory allocation."); }  
  
  /* compute LU decomposition */
  ludcmp(a,(*n),indx,&d);
  
  /*compute determinante*/
  (*det) = 1.0 ;
  for (  i=0 ; i<(*n) ; i++ ) { (*det) *= a[i*(*n)+i ] ; }


  /* compute inverse */   
  for ( j=0;j<(*n);j++) {     
    for ( i=0 ; i<(*n) ; i++ ) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,(*n),indx,col);
    for (i=0;i<(*n);i++) y[i*(*n)+j]=col[i] ;  
  }

  free (col) ;
  free (indx) ;
  
  return ( y ) ;
} /****End matinv ****/


double * matdet ( double *a , int *n , double * det ) {
//#define N
  double  d ;
  int i,*indx ;
 
  if ( ( indx = (int *) malloc((*n) * sizeof(int)) ) == NULL )
  { error("No memory allocation."); }  
  
  /* compute LU decomposition */
  ludcmp(a,(*n),indx,&d);
  /*compute determinante*/
  (*det) = 1.0 ;
  for (  i=0 ; i<(*n) ; i++ ) { (*det) *= a[i*(*n)+i ] ; }
  return ( det ) ;
} /** End of matdet **/