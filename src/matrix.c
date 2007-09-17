/*******************************************************************************

  matrix.c
  Project: GlobalAncova 
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/

#include <stdio.h>
#include <time.h> /* to need for srand() */
#include <R.h>
#include "matrix.h"
#include "ludcmp.h"

/***** Globale Variable *****/
int seed = 123456  ;
                                            
/***** Functions  *****/

/* sample: random permutation from [0,(*N-1)] */
void sample ( int *perm , int *N )
{
   int i , z1 , z2 , h  ; 
   time_t tim ;
   /* (*N) times:  swap of two elements  */ 
   for ( i=0;i<(*N);i++ ) perm[i]=i;
     /* change seed */  
  time(&tim)  ;
  seed += (int) tim ;
	srand ( (int) (seed) ) ;
                        
   for ( i=0 ; i<(*N) ; i++ )
   {
    z1 = ( (int) (rand() % (*N)) ) ;
    z2 = ( (int) (rand() % (*N)) ) ;
    h = perm[z1] ; 
    perm[z1] = perm[z2] ;
    perm[z2] = h ;
   }  
   /* Print */
   /* for ( i=0;i<*N;i++ ) Rprintf ( " %i " , perm[i]) ; 
      Rprintf ( "\n") ;
   */ 
} /***** End of sample *****/

void matprint ( double *a , int *za , int *sa ){

  int i,j ;
  Rprintf ( "\n" ) ;
  for ( i=0 ; i<*za ; i++ ) {
     for ( j=0 ; j<*sa ; j++ ) {
         Rprintf ("%f " , A_def(i,j) ) ;
      } 
     Rprintf ( "\n" ) ;
  }

} /* End matprint*/
void matprint_integer ( int *a , int *za , int *sa ){

  int i,j ;
  printf ( "\n" ) ;
  for ( i=0 ; i<*za ; i++ ) {
     for ( j=0 ; j<*sa ; j++ ) {
         Rprintf ("%i " , A_def(i,j) ) ;
      } 
     Rprintf ( "\n" ) ;
  }

} /* End matprint_integer*/
 
double * matadd ( double *a, double *b, int *x, int *y, double *c ) { 
  int i  ;  
  for ( i=0;i<(*x)*(*y);i++ ) { c[i] = a[i] + b[i] ; }
  return ( c ) ;
} /* End mat_add */

double * matsub ( double *a, double *b, int *x, int *y, double *c ) { 

  int i  ;  
  for ( i=0;i<(*x)*(*y);i++ ) { c[i] = a[i] - b[i] ; }
  return(c) ;
} /* End mat_add */

/* Multiplication of two matrizes
  (za,sa)-matrix A
  (zb,sb)-matrix B */
double * matmult ( double * a , int *za , int *sa ,
                   double * b , int *zb , int *sb ,
                   double * c , int *zc , int *sc )
{
  int  i , j , k  ;
  double  s=0.0 ;    
  if ( ((*sa)!=(*zb)) || ((*za)!=(*zc))|| ((*sa)!=(*sc))  ) { 
    error ( "matmult: Clash of dimension." ) ;
  }
  for (i=0; i<(*za) ; i++)
  {
    for (j=0; j<(*sb) ; j++)
    {
      s = 0.0 ;
      for (k=0; k<(*sa); k++)
      {
		      s += ( A_def(i,k) ) * ( B_def(k,j) ) ;
	    }
	    C_def(i,j) = s ;
	  }
  }
  return ( c ) ;
} /* End matmult */


double * mattrans ( double *a , int *za , int *sa , 
                    double *c , int *zc , int *sc ) {
   if ( ((*za)!=(*sc)) || ((*sa)!=(*zc)) ) {
     error ("Error in matskalar: Clash of Dimension") ;
   }
/* !!! */
  int i,j ;   
/*  for ( int i=0 ; i<(*za) ; i++ ) { 
    for ( int j=0 ; j<(*sa) ; j++) { C_def(j,i) = A_def(i,j) ;  }   */
  for ( i=0 ; i<(*za) ; i++ ) {
    for ( j=0 ; j<(*sa) ; j++) { C_def(j,i) = A_def(i,j) ;  }
/* !!! */
  }
  return ( c ) ;
} /* End mattrans*/


/* Multiplication of two matrizes
  (za,sa)-matrix A
  (zb,sb)-matrix B 

  Modifizierte Version von BLAS (June 20, 2007):
  Quelle: http://www.netlib.org/clapack/
   
  input :
  
  transa    integer  1=transpose matrix a ; 0=not transpose
  transb    integer  1=transpose matrix b ; 0=not transpose
  
*/
double * dgemm  ( double * a , int *za , int *sa ,
                  double * b , int *zb , int *sb ,
                  double * c , int *zc , int *sc ,
                  int *transa , int *transb  )
{
  int j , i__1 , i__2 , i__3 , i__ , l, n,m,k ;
  double temp;
  m = *zc ;
  n = *sc ;
  /*** C = A * B ***/
  if ( (*transa==0) && (*transb==0) ) {
    k = *sa ; 
    i__1 = n ;     
  	for (j =0; j < i__1; ++j) {
      i__2 = m;
      for (i__ = 0; i__ < i__2; ++i__) {
        C_def(i__,j) = 0.;
      }
      i__2 = k;
      for (l = 0; l < i__2; ++l) {
        if (B_def(l,j) != 0.) {
          temp = B_def(l,j);
          i__3 = m;
          for (i__ = 0; i__ < i__3; ++i__) {
            C_def(i__,j) = C_def(i__,j) + temp * A_def(i__,l);
          }
        }
      } 
    }
  } /*---------------------------------------------------*/
  /*** C = t(A) * B ***/  
  if ( (*transa==1) && (*transb==0) ) {
    k = *za ;   
    i__1 = n;
	  for (j = 0; j < i__1; ++j) {
  		i__2 = m;
  		for (i__ = 0; i__ < i__2; ++i__) {
    		temp = 0.;
    		i__3 = k;
    		for (l = 0; l < i__3; ++l) {
      		temp += A_def(l,i__) * B_def(l,j);
    		}
    		C_def(i__,j) = temp;
  		}                            
	  }   
	}  /*---------------------------------------------------*/
	
  /*** C = A * t(B) ***/  
  if ( (*transa==0) && (*transb==1) ) {
    k = *sb ; 
    i__1 = n;
    for (j = 0; j < i__1; ++j) {
      	i__2 = m;
      	for (i__ = 0; i__ < i__2; ++i__) {
        	C_def(i__,j) = 0.;
      	}
    	i__2 = k;
    	for (l = 0; l < i__2; ++l) {
      	if (B_def(j,l) != 0.) {
        	temp = B_def(j,l);
        	i__3 = m;
        	for (i__ = 0; i__ < i__3; ++i__) {
          	C_def(i__,j) = C_def(i__,j) + temp * A_def(i__,l);
        	}
        }
      }
    }
  } /*---------------------------------------------------*/
  /*** C = t(A) * t(B) ***/  
  if ( (*transa==1) && (*transb==1) ) {
    k = *sb ; 
    i__1 = n ;
    for (j = 0; j < i__1; ++j) {
      i__2 = m;
      for (i__ = 0; i__ < i__2; ++i__) {
        temp = 0.;
        i__3 = k;
        for (l = 0; l < i__3; ++l) {
          temp += A_def(l,i__) * B_def(j,l);
        }
          C_def(i__,j) = temp;
      }
    }
  }  /*---------------------------------------------------*/
  return ( c )  ;
} /** End of dgemm **/

 double * matskalar ( double * a , int *za , int *sa ,
                      double *skalar ,
                      double * c , int *zc , int *sc ) {
                      
   if ( ((*za)!=(*zc)) || ((*sa)!=(*sc)) ) {
     error ("Error in matskalar: Clash of Dimension") ;
   }
   int i ;
   for ( i=0 ; i<(*za)*(*sa) ; i++ ) { c[i] = (*skalar) * a[i] ; }
   return ( c ) ;                      
                      
 } /** End of matskalar **/

