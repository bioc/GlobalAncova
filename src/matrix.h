/*******************************************************************************

  matrix.h
  Project: GlobalAncova 
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/
 
  #define A_def(I,J)         a [(I)+((J)*(*za))]
  #define B_def(I,J)         b [(I)+((J)*(*zb))]
  #define C_def(I,J)         c [(I)+((J)*(*zc))]
  
  
  /**********************************************************
  *                                                        *
  * Functions                                              *
  *                                                        *
  **********************************************************/
  
  void sample ( int *perm , int *N ) ;
  
  /**********************************************************
  *                                                        *
  * matrix operations                                      *
  *                                                        *
  **********************************************************/
  /* Display Matrix A */
  void matprint ( double *a , int *za , int *sa ) ;
  void matprint_integer ( int *a , int *za , int *sa ) ;
  /* C = A + B */
  double * matadd ( double *a, double *b, int *x, int *y, double *c ) ;
  /* C = A - B */
  double * matsub ( double *a, double *b, int *x, int *y, double *c ) ;
  /* matmult is simple multiplication : C = A * B */
  double * matmult ( double * a , int *za , int *sa ,
                    double * b , int *zb , int *sb ,
                    double * c , int *zc , int *sc ) ;
  /* C = transpose(A)*/
  double * mattrans ( double *a , int *za , int *sa , 
                     double *c , int *zc , int *sc ) ;
  /* dgemm: better algorithm for matrix multiplication (CLAPACK) 
   C = A * B */
  double * dgemm  ( double * a , int *za , int *sa ,
                   double * b , int *zb , int *sb ,
                   double * c , int *zc , int *sc ,
                   int *transa , int *transb ) ;
  /* C = skalar * A */
  double * matskalar ( double * a , int *za , int *sa ,
                      double *skalar ,
                      double * c , int *zc , int *sc ) ;
                      
