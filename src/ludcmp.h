/*******************************************************************************

  ludcmp.h
  Project: GlobalAncova 
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/

  int ludcmp(double *a, int n, int *indx, double *d);
  void lubksb(double *a, int n, int *indx, double b[]) ;
  double * matinv ( double *a , int *n , double *y , double * det ) ;
  double * matdet ( double *a , int *n , double * det ) ;
