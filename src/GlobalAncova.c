/*******************************************************************************

  GlobalAncova.c
  Project: GlobalAncova
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/

#include <stdio.h>
#include <R.h>
#include "GlobalAncova.h"
#include "matrix.h"
#include "ludcmp.h"


/*** row.orth2d ***/
/* row.orth2d <- function(xx, D)
        xx %*% (diag(nrow(D)) - D %*% solve(t(D) %*% D) %*% t(D)) */
double * row_orth2d ( double *xx  , int *zxx  , int *sxx ,
                      double *D   , int *zD   , int *sD ,
                      double *res , int *zres , int *sres )
{

  if  ((*sxx!=(*zD)) || ((*zxx)!=(*zres)) || ((*zD)!=(*sres)) ) {
     error("row_orth2d: Clash of Dimension");
  }
   int  i , j,  trans1=1 , trans0=0 ;
   double  *tDD , *tmp , *tmp2,  *tmp3 , det=-1000.0  ;
   tDD = (double *) malloc ( (*sD)*(*sD)*sizeof(double) ) ;
   if ( tDD == NULL )
   { error("Warning in row_orth2d: Memory block of %d bytes unavailable",
      (*sD)*(*sD)*sizeof(double));
   }   
   tmp = (double *) malloc ( (*sD)*(*sD)*sizeof(double) ) ;
   if ( tmp  == NULL )
   { error("Warning in row_orth2d: Memory block of %d bytes unavailable",
      (*sD)*(*sD)*sizeof(double) );
   }
   tmp2 = (double *) malloc ( (*zD)*(*sD)*sizeof(double) ) ;
   if ( tmp2 == NULL )
   { error("Warning in row_orth2d: Memory block of %d bytes unavailable",
      (*zD)*(*sD)*sizeof(double)  );
   }    
   tmp3 = (double *) malloc ( (*zD)*(*zD)*sizeof(double) ) ;
   if ( tmp3 == NULL )
   { error("Warning in row_orth2d: Memory block of %d bytes unavailable",
      (*zD)*(*zD)*sizeof(double) );
   }   
   // R> tmp = solve(t(D) %*% D)
   tDD = dgemm ( D , zD , sD, D, zD , sD ,tDD , sD , sD , &trans1 , &trans0 ) ;
   tmp = matinv ( tDD , sD , tmp , &det ) ;
   if ( fabs(det) < 1e-10 ) {  error ("row_orth2d: System is singular. \n") ; }
   // R> tmp2 = D %*% [solve(t(D) %*% D)]
   tmp2 = dgemm ( D , zD , sD , tmp , sD , sD , tmp2 , zD , sD , &trans0 , &trans0 ) ;
     //  tmp3 = [ D %*% solve(t(D) %*% D) ] %*% t(D))
   tmp3 = dgemm ( tmp2 , zD , sD , D , zD , sD , tmp3 , zD, zD, &trans0 , &trans1 ) ;
  // R> (diag(nrow(D)) - D %*% solve(t(D) %*% D) %*% t(D))
  // Alternative (without unity matrix):
     for ( i=0 ; i<(*zD) ; i++ ) {
       for ( j=0 ; j<(*zD) ; j++ ) {
         if ( j==i) {
           tmp3[i+(j*(*zD))] = 1.0-tmp3[i+(j*(*zD))];
         } else {
           tmp3[i+(j*(*zD))] = -tmp3[i+(j*(*zD))];
         }
       }
     }
  // R> res = xx %*% [(diag(nrow(D)) - D %*% solve(t(D) %*% D) %*% t(D)) ]
   res = dgemm ( xx , zxx , sxx , tmp3 , zD , zD , res , zres , sres , &trans0 , &trans0 ) ;
   free(tDD) ;
   free(tmp) ;
   free(tmp2) ;
   free(tmp3);
   return ( res ) ;
} /** End of row_orth2d **/

void genewiseGA ( double *xx , int *zxx , int *sxx ,
                  double *D_full , int *zD_full , int *sD_full ,
                  double *D_red , int *zD_red , int *sD_red ,
                  double *SS_red_i , double *SS_full_i , double *SS_extra_i ) {

    double s , *R_full ;
    int i,j ;
    R_full = (double *) malloc ( (*zxx)*(*zD_full) * sizeof(double) );
    if ( R_full == NULL )
    { error("Warning in genewiseGA: Memory block of %d bytes unavailable",
      (*zxx)*(*zD_full) * sizeof(double) );
    }  

    // R> R.full <- row.orth2d(xx, D.full)
    row_orth2d ( xx,zxx,sxx, D_full,zD_full,sD_full,  R_full,zxx,zD_full) ;
      
    // R> SS_red_i == NULL (!)  hier: SS_red_i[0]==-1 soll identisch zu NULL sein.
    if ( SS_red_i[0] == -1 ) {
      double *R_red ;
      R_red = (double *) malloc ( (*zxx)*(*zD_red) * sizeof(double) );
      if ( R_red == NULL )
      { error("Warning in genewiseGA: Memory block of %d bytes unavailable",
        (*zxx)*(*zD_red) * sizeof(double) );
      }             
      //R> R.red  <- row.orth2d(xx, D.red)
      row_orth2d ( xx,zxx,sxx, D_red,zD_red,sD_red,  R_red,zxx,zD_red) ;
      //R> SS.red.i <- rowSums(R.red*R.red)
      for ( i=0 ; i<(*zxx) ; i++ ) {
        s = 0.0 ;
        for ( j=0 ; j<(*zD_red) ; j++ ) {
           s += R_red_def(i,j)*R_red_def(i,j) ;
        }
       SS_red_i[i] = s ;
      }
      free(R_red);
    } // if
     // R> # denominator: residual sum of squares in the full model
     // R> SS.full.i <- rowSums(R.full*R.full)
      for ( i=0 ; i<(*zxx) ; i++ ) {
        s = 0.0 ;
        for ( j=0 ; j<(*zD_full) ; j++ ) {
           s += R_full_def(i,j)*R_full_def(i,j) ;
        }
       SS_full_i[i] = s ;
      }
   // R> for ( i=0 ; i<*zxx ; i++ ) { Rprintf("%f  ",SS_full_i[i]);}
   // R>  # nominator: extra residual sum of squares
   // R>  SS.extra.i <- SS.red.i - SS.full.i
   for ( i=0 ; i<*zxx ; i++ ) { SS_extra_i[i] = SS_red_i[i] - SS_full_i[i] ; }
  // R>  return(cbind(nominator=SS.extra.i, denominator=SS.full.i))
  free (R_full);
} /** End of genewiseGA **/

  
void permut_withFperm ( double *D_full , int *zD_full , int *sD_full , double *D_full_perm,
              double *D_red , int *zD_red , int *sD_red ,
              int * N_Subjects , double *rr ,  int *zrr , int *srr ,
              double *SS_red_i , int *perm ,
              int *test_col , int *n_test_col , double *F_value,
              double *DF_full, double *DF_extra , int *ord ,
              int *test_genes , int *num_test_genes , int *n_test_genes ,
              int *count, int *n_singular , int *use_permMat , double *F_perm0 ) {

  int i,ii,j,k,m1,n1 , trans1=1, trans0=0 , np=0 ;
  double *SS_full_i ,  *SS_extra_i , det=0 ;
  double s1 , s2 , *F_perm , *tDD_perm;
  
  SS_full_i  = (double *) malloc ( (*zrr) * sizeof(double)) ;
  if ( SS_full_i == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
     (*zrr) * sizeof(double) );
  }
  SS_extra_i = (double *) malloc ( (*zrr) * sizeof(double)) ;
  if ( SS_extra_i == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*zrr) * sizeof(double));
  }   
  F_perm     = (double *) malloc ( (*n_test_genes) * sizeof(double));
  if ( F_perm == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*n_test_genes) * sizeof(double));
  }   
  tDD_perm   = (double *) malloc ( (*sD_full)*(*sD_full) * sizeof(double)) ;  
  if ( tDD_perm == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*sD_full)*(*sD_full) * sizeof(double));
  }   
  int *ord_perm ;
  if ( (*use_permMat) == 0 ) { 
    ord_perm = (int *) malloc ( (*N_Subjects) * sizeof(int) ) ;
    if ( ord_perm == NULL )
    { error("Warning in permut: Memory block of %d bytes unavailable",
      (*N_Subjects) * sizeof(int));
    }    
    for (i=0;i<(*N_Subjects);i++) ord_perm[i] = i ;
  }
  for ( i=0 ; i<(*perm) ; i++ ) {
    if ( (*use_permMat) == 0 ) {
       sample ( ord_perm , N_Subjects )   ;
    }
    // R> D.full.perm[,test.col] <- D.full[ord[,i]+1, test.col]
    for ( j=0 ; j<(*N_Subjects) ; j++ )  {
       for ( k=0 ; k<(*n_test_col) ; k++ ) {
          if ( (*use_permMat) == 1 ) {
            n1 =  ord_def(j,i) ;
          } else
          {
            n1 =  ord_perm[j] ;
          }
          m1 = test_col[k] ;
          D_full_perm_def(j,m1) = D_full_def(n1,m1) ;
      }
    } // End make D_full_perm
    // R> if ( abs(det(t(D.full.perm)%*%D.full.perm)) > 1e-10 ) {
    tDD_perm = dgemm ( D_full_perm , zD_full , sD_full ,
                       D_full_perm , zD_full , sD_full ,
                       tDD_perm , sD_full , sD_full , &trans1 , &trans0  ) ;
    matdet ( tDD_perm , sD_full , &det ) ;
   
    if ( fabs(det) > 1e-10 ) {
    // R>  genewSS.perm <- genewiseGA(rr, D.full.perm, SS.red.i=SS.red.i)
   genewiseGA ( rr , zrr , srr ,
                D_full_perm , zD_full , sD_full ,
                D_red , zD_red , sD_red ,
                SS_red_i , SS_full_i , SS_extra_i ) ;              
    // R> F.perm <- sapply(test.genes, function(x) sum(genewSS.perm[x,1]) /
    //                             sum(genewSS.perm[x,2])) / (DF.extra / DF.full)
    k = 0 ;
    for ( ii=0 ; ii<(*n_test_genes) ; ii++ ) {
      s1 = 0.0 ; s2 = 0.0 ;
      for ( j=0 ; j<num_test_genes[ii] ; j++ ) {
        s1 += SS_extra_i [ test_genes[k] ] ;
        s2 += SS_full_i  [ test_genes[k] ] ;
        k++ ;
      }
      F_perm[ii] = ( s1 / s2 ) / ( (*DF_extra) / (*DF_full) ) ;
      F_perm0[np] = F_perm[ii] ;
      np++ ;
      // R> count <- count + (F.perm > F.value)
      if ( F_perm[ii] > F_value[ii] )
        { count[ii] ++ ;   }
     }
     } else {
         Rprintf("Warning C: system is singular.\n") ;
         n_singular[0] ++ ;

     }
   } // for *perm
  if ( (*use_permMat) == 0 ) {  free(ord_perm) ; }
  free(tDD_perm);
  free(F_perm) ;
  free(SS_full_i);
  free(SS_extra_i);

} /** End of permut**/


void permut ( double *D_full , int *zD_full , int *sD_full , double *D_full_perm,
              double *D_red , int *zD_red , int *sD_red ,
              int * N_Subjects , double *rr ,  int *zrr , int *srr ,
              double *SS_red_i , int *perm ,
              int *test_col , int *n_test_col , double *F_value,
              double *DF_full, double *DF_extra , int *ord ,
              int *test_genes , int *num_test_genes , int *n_test_genes ,
              int *count, int *n_singular , int *use_permMat ) {

  int i,ii,j,k,m1,n1 , trans1=1 , trans0=0 ;
  double *SS_full_i ,  *SS_extra_i , det=0 ;
  double s1 , s2 , *F_perm , *tDD_perm;
  SS_full_i  = (double *) malloc ( (*zrr) * sizeof(double)) ;
  if ( SS_full_i == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
     (*zrr) * sizeof(double) );
  }
  SS_extra_i = (double *) malloc ( (*zrr) * sizeof(double)) ;
  if ( SS_extra_i == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*zrr) * sizeof(double));
  }   
  F_perm     = (double *) malloc ( (*n_test_genes) * sizeof(double));
  if ( F_perm == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*n_test_genes) * sizeof(double));
  }   
  tDD_perm   = (double *) malloc ( (*sD_full)*(*sD_full) * sizeof(double)) ;  
  if ( tDD_perm == NULL )
  { error("Warning in permut: Memory block of %d bytes unavailable",
    (*sD_full)*(*sD_full) * sizeof(double));
  }   
  int *ord_perm ;
  if ( (*use_permMat) == 0 ) { 
    ord_perm = (int *) malloc ( (*N_Subjects) * sizeof(int) ) ;
    if ( ord_perm == NULL )
    { error("Warning in permut: Memory block of %d bytes unavailable",
      (*N_Subjects) * sizeof(int));
    }    
    for (i=0;i<(*N_Subjects);i++) ord_perm[i] = i ;
  }

  for ( i=0 ; i<(*perm) ; i++ ) {
    if ( (*use_permMat) == 0 ) {
       sample ( ord_perm , N_Subjects )   ;
    }
    // R> D.full.perm[,test.col] <- D.full[ord[,i]+1, test.col]
    for ( j=0 ; j<(*N_Subjects) ; j++ )  {
       for ( k=0 ; k<(*n_test_col) ; k++ ) {
          if ( (*use_permMat) == 1 ) {
            n1 =  ord_def(j,i) ;
          } else
          {
            n1 =  ord_perm[j] ;
          }
          m1 = test_col[k] ;
          D_full_perm_def(j,m1) = D_full_def(n1,m1) ;
      }
    } // End make D_full_perm
    // R> if ( abs(det(t(D.full.perm)%*%D.full.perm)) > 1e-10 ) {
    tDD_perm = dgemm ( D_full_perm , zD_full , sD_full ,
                       D_full_perm , zD_full , sD_full ,
                       tDD_perm , sD_full , sD_full , &trans1 , &trans0  ) ;

    matdet ( tDD_perm , sD_full , &det ) ;
    if ( fabs(det) > 1e-10 ) {
    // R>  genewSS.perm <- genewiseGA(rr, D.full.perm, SS.red.i=SS.red.i)
   genewiseGA ( rr , zrr , srr ,
                D_full_perm , zD_full , sD_full ,
                D_red , zD_red , sD_red ,
                SS_red_i , SS_full_i , SS_extra_i ) ;
    // R> F.perm <- sapply(test.genes, function(x) sum(genewSS.perm[x,1]) /
    // sum(genewSS.perm[x,2])) / (DF.extra / DF.full)
    k = 0 ;
    for ( ii=0 ; ii<(*n_test_genes) ; ii++ ) {
      s1 = 0.0 ; s2 = 0.0 ;
      for ( j=0 ; j<num_test_genes[ii] ; j++ ) {
        s1 += SS_extra_i [ test_genes[k] ] ;
        s2 += SS_full_i  [ test_genes[k] ] ;
        k++ ;
      }
      F_perm[ii] = ( s1 / s2 ) / ( (*DF_extra) / (*DF_full) ) ;
      // R> count <- count + (F.perm > F.value)
      if ( F_perm[ii] > F_value[ii] )
        { count[ii] ++ ;   }
     }
     } else {
         Rprintf("Warning in permut: system is singular.\n") ;
         n_singular[0] ++ ;
     }
   } // for *perm
  if ( (*use_permMat) == 0 ) {  free(ord_perm) ; }
  free(tDD_perm);
  free(F_perm) ;
  free(SS_full_i);
  free(SS_extra_i);

} /** End of permut**/
