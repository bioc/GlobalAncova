/*******************************************************************************

  GlobalAncova.h
  Project: GlobalAncova 
  July 4, 2007
  LMU-IBE, Munich
  Sven Knüppel (MDC Berlin, Department Bioinformatics)

*******************************************************************************/

  #define R_red_def(I,J)        R_red      [(I)+((J)*(*zxx))]
  #define R_full_def(I,J)       R_full     [(I)+((J)*(*zxx))]
  #define D_full_perm_def(I,J)  D_full_perm[(I)+((J)*(*zD_full))]
  #define tD_full_perm_def(I,J) tD_full_perm[(I)+((J)*(*sD_full))]
  #define D_full_def(I,J)       D_full     [(I)+((J)*(*zD_full))]
  #define ord_def(I,J)          ord        [(I)+((J)*(*N_Subjects))]
  
  
  double * row_orth2d ( double *xx  , int *zxx  , int *sxx ,
                        double *D   , int *zD   , int *sD ,
                        double *res , int *zres , int *sres ) ;
                        
  void genewiseGA ( double *xx , int *zxx , int *sxx ,
                    double *D_full , int *zD_full , int *sD_full ,
                    double *D_red , int *zD_red , int *sD_red ,
                    double *SS_red_i , double *SS_full_i , double *SS_extra_i ) ;
  
void permut ( double *D_full , int *zD_full , int *sD_full , double *D_full_perm,
              double *D_red , int *zD_red , int *sD_red ,
              int * N_Subjects , double *rr ,  int *zrr , int *srr ,
              double *SS_red_i , int *perm ,
              int *test_col , int *n_test_col , double *F_value, 
              double *DF_full, double *DF_extra , int *ord , 
              int *test_genes , int *num_test_genes , int *n_test_genes ,
              int *count, int *n_singular , int *use_permMat );
              
void permut_withFperm ( double *D_full , int *zD_full , int *sD_full , double *D_full_perm,
              double *D_red , int *zD_red , int *sD_red ,
              int * N_Subjects , double *rr ,  int *zrr , int *srr ,
              double *SS_red_i , int *perm ,
              int *test_col , int *n_test_col , double *F_value, 
              double *DF_full, double *DF_extra , int *ord , 
              int *test_genes , int *num_test_genes , int *n_test_genes ,
              int *count, int *n_singular , int *use_permMat , double *F_perm0 ) ;

