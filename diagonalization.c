//#include <AbInitioCode.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#define TRUE 1
#define FALSE 0

/* --------------------------------------------------------------------------------------------------------- */
void diagonalization(
		     double **A,                  /* INPUT n x n matrix to diagonalize */
		     int n,	                         /* INPUT dimension of the matrix A */
		     double *lambda,		 /* OUTPUT eigenvalues vector, dim n */
		     double **y			 /* OUTPUT eigenvectors matrix, dim n x n */
		     ){
  int i,j,sym;sym=TRUE;
  double *data;   data=malloc(n*n*sizeof(double));
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(fabs(A[i][j]-A[j][i])>1.0e-6) {
	//if(sym==TRUE)     printf("%e-%e=%e\n",A[i][j],A[j][i],A[i][j]-A[j][i]);
	sym=FALSE;
      }
      data[i*n+j]=A[i][j];
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if(sym==FALSE){
    //printf("non Symmetric matrix\n");
    gsl_matrix_view m= gsl_matrix_view_array (data, n, n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);
    gsl_eigen_nonsymmv_workspace * w =gsl_eigen_nonsymmv_alloc (n);
    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free (w);
    gsl_eigen_nonsymmv_sort (eval, evec,  GSL_EIGEN_SORT_ABS_ASC);
    {
      int i, j;
      for (i = 0; i < n; i++)      {
	gsl_complex eval_i           = gsl_vector_complex_get (eval, i);
	gsl_vector_complex_view evec_i           = gsl_matrix_complex_column (evec, i);
	//printf ("eigenvalue = %g + %gi\n",	      GSL_REAL(eval_i), GSL_IMAG(eval_i));
	//printf ("eigenvector = \n");
	for (j = 0; j < n; ++j)	{
	  gsl_complex z =              gsl_vector_complex_get(&evec_i.vector, j);
	  //printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
	  //y[i][j]=GSL_REAL(z);
	}
      }
    }
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
  } else {			/* Symmetric matrix */
    //printf("Symmetric matrix\n");
    gsl_matrix_view m     = gsl_matrix_view_array (data, n, n);
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_symmv_workspace * w =     gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec,  GSL_EIGEN_SORT_VAL_ASC);
    {
      int i;
      for (i = 0; i < n; i++)      {
	lambda[i]=gsl_vector_get (eval, i);
	gsl_vector_view evec_i            = gsl_matrix_column (evec, i);
	for (j = 0; j < n; ++j)	y[i][j]=gsl_vector_get(&evec_i.vector, j);
      }
    }
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
  }
  free(data);
}
