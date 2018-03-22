#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0


double dot(int N,double *vec1,double *vec2);


#define GS_ZERO 1.0e-6


double **GramSchmidt(int N,double **set0,int nvec0,double *set1){
  double **out;out=malloc((nvec0+1)*sizeof(double*));
  int i,j,l;
  double alpha;
  for(j=0;j<nvec0+1;j++)  out[j]=malloc(N*sizeof(double));
  for(j=0;j<nvec0;j++){
    alpha=dot(N,set0[j],set1);
    for(l=0;l<N;l++) {
      set1[l]-=alpha*set0[j][l];
      out[j][l]=set0[j][l];
    }
    alpha=dot(N,set1,set1);  alpha=pow(alpha,.5);
    if(fabs(alpha)<GS_ZERO) {
      printf("Error in %s at line %d\n",__FILE__,__LINE__);
      exit(0);
    }
    for(l=0;l<N;l++)     out[nvec0][l]=set1[l]/alpha;
  }
  return out;
}
