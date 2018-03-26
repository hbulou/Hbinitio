#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define FALSE 0
#define TRUE 1
void SteepestDescent(double **A,double *b,int N,double *x);

int main(int nargument,char **argument){
  int i,j;
  FILE *out;
  int N=30;
  double L=10.0,h;
  h=L/N;
  double *R;R=malloc(N*sizeof(double));
  for(i=0;i<N;i++) R[i]=(i+1)*h;
  
  double *b;b=malloc(N*sizeof(double));
  double *x;x=malloc(N*sizeof(double));
  double **A;A=malloc(N*sizeof(double*));

  for(i=0;i<N;i++) A[i]=malloc(N*sizeof(double));

  for(i=0;i<N;i++) {
    A[i][i]=-2.0/h/h;
    if(i>0) A[i][i-1]=1.0/h/h;
    if(i<N-1) A[i][i+1]=1.0/h/h;
  }
  for(i=0;i<N;i++) b[i]=-4.0*R[i]*exp(-2.0*R[i]);
  b[N-1]+=-1.0/h/h;

  
  //Conjugate_Gradient();
  for(i=0;i<N;i++) x[i]=0.0;
  //  for(i=0;i<N;i++) printf("%12.6f ",x[i]);printf("\n");


  SteepestDescent(A,b,N,x);

  out=fopen("U.dat","w+");
  for(i=0;i<N;i++) fprintf(out,"%g %g %g\n",R[i],x[i],1.0-(R[i]+1.0)*exp(-2.0*R[i]));
  fclose(out);

  
  return 0;
}
