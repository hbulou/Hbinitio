#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define FALSE 0
#define TRUE 1
void ConjugateGradient(double **A,double *b,int N,double *x){

  double *Ax;Ax=malloc(N*sizeof(double));
  double nrj,dnrj,nrj_old;
  double *r;r=malloc(N*sizeof(double));
  double *p0;p0=malloc(N*sizeof(double));
  double *p1;p1=malloc(N*sizeof(double));
  double *Ar;Ar=malloc(N*sizeof(double));
  double *Ap0;Ap0=malloc(N*sizeof(double));
  double *Ap1;Ap1=malloc(N*sizeof(double));
  int loop,i,j;
  int loopmax=1000;
  int cvg=FALSE;
  double num_beta,den_beta,beta;
  double num_alpha,den_alpha,alpha;
  FILE *out;
  double tol=1.0e-6;

  
  while(loop<loopmax && cvg==FALSE){
    for(i=0;i<N;i++){
      Ax[i]=0.0;
      for(j=0;j<N;j++) Ax[i]+=A[i][j]*x[j];
    }

    if(loop>0) nrj_old=nrj;
    nrj=0.0;
    for(i=0;i<N;i++) nrj+=x[i]*Ax[i];
    nrj*=.5;
    for(i=0;i<N;i++) nrj+=x[i]*b[i];

    for(i=0;i<N;i++) r[i]=b[i]-Ax[i];

    for(i=0;i<N;i++){
      Ar[i]=0.0;
      for(j=0;j<N;j++) Ar[i]+=A[i][j]*r[j];
    }

    for(i=0;i<N;i++){
      Ap0[i]=0.0;
      for(j=0;j<N;j++) Ap0[i]+=A[i][j]*p0[j];
    }

    
    if(loop==0) {
      for(i=0;i<N;i++) p1[i]=r[i];
    } else {
      num_beta=0.0;
      for(i=0;i<N;i++) num_beta+=p0[i]*Ar[i];
      den_beta=0.0;
      for(i=0;i<N;i++) den_beta+=p0[i]*Ap0[i];
      beta=-num_beta/den_beta;
      for(i=0;i<N;i++) p1[i]=r[i]+beta*p0[i];
    }
    
    for(i=0;i<N;i++){
      Ap1[i]=0.0;
      for(j=0;j<N;j++) Ap1[i]+=A[i][j]*p1[j];
    }

    num_alpha=0.0;den_alpha=0.0;
    for(i=0;i<N;i++){
      num_alpha+=p1[i]*r[i];
      den_alpha+=p1[i]*Ap1[i];
    }

    if(loop>0) {
      dnrj=nrj_old-nrj;
      printf("loop %d -> nrj = %12.6f (%e)\n",loop,nrj,dnrj);
      out=fopen("debug.log","a+");      fprintf(out,"%d %g %e %e\n",loop,nrj,dnrj,num_alpha); fclose(out);
      if(fabs(dnrj)<tol) cvg=TRUE;
    }

    
    alpha=num_alpha/den_alpha;
    for(i=0;i<N;i++) x[i]+=alpha*p1[i];
    for(i=0;i<N;i++) p0[i]=p1[i];
    loop++;

    
  }


  if(cvg==TRUE)  printf("JOB DONE !\n");
  else {
    printf("!!!! WARNING !!! JOB not DONE !\n");
  }

  free(Ax);free(Ar);free(r);free(p0);free(p1);free(Ap0);free(Ap1);
  
}
