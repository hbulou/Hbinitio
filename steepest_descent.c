#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define FALSE 0
#define TRUE 1
void SteepestDescent(double **A,double *b,int N,double *x){
  int i,j,loop;
  double tol=1.0e-6;
  int loopmax=10000;
  FILE *out; out=fopen("debug.log","w+");fclose(out);
  loop=0;
  int cvg=FALSE;
  double nrj,dnrj,nrj_old;
  double num_alpha,den_alpha,alpha;
  double *r;r=malloc(N*sizeof(double));
  double *Ax;Ax=malloc(N*sizeof(double));
  double *Ar;Ar=malloc(N*sizeof(double));

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
    
    num_alpha=0.0;den_alpha=0.0;
    for(i=0;i<N;i++){
      num_alpha+=r[i]*r[i];
      den_alpha+=r[i]*Ar[i];
    }
    
    if(loop>0) {
      dnrj=nrj_old-nrj;
      printf("loop %d -> nrj = %12.6f (%e)\n",loop,nrj,dnrj);
      out=fopen("debug.log","a+");      fprintf(out,"%d %g %e %e\n",loop,nrj,dnrj,num_alpha); fclose(out);
      if(fabs(dnrj)<tol) cvg=TRUE;
    }


    alpha=num_alpha/den_alpha;
    for(i=0;i<N;i++) x[i]+=alpha*r[i];

    loop++;
    //    for(i=0;i<N;i++) printf("%12.6f ",x[i]);printf(" | <r|r>= %e\n",num_alpha);
  }




  if(cvg==TRUE)  printf("JOB DONE !\n");
  else {
    printf("!!!! WARNING !!! JOB not DONE !\n");
  }



}
