#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0

/* -------------------------------------------------------------------------------------------------------------- */
void diagonalization(double **A,int n,double *lambda,double **y);
double dot(int N,double *vec1,double *vec2);
double **GramSchmidt(int N,double **set0,int nvec0,double *set1,int *LI);
/* -------------------------------------------------------------------------------------------------------------- */
void davidson2D(int N,int M,double **v,double ax,double bx,double ay,double by,int nev,int first_ev){
  
  /* 
     OUTPUT: the eigenvectors are given in **v.
   */
  int NM;NM=N*M;
  int i,j,k,l;
  double tol=1.0e-4;
  double etol=1.0e-6;
  int nloopmax=5000;
  int cvg; 
  double norm;  
  double **T;			/* reduced matrix */
  double **y;
  double *lambda;
  double **Ritz;
  double **r;          r=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) r[i]=malloc(NM*sizeof(double));
  double *normr;   normr=malloc(nev*sizeof(double));
  double *ev;        ev=malloc(nev*sizeof(double));
  double *dev;      dev=malloc(nev*sizeof(double));  for(i=0;i<nev;i++) dev[i]=2*etol;

  double **t;         t=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) t[i]=malloc(NM*sizeof(double));
  double **q;        q=malloc(nev*sizeof(double*)); for(i=0;i<nev;i++) q[i]=malloc(NM*sizeof(double));
  int *q_insert;q_insert=malloc(nev*sizeof(int));
  double *beta;  if(nev>1) beta=malloc((nev-1)*sizeof(double));
  double *alpha;

  int loop;loop=1;
  FILE *log;log=fopen("debug.log","w+");fclose(log);
  int nq,iq,first_q;
  int nvecprev;

  int IL;
  int endloop=FALSE;
  int nvecini=  nvecini=first_ev+nev;
  int nvecmax; nvecmax=40;
  int nvec; nvec=nvecini;
  double **V;V=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) V[i]=malloc(NM*sizeof(double));
  for(j=0;j<nvecini;j++) for(i=0;i<NM;i++) V[j][i]=v[j][i];


  /* ----------------------
     Davidson loop 
  ------------------------ */
  while(endloop==FALSE){
    printf("----------------------------------------\n");
    /* ----------------------------------------------------------------------------------
       Computation of the reduced matrix T(nvec,nvec)=VT.A.V to be diagonalized in the 
       reduced subspace
       ------------------------------------------------------------------------------------ */
    T=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) T[i]=malloc(nvec*sizeof(double));
    for(i=0;i<nvec;i++){
      for(j=0;j<nvec;j++){
	T[i][j]=
	  V[i][0]*(ax+ay)*V[j][0]+V[i][0]*bx*V[j][1]+V[i][0]*by*V[j][N]+
	  V[i][N-1]*(ax+ay)*V[j][N-1]+V[i][N-1]*bx*V[j][N-2]+V[i][N-1]*by*V[j][2*N-1]+
	  V[i][N*(M-1)]*(ax+ay)*V[j][N*(M-1)]+V[i][N*(M-1)]*bx*V[j][N*(M-1)+1]+V[i][N*(M-1)]*by*V[j][N*(M-1)-N]+
	  V[i][N*M-1]*(ax+ay)*V[j][N*M-1]+V[i][N*M-1]*bx*V[j][N*M-2]+V[i][N*M-1]*by*V[j][N*M-1-N];
	for(l=1;l<=N-2;l++)
	  T[i][j]+=V[i][l]*(ax+ay)*V[j][l]+V[i][l]*bx*(V[j][l-1]+V[j][l+1])+V[i][l]*by*V[j][l+N];
	for(l=N*(M-1)+1;l<=N*M-2;l++)
	  T[i][j]+=V[i][l]*(ax+ay)*V[j][l]+V[i][l]*bx*(V[j][l-1]+V[j][l+1])+V[i][l]*by*V[j][l-N];
	for(l=N;l<=N*(M-2);l+=N)
	  T[i][j]+=V[i][l]*(ax+ay)*V[j][l]+V[i][l]*bx*V[j][l+1]+V[i][l]*by*(V[j][l-N]+V[j][l+N]);
	for(l=2*N-1;l<=N*(M-1)-1;l+=N)
	  T[i][j]+=V[i][l]*(ax+ay)*V[j][l]+V[i][l]*bx*V[j][l-1]+V[i][l]*by*(V[j][l-N]+V[j][l+N]);
	for(k=0;k<M-2;k++)
	  for(l=N+1+k*N;l<=2*N-2+k*N;l++)
	    T[i][j]+=V[i][l]*(ax+ay)*V[j][l]+V[i][l]*bx*(V[j][l-1]+V[j][l+1])+V[i][l]*by*(V[j][l-N]+V[j][l+N]);
      }
    }
    /* -------------------------------
       diagonalization of T(nvec,nvec)
       ------------------------------ */
    y=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) y[i]=malloc(nvec*sizeof(double));
    lambda=malloc(nvec*sizeof(double));
    diagonalization(T,nvec,lambda,y);
    for(i=0;i<nvec;i++) free(T[i]);free(T);
    for(i=0;i<nvec;i++) printf("%e ",lambda[i]);printf("\n");
    /* -----------------------------------------------------------------------
       expression of the |Ritz vectors> = A |y>
       nev is the number of eigenvalues to compute : nev=last_ev-first_ev+1
       There are nvec+nev Ritz vectors : 
           - nvec coming from the diagonalization of T; they are computed here
           - nev 
       ----------------------------------------------------------------------- */
    Ritz=malloc(nvec*sizeof(double*));  for(i=0;i<nvec;i++) Ritz[i]=malloc(NM*sizeof(double));
    for(k=0;k<nvec;k++){
      for(l=0;l<NM;l++) Ritz[k][l]=0.0;
      for(i=0;i<nvec;i++) 	for(l=0;l<NM;l++) Ritz[k][l]+=V[i][l]*y[k][i];
    }
    for(i=0;i<nvec;i++) free(y[i]);free(y);
    /* ------------------------------------------------------------------
       computation of the residuals |r_k>=A |Ritz_k> - lambda_k |Ritz_k>
       only yhr norm of the eigenvalues we search are computed
       ------------------------------------------------------------------ */
    for(k=0;k<nev;k++){
      r[k][0]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][0]+bx*Ritz[k+first_ev][1]+by*Ritz[k+first_ev][N];
      r[k][N-1]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][N-1]+bx*Ritz[k+first_ev][N-2]+by*Ritz[k+first_ev][2*N-1];
      r[k][N*(M-1)]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][N*(M-1)]+bx*Ritz[k+first_ev][N*(M-1)+1]+by*Ritz[k+first_ev][N*(M-1)-N];
      r[k][N*M-1]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][N*M-1]+bx*Ritz[k+first_ev][N*M-2]+by*Ritz[k+first_ev][N*M-1-N];
      for(l=1;l<=N-2;l++)
	r[k][l]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][l]+bx*(Ritz[k+first_ev][l-1]+Ritz[k+first_ev][l+1])+by*Ritz[k+first_ev][l+N];
      for(l=N*(M-1)+1;l<=N*M-2;l++)
	r[k][l]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][l]+bx*(Ritz[k+first_ev][l-1]+Ritz[k+first_ev][l+1])+by*Ritz[k+first_ev][l-N];
      for(l=N;l<=N*(M-2);l+=N)
	  r[k][l]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][l]+bx*Ritz[k+first_ev][l+1]+by*(Ritz[k+first_ev][l-N]+Ritz[k+first_ev][l+N]);	  
      for(l=2*N-1;l<=N*(M-1)-1;l+=N)
	  r[k][l]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][l]+bx*Ritz[k+first_ev][l-1]+by*(Ritz[k+first_ev][l-N]+Ritz[k+first_ev][l+N]);	  

      for(i=0;i<M-2;i++) {
	for(l=N+1+i*N;l<=2*N-2+i*N;l++){
	  r[k][l]=(ax+ay-lambda[k+first_ev])*Ritz[k+first_ev][l]+bx*(Ritz[k+first_ev][l-1]+Ritz[k+first_ev][l+1])+by*(Ritz[k+first_ev][l-N]+Ritz[k+first_ev][l+N]);
	}
      }
      normr[k]=0.0;
      for(l=0;l<NM;l++)     normr[k]+=r[k][l]*r[k][l];
      normr[k]=pow(normr[k],.5);
    }

    for(k=0;k<nev;k++) {
      if(loop>1) dev[k]=ev[k]-lambda[k+first_ev];
      ev[k]=lambda[k+first_ev];
    }
    
    /* -------------------------------------------------------------------------
     print the logfile
    --------------------------------------------------------------------------- */
    printf("loop %d nvecmax= %d\n",loop,nvecmax);
    if(loop>1){
      printf("dev= ");for(k=0;k<nev;k++)  printf("%e ",dev[k]);printf("\n");
    }
    printf("normr= ");for(k=0;k<nev;k++)  printf("%e ",normr[k]);printf("\n");
    printf("lambda= ");for(k=0;k<nev;k++) printf("%e ",lambda[k+first_ev]);printf("\n");
    log=fopen("debug.log","a+");
    fprintf(log,"loop= %d lambda= ",loop);
    for(k=0;k<nev;k++) fprintf(log,"%e ",lambda[k+first_ev]);//fprintf(log,"\n");
    fprintf(log,"normr= ");for(k=0;k<nev;k++)  fprintf(log,"%e ",normr[k]);fprintf(log,"\n");
    fclose(log);    
    /* -------------------------------------------------------------------------
       check the convergence
    --------------------------------------------------------------------------- */
   cvg=TRUE;

    if(loop>1){
      i=0;
      while(i<nev && cvg==TRUE){
	if(normr[i]>tol) cvg=FALSE;
	//if(fabs(dev[i])>etol) cvg=FALSE;
	i++;
      }
    } else {
      cvg=FALSE;
    }

    if(cvg==FALSE){
      if(nvec+nev<=nvecmax){	/* one computes nev additional vectors to expand the Krylov's subspace */
	/* --------------------------------------------------------------------------------------
	   computation of |t_k>=M|r_k>
	   M is the preconditioning matrix. It is an approximation of (A-lambdaI)^{-1}
	   Here M=(D-lambda I)^{-1}  where D is the main diagonale of A -> Jacobi Preconditionner
	   ------------------------------------------------------------------------------------ */

	for(i=0;i<nvec;i++) free(V[i]);free(V);
	for(k=0;k<nev;k++){
	  for(i=0;i<NM;i++) 	t[k][i]=r[k][i]/(ax+ay-lambda[k+first_ev]);
	  V=GramSchmidt(NM,Ritz,nvec,t[k],&IL);
	  if(IL==TRUE){
	    for(j=0;j<nvec;j++)      free(Ritz[j]);free(Ritz); nvec++;   Ritz=malloc(nvec*sizeof(double*));
	    for(j=0;j<nvec;j++){
	      Ritz[j]=malloc(NM*sizeof(double));
	      for(l=0;l<NM;l++) Ritz[j][l]=V[j][l];
	    }
	  } else {
	    printf("k=%d/%d IL= %d\n",k,nvec,IL);
	  }
	}
	//if(loop==240) exit(0);
	/* /\* PRINT BLOCK *\/ */
	log=fopen("debug.log","a+");
	fprintf(log,"===================================\n");
	fprintf(log,"# nvec = %d\n",nvec);
	/* for(l=0;l<NM;l++) { */
	/*   for(i=0;i<nvec;i++) printf("%12.6f ",V[i][l]); */
	/*   printf("\n"); */
	/* } */
	/* printf("\n"); */
	for(i=0;i<nvec;i++){
	  for(j=0;j<nvec;j++){
	    norm=0.0;
	    for(l=0;l<NM;l++) norm+=V[i][l]*V[j][l];
	    fprintf(log,"%12.6f ",norm);
	  }
	  fprintf(log,"\n");
	}
	fprintf(log,"===================================\n");
	fclose(log);
	/* exit(0); */
	/* /\* END OF PRINT BLOCK *\/ */
	


      } else { // nvec+nev>nvecmax
	nvecprev=nvec;
	nvec=nvecini;
	for(i=0;i<nvecprev;i++) free(V[i]);free(V);
	V=malloc(nvec*sizeof(double*)); for(i=0;i<nvec;i++) V[i]=malloc(NM*sizeof(double));
	for(k=0;k<nvec;k++) {
	  for(l=0;l<NM;l++) V[k][l]=Ritz[k][l];
	}
	for(i=0;i<nvecprev;i++) free(Ritz[i]);free(Ritz);
      }
    } /* if(cvg==FALSE) */ 
    loop++;
    free(lambda);
    endloop=FALSE;
    if(loop>=nloopmax || cvg==TRUE) {
      endloop=TRUE;
      for(k=0;k<nev;k++){
	for(l=0;l<NM;l++) {
	  v[k][l]=Ritz[k+first_ev][l];
	}
      }
    }
    /* end of Davidson loop */
  } 
  
}
